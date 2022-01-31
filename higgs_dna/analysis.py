import os
import time
import datetime
import copy
import json
import pickle
import dill
import numpy

import uproot
import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.samples.sample_manager import SampleManager
from higgs_dna.samples.sample import Sample
from higgs_dna.job_management.managers import LocalManager, CondorManager
from higgs_dna.job_management.task import Task
from higgs_dna.systematics.systematics_producer import SystematicsProducer
from higgs_dna.taggers.tag_sequence import TagSequence
from higgs_dna.utils.misc_utils import load_config, update_dict
from higgs_dna.constants import NOMINAL_TAG, CENTRAL_WEIGHT, TRIGGER
from higgs_dna.utils.metis_utils import do_cmd

def run_analysis(config, summary = {}):
    """
    Function to be run inside each job.
    """
    t_start = time.time()

    config = load_config(config)
    job_summary = { 
            "config" : config
    }

    # Load events
    t_start_load = time.time()

    events, sum_weights = AnalysisManager.load_events_TEMP(config["files"], config["branches"])
    job_summary["n_events"] = len(events)
    job_summary["sum_weights"] = sum_weights

    t_elapsed_load = time.time() - t_start_load

    # 1. Add relevant sample metadata to events
    t_start_samples = time.time()

    sample = Sample(**config["sample"])
    events = sample.prep(events)

    t_elapsed_samples = time.time() - t_start_samples

    # 2. Produce systematics
    t_start_syst = time.time()

    systematics_producer = SystematicsProducer(
        name = config["name"],
        options = config["systematics"],
        sample = sample
    )
    events = systematics_producer.produce(events)

    t_elapsed_syst = time.time() - t_start_syst

    # 3. Apply tag sequence
    t_start_taggers = time.time()

    tag_sequence = TagSequence(
        name = config["name"],
        tag_list = config["tag_sequence"],
        sample = sample
    )
    events, tag_idx_map = tag_sequence.run(events)

    t_elapsed_taggers = time.time() - t_start_taggers

    # 4. Compute remaining systematics that rely on tagger outputs
    t_start_apply_syst = time.time()

    events = systematics_producer.apply_remaining_weight_systs(events, tag_idx_map)

    t_elapsed_apply_syst = time.time() - t_start_apply_syst
    t_elapsed_syst += t_elapsed_apply_syst

    job_summary["n_events_selected"] = {}
    for syst_tag, syst_events in events.items():
        job_summary["n_events_selected"][syst_tag] = len(syst_events)

    # Write selected events
    os.system("mkdir -p %s" % config["output_dir"])
    if not config["remote_job"]: # a remote job that does not have write access to your home area (e.g. a T2 condor job, as opposed to an lxplus condor job)
        output_name = config["job_batch_output_dir"] + "/output_job_%d" % config["job_id"]
    else:
        output_name = "output_job_%d" % config["job_id"]
    job_summary["outputs"] = AnalysisManager.write_events_TEMP(events, config["save_branches"], output_name) 

    t_elapsed = time.time() - t_start

    job_summary["time"] = t_elapsed
    job_summary["time_frac_load"] = t_elapsed_load / t_elapsed
    job_summary["time_frac_samples"] = t_elapsed_samples / t_elapsed
    job_summary["time_frac_syst"] = t_elapsed_syst / t_elapsed
    job_summary["time_frac_taggers"] = t_elapsed_taggers / t_elapsed

    job_summary["successful"] = True
    summary["job_%d" % config["job_id"]] = job_summary

    if not config["remote_job"]:
        summary_file = config["summary_file"]
    else:
        summary_file = config["summary_file"].split("/")[-1]
    with open(summary_file, "w") as f_out:
        json.dump(job_summary, f_out, sort_keys = True, indent = 4)

    return job_summary


class AnalysisManager():
    """
    Manages the running of an entire analysis.
    """
    def __init__(self, config = {}, name = None, function = None, samples = None, tag_sequence = None, systematics = None, variables_of_interest = None, batch_system = "local", fpo = None, output_dir = "output", merge_outputs = False, resubmit_retired = False, n_cores = 6, **kwargs):
        # TODO: check that config dict has all required fields
        self.config = copy.deepcopy(load_config(config))

        host = do_cmd("hostname")
        if "t2.ucsd" in host:
            self.host = "UCSD"
        elif "lxplus" in host:
            self.host = "lxplus"

        self.user = do_cmd("whoami")

        self.output_dir = os.path.abspath(output_dir)

        if (batch_system == "condor" or batch_system == "HTCondor") and "batch_output_dir" not in kwargs and self.host in ["UCSD"]:
            self.batch_output_dir = "/hadoop/cms/store/user/%s/HiggsDNA/%s/" % (self.user, output_dir)
        else:
            self.batch_output_dir = os.path.abspath(kwargs.get("batch_output_dir", self.output_dir))

        self.merge_outputs = merge_outputs
        self.resubmit_retired = resubmit_retired

        # Check if a pkl file for this analysis manager is present and load previous progress if so
        self.summary_file = self.output_dir + "/summary.json"
        self.pickle_file = self.output_dir + "/analysis_manager.pkl"

        if os.path.exists(self.pickle_file):
            logger.info("[AnalysisManager : __init__] Found previous pickle file '%s' for this analysis, loading previous state and progress." % (self.pickle_file))
            logger.warning("[AnalysisManager : __init__] We are loading a saved pickle file '%s' -- please be sure this behavior is inteneded. If you want to run a completely new analysis, please specify a new output directory or run with the option TODO" % (self.pickle_file))
            with open(self.pickle_file, "rb") as f_in:
                saved_analysis_manager = dill.load(f_in)

            for attr, value in saved_analysis_manager.__dict__.items(): 
                logger.debug("[AnalysisManager : __init__] attribute '%s' : %s" % (attr, str(value)))
                setattr(self, attr, value) 

            # Overwrite the saved value for these with the values given from command line
            self.merge_outputs = merge_outputs
            self.resubmit_retired = resubmit_retired
            self.n_cores = n_cores
            
            if isinstance(self.jobs_manager, LocalManager):
                self.jobs_manager.n_cores = n_cores

            if self.resubmit_retired:
                logger.debug("[AnalysisManager : __init__] Unretiring jobs that failed up to the maximum number of retries.")
                remerge = False
                for task in self.jobs_manager.tasks:
                    unretired_jobs = task.unretire_jobs()
                    remerge = remerge or unretired_jobs

                if remerge: # if at any point we unretired jobs, we should remerge outputs
                    self.jobs_manager.remerge = True # however, we only want to update if we need to remerge
                    # This prevents the scenario of user ctrl+c-ing after running with --resubmit_retired, running again without that option, and us mistakenly thinking that we no longer need to remerge.

        # Otherwise, run normal constructor
        else:
            # First, check if the config passes these options
            fields = ["name", "function", "variables_of_interest", "tag_sequence", "systematics", "samples", "branches"]
            for field in fields:
                if field in self.config.keys():
                    setattr(self, field, self.config[field])

            # Second, check if they were provided manually.
            # If provided both through the config and explicitly, we take the explicitly provided argument.
            args = [name, function, variables_of_interest, tag_sequence, systematics, samples]
            for field, arg in zip(fields, args):
                if field == "samples":
                    for x in ["sample_list", "years"]:
                        sample_arg = kwargs.get(x, None)
                        if sample_arg is not None:
                            sample_arg = sample_arg.split(",")
                            logger.warning("[AnalysisManager : __init__] Samples argument '%s' was provided both through the config json and explicitly in the constructor. We will take the version provided by the constructor." % x)
                            self.samples[x] = sample_arg
                            self.config["samples"][x] = sample_arg


                elif arg is not None:
                    if hasattr(self, field):
                        logger.warning("[AnalysisManager : __init__] Argument '%s' was provided both through the config json and explicitly in constructor. We will take the version provided by the constructor." % field)
                    setattr(self, field, arg)
                    self.config[field] = arg # update internal config dict to properly document details for this analysis


            self.batch_system = batch_system
            self.fpo = fpo
            self.n_cores = n_cores

            for dir in [self.output_dir, self.batch_output_dir]:
                os.system("mkdir -p %s" % dir)

            if self.variables_of_interest is None:
                self.variables_of_interest = []

            self.save_branches = self.variables_of_interest + [CENTRAL_WEIGHT] #TODO: save other weight variations as well

            if self.name is None:
                self.name = "my_analysis"

            # Create jobs manager
            if self.batch_system == "local":
                self.jobs_manager = LocalManager(n_cores = self.n_cores)
            elif self.batch_system == "HTCondor" or self.batch_system == "condor":
                self.jobs_manager = CondorManager(output_dir = self.output_dir, batch_output_dir = self.batch_output_dir)

            # Create samples manager
            self.sample_manager = SampleManager(**self.config["samples"])

            logger.info("[AnalysisManager : __init__] Created AnalysisManager with following specifications:")
            for attr, val in self.__dict__.items():
                message = "\t '%s' : %s" % (attr, str(val))
                if attr == "config":
                    logger.debug(message)
                else:
                    logger.info(message)

            self.test_construction()
        
            self.prepared_analysis = False

    def test_construction(self):
        logger.info("[AnalysisManager: test_construction] Testing TagSequence and SystematicsProducer construction.")
        test_syst = SystematicsProducer(options = self.config["systematics"])
        test_tags = TagSequence(tag_list = self.config["tag_sequence"])
        logger.info("[AnalysisManager: test_construction] Constructed TagSeqeunce and SystematicsProducer successfully.")


    @staticmethod
    def create_chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        chunks = []
        for i in range(0, len(lst), n):
            chunks.append(lst[i:i+n])
        return chunks


    def run(self, max_jobs = None):
        """
        Create SystematicsProducer and TagSequence for each sample/year,
        split input files into jobs,
        submit & monitor jobs (through higgs_dna.job_management.JobsManager),
        and summarize.

        :param max_jobs: maximum number of total jobs to run, defaults to None. If None, will run all jobs. If an int > 0, will run at most that many jobs per sample.
        :type max_jobs: int
        """
        logger.debug("[AnalysisManager : run] Running analysis '%s'." % self.name)

        # Measure runtime
        start = time.time()

        # Load samples
        self.samples = self.sample_manager.get_samples()

        if not self.prepared_analysis:
            self.prepare_analysis(max_jobs)

        n_input_files = 0
        n_tasks = 0
        n_jobs = 0
        for task in self.jobs_manager.tasks:
            n_tasks += 1
            n_input_files += len(task.files)
            n_jobs += len(task.jobs)

        logger.info("[AnalysisManager : run] Running %d tasks with %d total input files split over %d total jobs." % (n_tasks, n_input_files, n_jobs))
        os.system("sleep 2s")

        summary = self.jobs_manager.submit_jobs(summarize = True)
        self.save()
        idx = 1
        while not self.jobs_manager.complete:
            os.system("sleep 1s")
            summarize = idx % 10 == 0
            summary = self.jobs_manager.submit_jobs(summarize = summarize)
            self.save()
            idx += 1


        if self.merge_outputs:
            self.jobs_manager.merge_outputs(self.output_dir)
            self.save()

        end = time.time() - start
        logger.info("[AnalysisManager : run] Finished running analysis '%s'. Elapsed time: %s (hours:minutes:seconds)." % (self.name, str(datetime.timedelta(seconds = end))))

        for task, info in summary.items():
            logger.debug("\n[AnalysisManager : run] Task '%s', PERFORMANCE summary" % (task))
            logger.debug("\t [PERFORMANCE : %s] Processed %d total events in %s (hours:minutes:seconds) of total runtime (%.2f Hz)." % (task, info["physics"]["n_events_initial"], str(datetime.timedelta(seconds = info["performance"]["time"])), float(info["physics"]["n_events_initial"]) / info["performance"]["time"] ))
            for portion in ["load", "syst", "taggers"]:
                logger.debug("\t [PERFORMANCE : %s] Fraction of runtime spent on %s: %.2f percent" % (task, portion, (100. * info["performance"]["time_%s" % portion]) / info["performance"]["time"]))

            logger.debug("[AnalysisManager : run] Task '%s', PHYSICS summary" % (task)) 
            for syst_tag, n_events in info["physics"]["n_events_selected"].items():
                logger.debug("\t [PHYSICS : %s] events set '%s' has eff. of %.2f percent (all taggers)" % (task, syst_tag, (float(n_events) * 100.) / float(info["physics"]["n_events_initial"])))
            if "scale1fb" in info["physics"].keys():
                logger.debug("\t [PHYSICS : %s] With a cross section times BF of %.9f pb and a sum of weights of %.9f, scale1fb for this sample is %.9f" % (task, info["physics"]["norm_factor"], info["physics"]["sum_weights"], info["physics"]["scale1fb"]))


        retired_jobs = []
        for task in self.jobs_manager.tasks:
            retired_jobs += [job for job in task.jobs if job.status == "retired"]

        if retired_jobs:
            logger.warning("[AnalysisManager : run] WARNING: there were %d retired jobs (meaning they failed up to the maximum number of retries). You may want to look into these:" % (len(retired_jobs)))
            for job in retired_jobs:
                logger.warning("\t %s" % job.name_full)

        self.summarize()


    def prepare_analysis(self, max_jobs):
        idx = 0
        for sample in self.samples:
            name = sample.name
            output_dir = self.output_dir + "/" + name + "/"
            batch_output_dir = self.batch_output_dir + "/" + name + "/"

            # 1. Branches to load/save
            if sample.is_data:
                task_branches = [x for x in self.branches if not ("gen" in x or "Gen" in x or "hadronFlavour" in x or "L1PreFiringWeight" in x or "Pileup_nTrueInt" in x)]
                task_save_branches = [x for x in self.save_branches if not ("weight" in x or "gen" in x or "Gen" in x)]
                for idx, x in enumerate(task_save_branches):
                    if isinstance(x, list):
                        for y in x:
                            if "weight" in y or "gen" in y or "Gen" in y or "hadronFlavour" in y:
                                task_save_branches.remove(x)
                task_branches += TRIGGER[sample.year]

            else:
                task_branches = [x for x in self.branches if not ("HLT" in x)]

                if sample.year == "2018":
                    task_branches = [x for x in task_branches if not ("L1PreFiringWeight" in x)]
                task_save_branches = self.save_branches

            

            # 2. Sample details
            task_sample = copy.deepcopy(sample)
            task_sample.files = [x.__dict__ for x in task_sample.files] # make compatible with json
            task_sample = task_sample.__dict__ # make compatible with json

            # 3. Systematics
            task_systematics = copy.deepcopy(self.config["systematics"])
            if sample.systematics is not None:
                task_systematics = update_dict(
                        original = task_systematics,
                        new = sample.systematics
                )

            # 4. Tag Sequence
            task_tag_sequence = copy.deepcopy(self.config["tag_sequence"])

            # 5. Wrapper function
            function = copy.deepcopy(self.config["function"])

            # 6. Job Details
            if self.fpo is not None:
                fpo = self.fpo
            else:
                if "Data" in name:
                    fpo = 10
                else:
                    fpo = 3
            task = Task(
                    name = name,
                    output_dir = output_dir,
                    batch_output_dir = batch_output_dir,
                    n_files_per_job = fpo, 
                    files = sample.files,
                    max_jobs = max_jobs,
                    config = {
                        "sample" : task_sample,
                        "systematics" : task_systematics,
                        "tag_sequence" : task_tag_sequence,
                        "function" : function,
                        "branches" : task_branches,
                        "save_branches" : task_save_branches
                    }
            )

            self.jobs_manager.add_task(task)

        self.prepared_analysis = True


    def save(self):
        with open(self.pickle_file.replace(".pkl", "_temp.pkl"), "wb") as f_out:
            dill.dump(self, f_out)

        os.system("mv %s %s" % (self.pickle_file.replace(".pkl", "_temp.pkl"), self.pickle_file))


    def summarize(self):
        """
        TODO
        """
        self.summary = {}
        self.summary["sample_id_map"] = self.sample_manager.process_id_map

        with open(self.summary_file, "w") as f_out:
            json.dump(self.summary, f_out, sort_keys = True, indent = 4)

        return


    def construct_analysis(self, sample):
        """
        Construct the SystematicsProducer and TagSequence objects specific to this sample/year

        :param sample: the sample (along with year) to create SystematicsProducer and TagSequence for
        :type sample: higgs_dna.samples.sample.Sample
        :return: SystematicsProducer and TagSequence objects specific to this sample/year
        :rtype: higgs_dna.systematics.systematics_producer.SystematicsProducer, higgs_dna.taggers.tag_sequence.TagSequence
        """
        # Create systematics producer specific to this sample (a specific sample/year may have additional systematics)
        sample_systematics_producer = SystematicsProducer(
                options = self.systematics_producer.options,
                name = sample.name,
                sample = sample
        )
        if sample.systematics is not None:
            sample_systematics_producer.add_systematics(sample.systematics)

        # Create tag sequence specific to this sample (a specific sample/year may have different selection options)
        #sample_tag_sequence = copy.deepcopy(self.tag_sequence)
        sample_tag_sequence = TagSequence(tag_list = self.tag_sequence)
        sample_tag_sequence.name = sample.name

        for attr in ["year", "process", "is_data"]:
            sample_tag_sequence.update_tagger_info(
                    attr = attr, 
                    value = getattr(sample, attr)
            )

        return sample_systematics_producer, sample_tag_sequence


    @staticmethod
    def load_events_TEMP(files, branches):
        """
        Temporary function: should be replaced by a dedicated class/function which optimally loads events from multiple files. 
        Maybe coffea features can help here? 
        """
        events = []
        sum_weights = 0
        for file in files:
            with uproot.open(file, timeout = 500) as f:
                runs = f["Runs"]
                if "genEventCount" in runs.keys() and "genEventSumw" in runs.keys():
                    sum_weights += numpy.sum(runs["genEventSumw"].array())
                elif "genEventCount_" in runs.keys() and "genEventSumw_" in runs.keys():
                    sum_weights += numpy.sum(runs["genEventSumw_"].array())
                tree = f["Events"]
                events.append(tree.arrays(branches, library = "ak", how = "zip"))

        events = awkward.concatenate(events)
        return events, sum_weights

    @staticmethod
    def write_events_TEMP(events, save_branches, name):
        """
        Temporary function: should be replaced by configurable way to write events in multiple output formats.
        Maybe coffea features can help here?
        """
        outputs = {}

        # Make any missing dirs/subdirs
        #dirs = "/".join(name.split("/")[:-1])
        #os.system("mkdir -p %s" % dir)

        for syst_tag, syst_events in events.items():
            save_map = {}
            for branch in save_branches:
                if isinstance(branch, tuple) or isinstance(branch, list):
                    save_name = "_".join(branch)
                    if isinstance(branch, list):
                        branch = tuple(branch)
                else:
                    save_name = branch
                save_map[save_name] = syst_events[branch]

            for field in syst_events.fields:
                if "weight_" in field and not field in save_map.keys():
                    save_map[field] = syst_events[field]

            syst_events = awkward.zip(save_map)
            out_name = "%s_%s.parquet" % (name, syst_tag)

            logger.debug("[AnalysisManager : write_events_TEMP] Writing output file '%s'." % (out_name))
            awkward.to_parquet(syst_events, out_name) 
            outputs[syst_tag] = out_name

        return outputs
