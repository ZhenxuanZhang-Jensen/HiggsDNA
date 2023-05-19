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
from higgs_dna.utils.misc_utils import load_config, update_dict, is_json_serializable
from higgs_dna.constants import NOMINAL_TAG, CENTRAL_WEIGHT, BRANCHES
from higgs_dna.utils.metis_utils import do_cmd


def run_analysis(config):
    """
    Function that gets run for each individual job. Performs the following:
        1. Load events from input files
        2. Add any relevant ``Sample`` metadata to events
        3. Produce all ``WeightSystematic``s and ``SystematicWithIndependentCollection``s that **do not** rely on selections performed by or fields added by ``Tagger``s.
        4. Apply ``TagSequence``
        5. Produce all remaining ``WeightSystematic``s and ``SystematicWithIndependentCollection``s that **do** rely on selections performed by or fields added by ``Tagger``s.
        6. Write output events 

    In general, you should not have to call this function yourself.
    It will be configured by the ``AnalysisManager`` for each job.

    :param config: dictionary/json that specifies the physics content (samples, tag sequence, systematics) of a given job
    :type config: dict or str
    """
    t_start = time.time()
    config = load_config(config)
    job_summary = { 
            "config" : config
    }
    ### 1. Load events ###
    t_start_load = time.time()
    events, sum_weights = AnalysisManager.load_events(config["files"], config["branches"])

    # Record n_events and sum_weights for scale1fb calculation
    job_summary["n_events"] = len(events)
    job_summary["sum_weights"] = sum_weights
    t_elapsed_load = time.time() - t_start_load

    ### 2. Add relevant sample metadata to events ###
    t_start_samples = time.time()
    sample = Sample(**config["sample"])
    events = sample.prep(events)
    t_elapsed_samples = time.time() - t_start_samples

    ### 3. Produce systematics ###
    t_start_syst = time.time()
    systematics_producer = SystematicsProducer(
        name = config["name"],
        options = config["systematics"],
        sample = sample
    )
    events = systematics_producer.produce(events)
    t_elapsed_syst = time.time() - t_start_syst

    ### 4. Apply tag sequence ###
    t_start_taggers = time.time()
    tag_sequence = TagSequence(
        name = config["name"],
        tag_list = config["tag_sequence"],
        sample = sample,
        output_dir=config["output_dir"]
    )
    events, tag_idx_map = tag_sequence.run(events)
    t_elapsed_taggers = time.time() - t_start_taggers

    ### 5. Compute remaining systematics that rely on tagger outputs ###
    t_start_apply_syst = time.time()
    events = systematics_producer.apply_remaining_weight_systs(events, tag_idx_map)
    t_elapsed_apply_syst = time.time() - t_start_apply_syst
    t_elapsed_syst += t_elapsed_apply_syst

    # Record number of selected events for each SystematicWithIndependentCollection
    job_summary["n_events_selected"] = {}
    for syst_tag, syst_events in events.items():
        job_summary["n_events_selected"][syst_tag] = len(syst_events)

    ### 6. Write selected events ###
    # Sometimes this may be running on a remote node that does not have access to the full host filesystem.
    # So, check if the relevant directories exist to save outputs in their full path and save them to current dir if not.
    if not os.path.exists(config["dir"]):
        output_dir = ""
        config["summary_file"] = os.path.split(config["summary_file"])[-1] 
    else: 
        output_dir = os.path.abspath(config["dir"]) + "/"
    output_name = output_dir + config["output_name"]
    
    job_summary["outputs"] = AnalysisManager.write_events(events, config["variables_of_interest"], output_name) 
    t_elapsed = time.time() - t_start

    # Calculate performance metrics
    job_summary["time"] = t_elapsed
    job_summary["time_frac_load"] = t_elapsed_load / t_elapsed
    job_summary["time_frac_samples"] = t_elapsed_samples / t_elapsed
    job_summary["time_frac_syst"] = t_elapsed_syst / t_elapsed
    job_summary["time_frac_taggers"] = t_elapsed_taggers / t_elapsed
    job_summary["successful"] = True

    # Dump json summary
    with open(config["summary_file"], "w") as f_out:
        json.dump(job_summary, f_out, sort_keys = True, indent = 4)
    return job_summary


class AnalysisManager():
    """
    Manages the running of an entire analysis.
    See HiggsDNA tutorial for more details:
        - https://sam-may.github.io/higgs_dna_tutorial.github.io/

    The ``AnalysisManager`` class is designed such that you can launch an instance through ``scripts/run_analysis.py``, kill the script, and then relaunch to resume running your analysis.
    ``AnalysisManager`` will periodically pickle itself to enable saving/loading its progress.
    This is safe to do and you don't need to worry about jobs being forgotten about and/or being submitted multiple times.
    Most attributes are not allowed to change when loading to prevent things like inconsistent mappings between input files and jobs between runnings. 
    
    If you are running an entirely new analysis, you must specify a ``config`` json or dict with the fields listed in ``REQUIRED_FIELDS``.

    If you are resuming a previously paused/killed instance, it is enough to just pass the ``output_dir`` of your previous run.
    When resuming a previously paused/killed instance (i.e. same ``output_dir``), it is possible to change the attributes listed in ``OVERWRITABLE_FIELDS``.
    All other attributes will be taken from the pickled instance.

    :param output_dir: path to output directory
    :type output_dir: str
    :param config: json/dictionary that specifies the physics content (samples, tag sequence, systematics) of an entire analysis
    :type config: dict or str, optional if and only if you specify an ``output_dir`` that has a pickled ``AnalysisManager`` instance from a previous run
    """
    REQUIRED_FIELDS = ["variables_of_interest", "tag_sequence", "systematics", "samples", "branches"] # these fields **must** be present in config file
    OVERWRITABLE_FIELDS = ["merge_outputs", "unretire_jobs", "retire_jobs", "reconfigure_jobs", "n_cores", "log-level", "log-file", "short"] # these fields can be safely changed between different runs
    # TODO: add functionality for switching between batch_system = "local" and batch_system = "condor" between runs
    DEFAULTS = { # kwargs with default values
            "name" : "my_analysis", # doesn't affect actual code, just for more informative printouts
            "function" : {
                "module_name" : "higgs_dna.analysis",
                "function_name" : "run_analysis"
            },
            "batch_system" : "local",
            "fpo" : None, # number of input files per output file (i.e. per job)
            "n_cores" : 4, # number of cores for local running
            "use_xrdcp" : False, # xrdcp to local or not
            "merge_outputs" : False,
            "unretire_jobs" : False,
            "retire_jobs" : False,
            "reconfigure_jobs" : False,
            "short" : False # run only 1 job per sample x year
    }

    def __init__(self, output_dir = "output", config = {}, **kwargs):
        self.output_dir = os.path.abspath(output_dir)
        self.pickle_file = self.output_dir + "/analysis_manager.pkl"
        self.summary_file = self.output_dir + "/summary.json"

        # Resuming a previous run?
        #   Check if there is a pickled AnalysisManager instance. If so, you only need to properly specify ``output_dir`` and everything else will be loaded from the pkl.
        if os.path.exists(self.pickle_file): # resuming a previous run
            logger.warning("[AnalysisManager : __init__] We are loading a saved pickle file '%s' -- please be sure this behavior is inteneded. If you want to run a completely new analysis, please specify a new ``output_dir``." % (self.pickle_file))
            if config:
                logger.warning("[AnalysisManager : __init__] A non-empty ``config`` arg was passed. This will be ignored and the ``AnalysisManager`` state will be loaded from the pkl file. If you are passing the same ``config`` that was used previously, it is safe to ignore this warning.")
            
            # Load pickled ``AnalysisManager`` instance
            with open(self.pickle_file, "rb") as f_in:
                saved_analysis_manager = dill.load(f_in)
            for attr, value in saved_analysis_manager.__dict__.items():
                logger.debug("[AnalysisManager : __init__] Setting attribute '%s' as : %s" % (attr, str(value)))
                setattr(self, attr, value)
            
            # Update overwritable kwargs and warn about any others
            for k,v in kwargs.items():
                if k in self.OVERWRITABLE_FIELDS:
                    # Check if this is already present in pickled instance, and if so, let the user know if it is being changed from its previous value.
                    if hasattr(self, k):
                        if not v == getattr(self, k):
                            logger.info("[AnalysisManager : __init__] Attribute '%s' is being updated from its previous value of '%s' -> '%s'." % (k, getattr(self, k), v))
                    logger.debug("[AnalysisManager : __init__] Setting attribute '%s' as : '%s'." % (k, v)) 
                    setattr(self, k, v)
                elif hasattr(self, k):
                    logger.warning("[AnalysisManager : __init__] kwarg '%s' with value '%s' was given to constructor, but will be ignored in favor of the pickled value '%s'." % (k, v, getattr(self, k)))
                else:
                    logger.warning("[AnalysisManager : __init__] Not sure what to do with kwarg '%s' with value '%s'." % (k, v))

            # Modify ``JobsManager`` and ``Task`` objects based on any updated values for ``OVERWRITABLE_FIELDS``
            self.modify_jobs()


        # Running a new analysis
        else:
            self.config = copy.deepcopy(load_config(config))
            logger.info("[AnalysisManager : __init__] Initializing a new AnalysisManager instance, as no previous pkl file was found to load state from.")

            # Check for required fields in config
            if any([x not in config.keys() for x in self.REQUIRED_FIELDS]):
                logger.exception("[AnalysisManager : __init__] The fields '%s' are required to be present in the ``config`` json/dict, but one or more were not found. The found keys were : '%s'." % (str(self.REQUIRED_FIELDS), str(config.keys())))
                raise ValueError()
            # Set config and kwargs as attributes
            for dicts in [config, kwargs]:
                for k,v in dicts.items():
                    logger.debug("[AnalysisManager : __init__] Setting attribute '%s' as : '%s'." % (k, v))
                    setattr(self, k, v)
            # Set any remaining attributes from default values
            for k,v in self.DEFAULTS.items():
                if not hasattr(self, k):
                    setattr(self, k, v)

            self.update_samples()

            # make output dir
            os.system("mkdir -p %s" % (self.output_dir))

            if self.variables_of_interest is None:
                self.variables_of_interest = []
            
            
            # Create jobs manager
            if self.batch_system.lower() in ["condor", "htcondor"]:
                self.jobs_manager = CondorManager(output_dir = self.output_dir)
            elif self.batch_system.lower() == "local":
                self.jobs_manager = LocalManager(output_dir = self.output_dir, n_cores = self.n_cores)

            # Create samples manager
            self.sample_manager = SampleManager(**self.samples)

            # Test construction of tag sequence and systematics producer
            self.test_construction()

            self.prepared_analysis = False


    def update_samples(self):
        """
        Propagate correct behavior of --sample_list and --years args from command line.
        """
        # Check for command line updates to samples
        if hasattr(self, "years"):
            years = self.years.split(",")
            if not years == self.samples["years"]:
                logger.warning("[AnalysisManager : update_samples] Years were provided through the config as '%s', but were also specified from the command line as '%s', which is the version we will use." % (str(self.samples["years"]), str(years)))
                self.samples["years"] = years

        if hasattr(self, "sample_list"):
            samples = self.sample_list.split(",")
            if not samples == self.samples["sample_list"]:
                logger.warning("[AnalysisManager : update_samples] Sample list was provided through the config as '%s', but was also specified from the command line as '%s', which is the version we will use." % (str(self.samples["sample_list"]), str(samples)))
                self.samples["sample_list"] = samples


    def modify_jobs(self):
        """
        Modify the behavior of ``JobsManager`` and ``Task`` instances based on new values given through kwargs.
        """
        # Check if we switched batch system
        if self.batch_system == "local" and not isinstance(self.jobs_manager, LocalManager):
            logger.info("[AnalysisManager : modify_jobs] Converting all unfinished jobs from CondorJob -> LocalJob.")
            self.jobs_manager = self.jobs_manager.convert_to_local() # convert from CondorManager -> LocalManager
        elif self.batch_system.lower() in ["condor", "HTCondor"] and not isinstance(self.jobs_manager, CondorManager):
            logger.info("[AnalysisManager : modify_jobs] Converting all unfinished jobs from LocalJob -> CondorJob.")
            self.jobs_manager = self.jobs_manager.convert_to_condor() # convert from LocalManager -> CondorManager
            

        # Update n_cores for ``LocalManager``
        if isinstance(self.jobs_manager, LocalManager):
            self.jobs_manager.n_cores = self.n_cores

        # Check if we were previously running with --short, but was removed for this run
        if not self.short:
            prev_short = any([task.max_jobs >= 0 for task in self.jobs_manager.tasks])
            if prev_short:
                logger.info("[AnalysisManager : modify_jobs] It appears you previously ran with the ``--short`` option but have now removed it. Will submit the full set of jobs for each task.")
                for task in self.jobs_manager.tasks:
                    task.max_jobs = -1
                    task.create_jobs() # create the new jobs (old ones will not be overwritten)
                    self.jobs_manager.customized_jobs = False # need to configure the new jobs
                self.jobs_manager.remerge = self.merge_outputs

        # Reconfigure jobs
        if self.reconfigure_jobs:
            logger.info("[AnalysisManager : modify_jobs] Forcing reconfiguration of jobs (rewriting python configs, executables and condor_submit files).")
            self.jobs_manager.customized_jobs = False
            self.jobs_manager.submit_jobs(summarize = False, dry_run = True)

        # Retire jobs
        if self.retire_jobs:
            logger.warning("[AnalysisManager : modify_jobs] Retiring all unfinished jobs.")
            for task in self.jobs_manager.tasks:
                task.retire_jobs()

        # If ``merge_outputs`` and ``unretire_jobs`` selected, check if any jobs actually get unretired. If so, we need to remerge.
        if self.unretire_jobs:
            logger.info("[AnalysisManager : modify_jobs] Unretiring jobs that failed up to the maximum number of retries.")
            unretired_jobs = False
            for task in self.jobs_manager.tasks:
                task_unretired_jobs = task.unretire_jobs()
                unretired_jobs = unretired_jobs or task_unretired_jobs 

            self.jobs_manager.remerge = unretired_jobs and self.merge_outputs


    def test_construction(self):
        """
        Construct, but do not run, TagSequence and SystematicsProducer.
        """
        logger.info("[AnalysisManager: test_construction] Testing TagSequence and SystematicsProducer construction.")
        test_syst = SystematicsProducer(options = self.systematics)
        test_tags = TagSequence(tag_list = self.tag_sequence)
        logger.info("[AnalysisManager: test_construction] Constructed TagSeqeunce and SystematicsProducer successfully.")


    def run(self):
        """
        Prepare and run analysis until finished, then merge outputs if requested and print out some diagnostic info.
        The AnalysisManager frequently pickles itself throughout this function, so it can easily pick up right where it
        left off if the script is paused or killed.
        """
        logger.debug("[AnalysisManager : run] Running analysis '%s'." % self.name)

        # Measure runtime
        start = time.time()

        # Load samples
        self.samples = self.sample_manager.get_samples()
        self.save()

        if not self.prepared_analysis:
            self.prepare_analysis()
            self.save()

        logger.info("[AnalysisManager : run] Running %d tasks with %d total input files split over %d total jobs." % (len(self.jobs_manager.tasks), sum([len(x.files) for x in self.jobs_manager.tasks]), sum([len(x.jobs) for x in self.jobs_manager.tasks])))

        summary = self.jobs_manager.submit_jobs()
        while not self.jobs_manager.complete():
            self.save()
            summary = self.jobs_manager.submit_jobs()

        if self.merge_outputs:
            self.jobs_manager.merge_outputs(self.output_dir)
            self.save()

        self.jobs_manager.summarize()

        end = time.time() - start
        logger.info("[AnalysisManager : run] Finished running analysis '%s'. Elapsed time: %s (hours:minutes:seconds)." % (self.name, str(datetime.timedelta(seconds = end))))

        for task, info in summary.items():
            logger.debug("\n[AnalysisManager : run] Task '%s', PERFORMANCE summary" % (task))
            print('info["performance"]["time"]',info["performance"]["time"])
            if (info["performance"]["time"] != 0):
                logger.debug("\t [PERFORMANCE : %s] Processed %d total events in %s (hours:minutes:seconds) of total runtime (%.2f Hz)." % (task, info["physics"]["n_events_initial"], str(datetime.timedelta(seconds = info["performance"]["time"])), float(info["physics"]["n_events_initial"]) / info["performance"]["time"] ))
            for portion in ["load", "syst", "taggers"]:
                if not info["performance"]["time"] > 0:
                    continue
                logger.debug("\t [PERFORMANCE : %s] Fraction of runtime spent on %s: %.2f percent" % (task, portion, (100. * info["performance"]["time_%s" % portion]) / info["performance"]["time"]))

            logger.debug("[AnalysisManager : run] Task '%s', PHYSICS summary" % (task)) 
            for syst_tag, n_events in info["physics"]["n_events_selected"].items():
                if not info["physics"]["n_events_initial"] > 0:
                    continue
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


    def prepare_analysis(self):
        for sample in self.samples:
            task_branches = self.branches
            # Add extra branches which are always loaded, even if they are not explicitly specified in config
            if sample.is_data:
                task_branches += BRANCHES["data"][sample.year] + BRANCHES["data"]["any"]
            else:
                task_branches += BRANCHES["mc"][sample.year] + BRANCHES["mc"]["any"]

            # Make Sample instance json-able so we can save it in job config files
            jsonable_sample = copy.deepcopy(sample)
            jsonable_sample.files = [x.__dict__ for x in jsonable_sample.files]
            jsonable_sample = jsonable_sample.__dict__
            config = {
                "sample" : copy.deepcopy(jsonable_sample),
                "branches" : task_branches
            }
            for x in ["systematics", "tag_sequence", "function", "variables_of_interest"]:
                config[x] = copy.deepcopy(getattr(self, x))
                    
            self.jobs_manager.add_task(
                    Task(
                        name = sample.name,
                        output_dir = self.output_dir + "/" + sample.name,
                        files = sample.files,
                        config = config,
                        fpo = self.fpo if self.fpo is not None else sample.fpo, 
                        max_jobs = 1 if self.short else -1
                    )
            )
        self.prepared_analysis = True


    def save(self):
        with open(self.pickle_file.replace(".pkl", "_temp.pkl"), "wb") as f_out:
            dill.dump(self, f_out)

        os.system("mv %s %s" % (self.pickle_file.replace(".pkl", "_temp.pkl"), self.pickle_file))


    def summarize(self):
        """
        Save map of sample_name : process_id so different samples can be distinguished in merged outputs.
        TODO: aggregate diagnostic info from Taggers and Systematics 
        """
        self.summary = {
            "sample_id_map" : self.sample_manager.process_id_map,
            "config" : self.config
        }

        for k,v in vars(self).items():
            if k == "summary":
                continue
            if k in self.summary["config"].keys():
                continue # don't double save things
            if is_json_serializable(v):
                self.summary[k] = v

        with open(self.summary_file, "w") as f_out:
            json.dump(self.summary, f_out, sort_keys = True, indent = 4)


    @staticmethod
    def load_events(files, branches):
        """
        Load all branches in ``branches`` from "Events" tree from all nanoAODs in ``files`` into a single zipped ``awkward.Array``.
        Also calculates and returns the sum of weights from nanoAOD "Runs" tree.        

        :param files: list of files
        :type files: list of str
        :param branches: list of branches. If any branches are requested that are not present in the input nanoAOD, they will be silently omitted.
        :type branches: list of str or tuple
        :returns: array of events, sum of weights
        :rtype: awkward.Array, float
        """
        events = []
        sum_weights = 0
        use_xrdcp = False
        for file in files:
            if use_xrdcp:
                local_file_name = file.split("/")[-1]
                # local_file_name = file.replace("/","_")
                time.sleep(10)
                logger.debug("sleep 10 secs before xrdcp")
                os.system("xrdcp %s %s" % (file, local_file_name))
                file = local_file_name
                logger.debug("local file name: %s" %file )
                logger.debug("cp file to local")
            with uproot.open(file, timeout = 1800) as f:
                #attention new block to read lumi and save in the txt file to read whold data lumi to make sure we have enough lumi
                lumi = f['LuminosityBlocks'].arrays(['run','luminosityBlock'])
                Runs = f['Runs'].arrays(['run'])
                import itertools
                for i in range(len(Runs.run)):
                    list_lumi = lumi.luminosityBlock[lumi.run == Runs.run[i]].to_list()
                    range_list = [[t[0][1], t[-1][1]] for t in (tuple(g[1]) for g in itertools.groupby(enumerate(list_lumi), lambda list_lumi: list_lumi[1]-list_lumi[0]))]
                    # with open("/eos/user/z/zhenxuan/brilws/lumi_cal.txt","a") as ftxt:
                    #     ftxt.write("\n")
                    #     ftxt.write('"'+ str(Runs.run[i]) + '"')
                    #     ftxt.write(":")
                    #     ftxt.write(str(range_list))
                    #     ftxt.write(",")
                ##################################
                runs = f["Runs"]
                if "genEventCount" in runs.keys() and "genEventSumw" in runs.keys():
                    sum_weights += numpy.sum(runs["genEventSumw"].array())
                elif "genEventCount_" in runs.keys() and "genEventSumw_" in runs.keys():
                    sum_weights += numpy.sum(runs["genEventSumw_"].array())
                tree = f["Events"]
                trimmed_branches = [x for x in branches if x in tree.keys()]
                events_file = tree.arrays(trimmed_branches, library = "ak", how = "zip")
                events.append(events_file)

                logger.debug("[AnalysisManager : load_events] Loaded %d events from file '%s'." % (len(events_file), file))
            if use_xrdcp:
                os.system("rm %s" % file)
                logger.debug("remove the local cp file: %s" %file)



        events = awkward.concatenate(events)
        return events, sum_weights


    @staticmethod
    def write_events(events, save_branches, name):
        """
        For each set of events in ``events``, saves all fields in ``save_branches`` to a .parquet file.

        :param events: dictionary of systematic variations with independent collections : array of events
        :type events: dict
        :param save_branches: list of fields to save in output file
        :type save_branches: list of str or tuple or list
        :param name: name template for output files, which will be updated based on each key of ``events``
        :type name: str
        :returns: dictionary of keys from ``events`` : parquet file
        :rtype: dict
        """
        outputs = {}

        for syst_tag, syst_events in events.items():
            save_map = {}
            for branch in save_branches:
                if isinstance(branch, tuple) or isinstance(branch, list):
                    save_name = "_".join(branch)
                    if isinstance(branch, list):
                        branch = tuple(branch)
                else:
                    save_name = branch
                if isinstance(branch, tuple):
                    present = branch[1] in syst_events[branch[0]].fields
                else:
                    present = branch in syst_events.fields
                if not present:
                    logger.warning("[AnalysisManager : write_events] Branch '%s' was not found in events array. This may be expected (e.g. gen info for a data file), but please ensure this makes sense to you." % str(branch))
                    continue
                save_map[save_name] = syst_events[branch]

            for field in syst_events.fields:
                if "weight_" in field and not field in save_map.keys():
                    save_map[field] = syst_events[field]

            syst_events = awkward.zip(save_map,depth_limit=1) #attention : change depth_limit to 1
            out_name = "%s_%s.parquet" % (name, syst_tag)

            logger.debug("[AnalysisManager : write_events] Writing output file '%s'." % (out_name))
            awkward.to_parquet(syst_events, out_name) 
            outputs[syst_tag] = out_name

        return outputs
