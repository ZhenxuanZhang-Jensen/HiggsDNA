import os
import importlib
import time
import awkward
import numpy
import json
import glob
import pyarrow
import datetime
import sys
from tqdm import tqdm

import logging
logger = logging.getLogger(__name__)

from higgs_dna.job_management.jobs import Job, LocalJob, CondorJob
from higgs_dna.utils.metis_utils import num_to_ordinal_string, do_cmd, do_cmd_timeout
from higgs_dna.utils import awkward_utils
from higgs_dna.utils.misc_utils import check_proxy, create_chunks, get_HiggsDNA_base, get_conda, expand_path
from higgs_dna.job_management.constants import CONDOR_STATUS_FLAGS, HOST_PARAMS, CONDA_TARFILE

import enlighten
en_manager = enlighten.get_manager() # for job progress status bar

class JobsManager():
    """

    """
    def __init__(self): 
        self.tasks = []
        self.prepared_inputs = False
        self.summary = {}
        self.job_type = Job
        self.jobs = []
        self.remerge = False
        self.customized_jobs = False
        self.copied_tars = False
        self.first_submit = True
        self.dashboard = []
        self.scroll_idx = 0

        # Determine host and user
        self.user = do_cmd("whoami")
        host = do_cmd("hostname")
        if "t2.ucsd" in host:
            self.host = "UCSD"
        elif "lxplus" in host:
            self.host = "lxplus"
        else:
            self.host = host

        self.conda_path = get_conda()
        self.higgs_dna_path = get_HiggsDNA_base() 

        self.starttime = time.time()


    def add_task(self, task):
        self.tasks.append(task)


    def get_n_jobs(self):
        if not self.tasks:
            return 0
        else:
            return sum([len(task.jobs) for task in self.tasks])


    def complete(self):
        """Check if all Tasks are complete."""
        return all([(task.complete and task.merged_output_files) for task in self.tasks])


    def submit_jobs(self, summarize = True, dry_run = False):
        """

        """
        if not self.prepared_inputs:
            self.prepare_inputs()

        if self.first_submit:
            logger.info("[JobsManager : submit_jobs] For the first round of submitting jobs, we need to configure the config and executable files for all jobs. This may take a while if you have a large number of jobs.")
            self.first_submit = False

        if not self.customized_jobs:
            self.customize_jobs()

        self.submit_to_batch(dry_run = dry_run)

        if summarize:
            self.summarize()

        return self.summary


    def submit_to_batch(self, dry_run = False):
        """

        """
        pass


    def get_runtime(self):
        runtime = time.time() - self.starttime
        h, m, s = str(datetime.timedelta(seconds = runtime)).split(":")
        result = "%.0f sec" % float(s)
        if float(m) > 0:
            result = ("%.0f min, " % float(m)) + result
        if float(h) > 0:
            result = ("%.0f hr, " % float(h)) + result
        return result


    def __getstate__(self):
        state = self.__dict__.copy()
        del state["dashboard"]
        return state


    def __setstate__(self, state):
        self.__dict__.update(state)
        self.dashboard = []


    def summarize(self):
        """

        """

        columns, lines = os.get_terminal_size()

        if len(self.tasks) > int(lines/2): # don't take up more than half the terminal window
            scrolling = True
        else:
            scrolling = False

        if not self.dashboard:
            for idx, task in enumerate(self.tasks):
                if idx >= int(lines/2):
                    continue
                self.dashboard.append(en_manager.status_bar(position = idx+1))

        if scrolling:
            task_sets = create_chunks(self.tasks, int(lines/2))
            if self.scroll_idx >= len(task_sets):
                self.scroll_idx = 0
            tasks = task_sets[self.scroll_idx]
            self.scroll_idx = (self.scroll_idx + 1) % len(task_sets)
        else:
            tasks = self.tasks

        for idx, task in enumerate(tasks):
            bar = task.pbar.bar
            if len(bar) > columns:
                bar = bar[:columns]
            self.dashboard[idx].update(bar)


    def merge_outputs(self, dir):
        """
        Merge merged_output files from all tasks into a single file for each systematic with independent collection.
        Input files which are missing fields with respect to other input files have those fields added with a dummy value of 1.
        A value of 1 is chosen since this should ordinarily only be the case for outputs from data which do not have fields for weights (or gen quantities).
        """

        self.outputs = {} # dictionary to order outputs by syst_tag (systematics with independent collections)

        for task in self.tasks:
            task.add_process_id()
            task.add_year()
            for syst_tag, merged_output in task.merged_outputs.items():
                if syst_tag not in self.outputs.keys():
                    self.outputs[syst_tag] = []
                self.outputs[syst_tag] += [merged_output]

        for syst_tag, outputs in self.outputs.items():
            merged_file = dir + "/merged_%s.parquet" % (syst_tag)
            if os.path.exists(merged_file) and not self.remerge:
                logger.debug("[JobsManager : merge_outputs] For syst_tag '%s' merged output file '%s' already exists, will not overwrite it." % (syst_tag, merged_file))
                continue
            elif os.path.exists(merged_file) and self.remerge:
                logger.debug("[JobsManager : merge_outputs] For syst_tag '%s' merged output file '%s' already exists, but some jobs were unretired or --short option was removed, so we will remerge inputs and overwrite this file." % (syst_tag, merged_file))

            merged_events = []
            for output in outputs:
                merged_events.append(awkward.from_parquet(output))
            
            # Add dummy values for missing fields (usually weight branches in data)
            branches = []
            for events in merged_events:
                branches += events.fields
            branches = set(branches)

            for idx, events in enumerate(merged_events):
                missing_branches = [b for b in branches if b not in events.fields]
                for b in missing_branches:
                    awkward_utils.add_field(
                            events = merged_events[idx],
                            name = b,
                            data = numpy.ones(len(merged_events[idx]))
                    )

            logger.info("[JobsManager : merge_outputs] For syst_tag '%s' merging %d files into merged file '%s'." % (syst_tag, len(outputs), merged_file))

            merged_events = awkward.concatenate(merged_events)
            
            logger.debug("\t merging %d events" % len(merged_events))
            awkward.to_parquet(merged_events, merged_file)


    def prepare_inputs(self):
        """

        """
        pass 


    def customize_jobs(self):
        """

        """
        pass


    def can_submit(self, job):
        """
        Whether a job can currently be submitted
        For local submission, this amounts to checking if the number of currently running jobs is less than the chosen number of cores.
        """
        return True

    
class LocalManager(JobsManager):
    """

    """
    def __init__(self, output_dir, **kwargs):
        super(LocalManager, self).__init__()

        self.n_cores = kwargs.get("n_cores", 6)
        self.output_dir = os.path.abspath(output_dir)
        self.n_running_jobs = 0
        self.job_type = LocalJob


    def customize_jobs(self):
        """

        """
        with tqdm(total = sum([len(task.jobs) for task in self.tasks])) as prog:
            for task in self.tasks:
                for job in task.jobs:
                    job.__class__ = self.job_type
                    job.output_dir = job.dir
                    job.set_files() # update the locations for the various files associated with this job
                    prog.update(1)
        self.customized_jobs = True


    def submit_to_batch(self, dry_run = False):
        """

        """
        for task in self.tasks:
            task.process()

            jobs_to_submit = [job for job in task.jobs if (job.status == "waiting" or job.status == "failed")]
            for job in jobs_to_submit:
                if self.can_submit(job):
                    submitted = job.submit(dry_run = dry_run)
                    if submitted and not dry_run:
                        logger.info("[LocalManager : submit_jobs] : %s submitted job '%s_%d' (for the %s time)" % (str(type(self)), job.name, job.idx, num_to_ordinal_string(job.n_attempts)))
                        os.system("sleep 0.1s")

            task.process()

            self.summary[task.name] = task.summary


    def can_submit(self, job):
        self.n_running_jobs = 0
        for task in self.tasks:
            for job in task.jobs:
                if job.status == "running":
                    self.n_running_jobs += 1
        return self.n_running_jobs < self.n_cores


    def convert_to_condor(self):
        """
        Convert itself into a CondorManager
        """
        raise NotImplementedError()
        """
        condor_manager = CondorManager(output_dir = self.output_dir)
        for k,v in self.__dict__.items():
            if k in ["job_type"]:
                continue
            setattr(condor_manager, k, v)

        condor_manager.customize_jobs()

        return condor_manager
        """

class CondorManager(JobsManager):
    """

    """
    def __init__(self, output_dir, **kwargs):
        super(CondorManager, self).__init__()

        if not self.host in HOST_PARAMS.keys():
            logger.warning("[CondorManager : __init__] Host site was identified as '%s', which is not in the list of currently supported host sites for HiggsDNA (supported sites : '%s'). We will default to `lxplus`-specific HTCondor options, but job submission with HTCondor will not be guaranteed to work! Local job submission should be unaffected." % (self.host, str(HOST_PARAMS.keys())))
            self.host_params = HOST_PARAMS["lxplus"]
        else:
            self.host_params = HOST_PARAMS[self.host] 

        self.job_type = CondorJob
        self.output_dir = os.path.abspath(output_dir)
        self.submit_all_file = self.output_dir + "/submit_all.txt"
        self.job_map = {} # run condor_q only once and store results in this map, rather than running condor_q for each job individually. self.job_map then gets passed to the Tasks, where they do the rest of the bookkeeping.

        if self.output_dir.startswith("/eos") and self.host == "lxplus":
            self.log_output_dir = self.higgs_dna_path + "/eos_logs/" + os.path.split(self.output_dir)[-1]
            os.system("mkdir -p %s" % self.log_output_dir) 
            logger.warning("[CondorManager : __init__] It appears you are running on lxplus and have specified an /eos area as your output dir: '%s'. Since log files for condor jobs are not allowed to be placed in the /eos area, we will write your log files to a separate area : '%s'. See https://batchdocs.web.cern.ch/troubleshooting/eos.html for more details." % (self.output_dir, self.log_output_dir))

        # On some host sites, like T2s, condor jobs cannot directly access the local cluster and jobs need to be output to a special area, e.g. /hadoop 
        if "condor_base_path" in self.host_params.keys():
            self.batch_output_dir = os.path.abspath( # always do abspaths to be very safe
                    self.host_params["condor_base_path"].replace("USERNAME_INITIAL", self.user[0]).replace("USERNAME", self.user).replace("DIRNAME", os.path.basename(self.output_dir))
            )
            if os.path.exists(self.batch_output_dir):
                logger.warning("[CondorManager : __init__] Output files from jobs will be written to directory '%s' which appears to already exist. If there are output files from an old run of HiggsDNA here, this may result in unintended behavior! We will assume this is the same analysis and that any output files there should be counted as 'done', rather than overwriting them. If this is not what you want, please delete the old files or specify a new --output_dir." % (self.batch_output_dir))
        else: 
            self.batch_output_dir = self.output_dir

        for redirector in ["xrd_redirector", "gfal_redirector"]:
            if redirector in self.host_params.keys():
                to_replace, replace_with = self.host_params[redirector]
                setattr(
                        self,
                        redirector.split("_")[0] + "_batch_output_dir",
                        self.batch_output_dir.replace(to_replace, replace_with)
                )

        for dir in [self.output_dir, self.batch_output_dir]:
            os.system("mkdir -p %s" % dir)


    def customize_jobs(self):
        """

        """
        with tqdm(total = sum([len(task.jobs) for task in self.tasks])) as prog:
            for task in self.tasks:
                for job in task.jobs:
                    job.__class__ = self.job_type
                    job.conda_tarfile = self.conda_tarfile
                    job.analysis_tarfile = self.analysis_tarfile
                    job.proxy = self.proxy    
                    job.host = self.host
                    job.host_params = self.host_params
                    job.set_output_dir(base = self.batch_output_dir)

                    if hasattr(self, "log_output_dir"):
                        job.set_log_output_dir(base = self.log_output_dir)

                    # If we are copying tarfiles via xrd, pass these paths to the jobs
                    if self.host_params["needs_tar"] and "xrd_redirector" in self.host_params.keys():
                        job.xrd_conda_tarfile = self.xrd_conda_tarfile
                        job.xrd_analysis_tarfile = self.xrd_analysis_tarfile
                    job.set_files() # update the locations for the various files associated with this job
                    prog.update(1)
        self.customized_jobs = True


    def submit_to_batch(self, dry_run = False):
        """

        """
        self.monitor_jobs()

        jobs_to_submit = []
        for task in self.tasks:
            task.process(job_map = self.job_map)
            jobs_to_submit += [job for job in task.jobs if ((job.status == "waiting" or job.status == "failed") and job.submit(dry_run = True))]

        # merge all individual submit files into a single one
        if dry_run:
            return 

        if jobs_to_submit: 
            submit_files = [job.condor_submit_file for job in jobs_to_submit]
            submit_chunks = create_chunks(submit_files, 100) # submit in batches of 100

            idx = 0
            for i, chunk in enumerate(submit_chunks):
                submit_file = self.submit_all_file.replace(".txt", "_%d.txt" % i)
                command = "cat %s > %s" % (" ".join(chunk), submit_file) 
                os.system(command)
                
                # submit all jobs
                results = do_cmd("condor_submit %s" % submit_file).split("\n")
                n_submitted = 0
                for line in results:
                    if "1 job(s) submitted to cluster" in line:
                        n_submitted += 1

                # Check if they were submitted successfully
                if n_submitted == len(chunk):
                    logger.info("[CondorManager : submit_to_batch] Submitted %d jobs." % n_submitted)

                else:
                    logger.exception("[CondorManager : submit_to_batch] We found %d jobs to submit, but only %d were successfully submitted." % (len(jobs_to_submit), n_submitted))
                    raise RuntimeError()

                # Assign cluster id's to jobs
                for line in results:
                    if "1 job(s) submitted to cluster" in line:
                        jobs_to_submit[idx].cluster_id = line.split("cluster")[-1][:-1].strip()
                        jobs_to_submit[idx].n_attempts += 1
                        self.job_map[jobs_to_submit[idx].cluster_id] = None
                        idx += 1


        self.monitor_jobs()

        for task in self.tasks:
            task.process(job_map = self.job_map)
            self.summary[task.name] = task.summary 


    def monitor_jobs(self):
        """

        """
        self.job_map = {}
        monitor_results = do_cmd("condor_q --json -attributes ClusterId,JobStatus")

        if not monitor_results: # no condor jobs were found
            return

        success = False
        while not success:
            try:
                monitor_results = json.loads(monitor_results)
                success = True
            except: # sometimes condor_q just intermittently fails, just try again later if this happens
                monitor_results = do_cmd("condor_q --json -attributes ClusterId,JobStatus")

        for x in monitor_results:
            self.job_map[str(x["ClusterId"])] = CONDOR_STATUS_FLAGS[x["JobStatus"]]


    def prepare_inputs(self):
        self.conda_tarfile = self.output_dir + "/" + "higgs-dna.tar.gz"
        self.analysis_tarfile = self.output_dir + "/" + "higgs_dna.tar.gz"

        # Make tar files only if we need them and they don't already exist
        if self.host_params["needs_tar"]:
            self.needs_conda_tarfile = not os.path.exists(self.conda_tarfile)
            self.needs_analysis_tarfile = not os.path.exists(self.analysis_tarfile)

            if not self.needs_conda_tarfile:
                tar_size = os.path.getsize(self.conda_tarfile) * (1. / 1024.)**3
                logger.warning("[CondorManager : prepare_inputs] conda pack '%s' of size %.3f GB already exists, not remaking." % (self.conda_tarfile, tar_size))
            elif os.path.exists((CONDA_TARFILE)):
                self.conda_tarfile = (CONDA_TARFILE)
                tar_size = os.path.getsize(self.conda_tarfile) * (1. / 1024.)**3
                logger.warning("[CondorManager : prepare_inputs] Found old conda pack '%s' of size %.3f GB, will use this one. If you have installed new packages since this was made, you should delete it and rerun (and a new one will be made automatically). If you have only modified code under HiggsDNA/ it is probably not necessary to remake." % (self.conda_tarfile, tar_size))
            else:
                conda_pack_command = "conda pack -n higgs-dna --ignore-editable-packages --ignore-missing-files -o %s --compress-level 5 --n-threads 12" % self.conda_tarfile # compression level of 5 seems like a good balance of speed and size reduction: compression of 1 gives a tarfile of size 463MB, while 5 gives 413MB in 30s, while 9 gives 409MB in 91s 
                logger.info("[CondorManager : prepare_inputs] Making conda pack '%s' with command '%s'" % (self.conda_tarfile, conda_pack_command))
                logger.info("Note: the ``conda pack`` command can sometimes take a long time. If you are getting annoyed with how long it takes, you can always manually copy old tar files to your ``--output_dir`` area and HiggsDNA will not remake them. This will only work if you have not installed new python packages to your conda env.")
                os.system(conda_pack_command)
                logger.info(self.conda_tarfile)
                logger.info((CONDA_TARFILE))
                logger.info("Copying conda tarfile to area '%s' to be reused for future runs. If you install new packages and need to remake the tarfile, simply delete this file and a new conda tarfile will be made for the next run." % ((CONDA_TARFILE)))
                os.system("cp %s %s" % (self.conda_tarfile, (CONDA_TARFILE)))


            if not self.needs_analysis_tarfile:
                tar_size = os.path.getsize(self.analysis_tarfile) * (1. / 1024.)**3
                logger.warning("[CondorManager : prepare_inputs] Analysis tarfile '%s' of size %.4f GB already exists, not remaking." % (self.analysis_tarfile, tar_size))

            else:
                tar_options = ""
                tar_command = "XZ_OPT='-1e -T12' tar -Jc %s -f %s -C %s higgs_dna jsonpog-integration data metadata" % (tar_options, self.analysis_tarfile, self.higgs_dna_path)
                logger.info("[CondorManager : prepare_inputs] Making analysis tarfile '%s' with command '%s'." % (self.analysis_tarfile, tar_command))

                t_start_tar = time.time()
                os.system(tar_command)
                t_elapsed_tar = time.time() - t_start_tar
                
                tarfile_size = os.path.getsize(self.analysis_tarfile) * (1. / 1024.)**3
                logger.debug("[CondorManager : prepare_inputs] Made analysis tarfile '%s' of size %.3f GB in %.1f s" % (self.analysis_tarfile, tarfile_size, t_elapsed_tar))

            # If an xrootd redirector is specified for this host, we can copy the tarfiles to /hadoop and xrdcp them into jobs
            if "xrd_redirector" in self.host_params.keys():
                self.batch_conda_tarfile = self.batch_output_dir + "/" + "higgs-dna.tar.gz"
                self.batch_analysis_tarfile = self.batch_output_dir + "/" + "higgs_dna.tar.gz"

                if not self.copied_tars: # the output dir is already the same as the batch output dir
                    logger.debug("[CondorManager : prepare_inputs] Transferring tarfiles to directory '%s'." % self.batch_output_dir)
                    for x in [self.batch_conda_tarfile, self.batch_analysis_tarfile]:
                        if os.path.exists(x): # delete any old versions before copying, giving priority to the freshly made ones
                            logger.warning("[CondorManager : prepare_inputs] Found existing tarfile '%s' in batch output directory, will not overwrite this. Manually delete this file if you want it to be overwritten." % (x))

                            #os.system("rm %s" % x)

                    if self.host_params["copy_tar"] == "cp":
                        if not os.path.exists(self.batch_conda_tarfile):
                            os.system("cp %s %s" % (self.conda_tarfile, self.batch_conda_tarfile))
                        if not os.path.exists(self.batch_analysis_tarfile):
                            os.system("cp %s %s" % (self.analysis_tarfile, self.batch_analysis_tarfile))
                        os.system("chmod 755 %s" % (self.batch_output_dir + "/*.tar.gz"))

                    elif self.host_params["copy_tar"] == "xrd":
                        x,y = self.host_params["xrd_redirector"]
                        if not os.path.exists(self.batch_conda_tarfile):
                            os.system("xrdcp %s %s" % (self.conda_tarfile, self.batch_conda_tarfile.replace(x,y)))
                        if not os.path.exists(self.batch_analysis_tarfile):
                            os.system("xrdcp %s %s" % (self.analysis_tarfile, self.batch_analysis_tarfile.replace(x,y)))

                    self.copied_tars = True

                if "xrd_redirector" in self.host_params.keys():
                    to_replace, replace_with = self.host_params["xrd_redirector"]
                    self.xrd_conda_tarfile = self.batch_conda_tarfile.replace(to_replace, replace_with)
                    self.xrd_analysis_tarfile = self.batch_analysis_tarfile.replace(to_replace, replace_with)
                    logger.debug("[CondorManager : prepare_inputs] Will access tarfiles via xrootd with mapping:")
                    logger.debug("\t '%s' -> '%s'" % (self.batch_conda_tarfile, self.xrd_conda_tarfile))
                    logger.debug("\t '%s' -> '%s'" % (self.batch_analysis_tarfile, self.xrd_analysis_tarfile))

        # Check grid proxy
        self.proxy = check_proxy()
        if self.proxy is None:
            logger.exception("[CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found to be expired.")
            raise RuntimeError()

        if self.host_params["needs_copy_proxy"]:
            logger.debug("[CondorManager : prepare_inputs] Copying proxy %s to home area %s in order to copy to jobs." % (self.proxy, self.higgs_dna_path))
            os.system("cp %s %s" % (self.proxy, self.higgs_dna_path))

        self.prepared_inputs = True


    def convert_to_local(self):
        """
        Convert itself into a LocalManager
        """

        raise NotImplementedError()
        """
        local_manager = LocalManager(output_dir = self.output_dir)
        for k,v in self.__dict__.items():
            if k in ["job_type"]:
                continue
            setattr(local_manager, k, v)

        local_manager.customize_jobs()
        local_manager.first_submit = True

        return local_manager
        """
