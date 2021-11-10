import os
import importlib
import time
import awkward
import numpy
import json
import glob

import logging
logger = logging.getLogger(__name__)

from higgs_dna.job_management.jobs import Job
from higgs_dna.utils.metis_utils import num_to_ordinal_string
from higgs_dna.utils import awkward_utils

class JobsManager():
    """

    """
    def __init__(self): 
        self.tasks = []

        self.prepared_inputs = False
        self.complete = False
        self.summary = {}

        self.job_type = Job
        self.jobs = []
        
        self.remerge = False

    def add_task(self, task):
        self.tasks.append(task)


    def submit_jobs(self, summarize):
        """

        """

        if not self.prepared_inputs:
            self.prepare_inputs()

        self.customize_jobs()

        self.submit_to_batch()

        """
        for task in self.tasks:
            task.process()

            jobs_to_submit = [job for job in task.jobs if (job.status == "waiting" or job.status == "failed")]
            for job in jobs_to_submit:
                if self.can_submit(job):
                    submitted = job.submit()
                    if submitted:
                        logger.info("[JobsManager : submit_jobs] : %s submitted job '%s_%d' (for the %s time)" % (str(type(self)), job.name, job.idx, num_to_ordinal_string(job.n_attempts)))
                        os.system("sleep 0.1s")

            task.process()

            self.summary[task.name] = task.summary 
        """

        self.complete = all([task.complete for task in self.tasks])
        if summarize:
            self.summarize()

        return self.summary


    def submit_to_batch(self):
        """

        """
        pass


    def summarize(self):
        """

        """
        logger.info("\n\n[JobsManager : summarize] Summarizing task progress.\n")
        for task, info in self.summary.items():
            task_info = info["jobs"]
            message = "[JobsManager : summarize] Task '%s' : \n" % task
            for status_set in [["completed"], ["waiting", "idle"], ["running"], ["failed"], ["retired"]]:
                n_jobs = 0
                for status in status_set:
                    n_jobs += task_info[status]
                if n_jobs == 0: # don't show e.g. "0/X jobs completed"
                    continue
                frac = float(n_jobs) / float(task_info["all"])
                status_str = status_set[0] if len(status_set) == 1 else "/".join(status_set)
                message += "\t%d/%d (%.2f percent) jobs %s\n" % (n_jobs, task_info["all"], frac * 100., status_str)

            logger.info(message)


    def merge_outputs(self, dir):
        """

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
                logger.debug("[JobsManager : merge_outputs] For syst_tag '%s' merged output file '%s' already exists, but some jobs were unretired, so we will remerge inputs and overwrite this file." % (syst_tag, merged_file))

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
            for f in outputs:
                logger.debug("\t %s" % f)
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

        """
        return True


from higgs_dna.job_management.jobs import LocalJob

class LocalManager(JobsManager):
    """

    """
    def __init__(self, **kwargs):
        super(LocalManager, self).__init__()

        self.n_cores = kwargs.get("n_cores", 4)
        self.n_running_jobs = 0
        self.job_type = LocalJob


    def customize_jobs(self):
        """

        """
        for task in self.tasks:
            for job in task.jobs:
                if isinstance(job, self.job_type):
                    continue
                job.__class__ = self.job_type


    def submit_to_batch(self):
        """

        """
        for task in self.tasks:
            task.process()

            jobs_to_submit = [job for job in task.jobs if (job.status == "waiting" or job.status == "failed")]
            for job in jobs_to_submit:
                if self.can_submit(job):
                    submitted = job.submit()
                    if submitted:
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


from higgs_dna.job_management.jobs import CondorJob
from higgs_dna.utils.metis_utils import do_cmd
from higgs_dna.utils.misc_utils import check_proxy
from higgs_dna.constants import CONDOR_STATUS_FLAGS

class CondorManager(JobsManager):
    """

    """
    def __init__(self, **kwargs):
        super(CondorManager, self).__init__()

        self.job_type = CondorJob
        self.output_dir = kwargs.get("output_dir")
        self.batch_output_dir = kwargs.get("batch_output_dir")
        self.submit_all_file = self.output_dir + "/submit_all.txt"

        host = do_cmd("hostname")
        if "t2.ucsd" in host:
            self.host = "UCSD"

        elif "lxplus" in host:
            self.host = "lxplus"

        self.job_map = {} # run condor_q only once and store results in this map, rather than running condor_q for each job individually. self.job_map then gets passed to the Tasks, where they do the rest of the bookkeeping.


    def customize_jobs(self):
        """

        """
        for task in self.tasks:
            for job in task.jobs:
                if isinstance(job, self.job_type):
                    continue
                job.__class__ = self.job_type
                job.conda_tarfile = self.conda_tarfile
                job.analysis_tarfile = self.analysis_tarfile
                job.proxy = self.proxy    
                job.host = self.host
                if self.host in ["UCSD"]:
                    job.config["remote_job"] = True 


    @staticmethod
    def create_chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        chunks = []
        for i in range(0, len(lst), n):
            chunks.append(lst[i:i+n])
        return chunks 


    def submit_to_batch(self):
        """

        """
        self.monitor_jobs()

        jobs_to_submit = []
        for task in self.tasks:
            task.process(job_map = self.job_map)
            jobs_to_submit += [job for job in task.jobs if ((job.status == "waiting" or job.status == "failed") and job.submit(dry_run = True))]

        # merge all individual submit files into a single one
        if jobs_to_submit: 
            submit_files = [job.batch_submit_file for job in jobs_to_submit]
            submit_chunks = self.create_chunks(submit_files, 100) # submit in batches of 100

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
        monitor_results = do_cmd("condor_q --json")

        if not monitor_results: # no condor jobs were found
            return

        else:
            monitor_results = json.loads(monitor_results)

        for x in monitor_results:
            self.job_map[str(x["ClusterId"])] = CONDOR_STATUS_FLAGS[x["JobStatus"]]


    def prepare_inputs(self):
        self.conda_path = do_cmd("echo $CONDA_PREFIX").replace("/higgs-dna", "")
        self.higgs_dna_path = do_cmd("pwd").split("HiggsDNA")[0] + "HiggsDNA/"

        self.conda_tarfile = self.output_dir + "/" + "higgs-dna.tar.gz"
        self.analysis_tarfile = self.output_dir + "/" + "higgs_dna.tar.gz"
        self.tar_options = ""

        if os.path.exists(self.conda_tarfile) and os.path.exists(self.analysis_tarfile):
            tarfile_size = os.path.getsize(self.conda_tarfile) * (1. / 1024.)**3
            analysis_tarfile_size = os.path.getsize(self.analysis_tarfile) * (1. / 1024.)**3
            logger.warning("[CondorManager : prepare_inputs] conda pack '%s' of size %.3f GB and analysis tarfile '%s' of size %.3f GB already exists, not remaking." % (self.conda_tarfile, tarfile_size, self.analysis_tarfile, analysis_tarfile_size))

        else:
            if not os.path.exists(self.conda_tarfile):
                conda_pack_command = "conda pack -n higgs-dna --ignore-editable-packages -o %s --compress-level 5 --n-threads 12" % self.conda_tarfile # compression level of 5 seems like a good balance of speed and size reduction: compression of 1 gives a tarfile of size 463MB, while 5 gives 413MB in 30s, while 9 gives 409MB in 91s 
                logger.info("[CondorManager : prepare_inputs] Making conda pack '%s' with command '%s'" % (self.conda_tarfile, conda_pack_command))
                os.system(conda_pack_command)

            if not os.path.exists(self.analysis_tarfile):
                tar_command = "XZ_OPT='-1e -T12' tar -Jc %s -f %s -C %s higgs_dna jsonpog-integration" % (self.tar_options, self.analysis_tarfile, self.higgs_dna_path)
                logger.info("[CondorManager : prepare_inputs] Making analysis tarfile '%s' with command '%s'." % (self.conda_tarfile, tar_command))

                t_start_tar = time.time()
                os.system(tar_command)
                t_elapsed_tar = time.time() - t_start_tar
                
                tarfile_size = os.path.getsize(self.analysis_tarfile) * (1. / 1024.)**3
                logger.debug("[CondorManager : prepare_inputs] Made tarfile '%s' of size %.3f GB in %.1f s" % (self.analysis_tarfile, tarfile_size, t_elapsed_tar))


        if self.host in ["UCSD"]:
            self.conda_tarfile_batch = self.batch_output_dir + "/" + "higgs-dna.tar.gz"
            self.analysis_tarfile_batch = self.batch_output_dir + "/" + "higgs_dna.tar.gz"

            #if not (os.path.exists(self.conda_tarfile_batch) or os.path.exists(self.analysis_tarfile_batch)):
            logger.debug("[CondorManager : prepare_inputs] Transferring tarfiles to hadoop directory '%s' so they may be copied with xrd to reduce I/O load on local cluster." % self.batch_output_dir)
            os.system("cp %s %s" % (self.conda_tarfile, self.conda_tarfile_batch))
            os.system("cp %s %s" % (self.analysis_tarfile, self.analysis_tarfile_batch))

            logger.debug("[CondorManager : prepare_inputs] Setting replication factor to 30 to increase transfer speed.")
            os.system("hadoop fs -setrep -R 30 %s" % (self.conda_tarfile_batch.replace("/hadoop","")))
            os.system("hadoop fs -setrep -R 30 %s" % (self.analysis_tarfile_batch.replace("/hadoop","")))
 

        # Check grid proxy
        self.proxy = check_proxy()
        if self.proxy is None:
            logger.exception("[CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found to be expired.")
            raise RuntimeError()

        if self.host == "lxplus":
            logger.debug("[CondorManager : prepare_inputs] Copying proxy %s to home area %s in order to copy to jobs." % (self.proxy, self.higgs_dna_path))
            os.system("cp %s %s" % (self.proxy, self.higgs_dna_path))

        self.prepared_inputs = True
