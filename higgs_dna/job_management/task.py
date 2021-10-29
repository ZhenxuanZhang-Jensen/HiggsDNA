"""
Much of the code is taken from:
    ProjectMetis metis/Task.py: https://github.com/aminnj/ProjectMetis/blob/f9e71556cb84496731fa71dcab2dfc82b6e3022f/metis/Task.py
    Author: Nick Amin
"""
import copy
import os
import math
import json
import awkward
import numpy
from tqdm import tqdm

import logging
logger = logging.getLogger(__name__)

from higgs_dna.job_management.jobs import Job 
from higgs_dna.utils import awkward_utils
from higgs_dna.constants import CENTRAL_WEIGHT

class Task():
    """

    """
    def __init__(self, name, output_dir, batch_output_dir, files, config, **kwargs):
        self.name = name
        self.output_dir = output_dir
        self.batch_output_dir = batch_output_dir
        self.files = files
        self.config = config
        self.job_config = copy.deepcopy(config)
        self.job_config["sample"]["files"] = None # lightweight config to send to jobs

        self.min_completion_frac = kwargs.get("min_completion_frac", 1.0)
        self.n_files_per_job = kwargs.get("n_files_per_job", 3)
        self.max_jobs = kwargs.get("max_jobs", None) 
        if self.max_jobs is None:
            self.max_jobs = -1

        for key, val in kwargs.items():
            setattr(self, key, val)

        for dir in [self.output_dir, self.batch_output_dir]:
            os.system("mkdir -p %s" % dir)

        self.complete = False
        self.wrote_process_ids = False
        self.wrote_years = False

        self.jobs = []
        self.create_jobs()

        self.summary = {}
        self.outputs = {}
        self.merged_outputs = {}
        self.summarize()


    def create_jobs(self):
        """

        """
        logger.info("[Task: create_jobs] Task '%s' : splitting %d input files into %d jobs" % (self.name, len(self.files), math.ceil(float(len(self.files)) / float(self.n_files_per_job)))) 
        file_splits = create_chunks(self.files, self.n_files_per_job)
        
        for idx, file_split in enumerate(tqdm(file_splits)):
            if self.max_jobs >= 0:
                if idx >= self.max_jobs:
                    continue
            job_config = copy.deepcopy(self.job_config)

            job = Job(
                name = self.name,
                function = copy.deepcopy(self.config["function"]),
                inputs = file_split,
                output_dir = self.output_dir,
                batch_output_dir = self.batch_output_dir,
                idx = idx + 1,
                config = job_config 
            )
            self.jobs.append(job)


    def process(self, job_map = None):
        """

        """
        if self.complete:
            return
        
        # Check status of all jobs
        for job in self.jobs:
            if job_map is not None:
                if os.path.exists(job.summary_file):
                    job.status = "completed"
                elif job.cluster_id in job_map.keys():
                    if job_map[job.cluster_id] is not None:
                        job.status = job_map[job.cluster_id]
                        continue
                elif job.status == "retired":
                    continue
                else:
                    job.status = "failed"
            else:
                job.query_status()

        # Are we done?
        self.n_completed_jobs = len([job for job in self.jobs if job.status == "completed"]) 
        self.n_retired_jobs = len([job for job in self.jobs if job.status == "retired"]) 
        self.completion_frac = float(self.n_completed_jobs) / float(len(self.jobs))
        self.completion_or_retired_frac = float(self.n_completed_jobs + self.n_retired_jobs) / float(len(self.jobs))

        # Did we successfully process minimum fraction of completed jobs?
        if self.completion_frac >= self.min_completion_frac:
            logger.info("[Task : process] Task '%s' COMPLETED : %d/%d (%.2f percent) of jobs completed, which is >= the minimum job completion fraction for this task (%.2f percent)." % (self.name, self.n_completed_jobs, len(self.jobs), 100. * self.completion_frac, 100. * self.min_completion_frac))
            self.complete = True # Done

        # Otherwise, are the uncompleted jobs "retired" (meaning they failed up to the maximum number of retries)?
        # If so, we will mark this as done but give you a warning about the retired jobs.
        # You can resubmit these jobs with the option --resubmit_retired in run_analysis.py
        elif self.completion_or_retired_frac >= self.min_completion_frac:
            logger.info("[Task : process] Task '%s' COMPLETED : %d/%d (%.2f percent) of jobs completed/retired which is >= the minimum job completion fraction for this task (%.2f percent)." % (self.name, self.n_completed_jobs + self.n_retired_jobs, len(self.jobs), 100. * self.completion_or_retired_frac, 100. * self.min_completion_frac))
            retired_jobs = [job for job in self.jobs if job.status == "retired"]
            for job in retired_jobs:
                logger.warning("[Task : process] WARNING: Task '%s' had to retire job '%s' since it ran unsuccessfully for %d straight times. If this is an MC sample, this will just reduce your statistics. If this is a data job, you have processed less events than you intended!" % (self.name, job.name_full, job.n_attempts))
            self.complete = True

        # If neither of the first two, we are not done yet
        else:
            self.summarize()
            return

        # Clean up: kill any idle or running jobs
        n_killed = 0
        jobs_to_kill = [job for job in self.jobs if (job.status == "idle" or job.status == "running")]
        for job in jobs_to_kill:
            logger.debug("[Task : process] Task '%s' : since Task is COMPLETED, we are killing job '%s' with status '%s'" % (self.name, job.name_full, job.status))
            job.kill()
        if n_killed > 0:
            logger.info("[Task : process] Task '%s' : since Task is COMPLETED, we killed all idle and running jobs (%d jobs killed)" % (self.name, n_killed))
        self.complete = True

        self.summarize()
            
    
    def unretire_jobs(self):
        """

        """
        retired_jobs = [job for job in self.jobs if job.status == "retired"]
        if len(retired_jobs) == 0:
            return False

        logger.info("[Task : unretire_jobs] Task '%s' : since you are re-running with option --resubmit_retired, we are resubmitting %d jobs which were previously retired (meaning they repeatedly failed up to the number of max attempts)." % (self.name, len(retired_jobs)))
        for job in retired_jobs:
            logger.debug("[Task : unretire_jobs] Task '%s' : resubmitting job '%s' with status '%s' and %d previously failed attempts." % (self.name, job.name_full, job.status, job.n_attempts))
            job.n_attempts = 0
            job.status = "waiting"
    
        self.complete = False
        self.wrote_years = False
        self.wrote_process_ids = False
        return True


    def summarize(self):
        """

        """
        # Summarize jobs
        job_summary = {}
        job_summary["all"] = len(self.jobs)
        for status in ["waiting", "idle", "running", "failed", "completed", "retired"]:
            job_summary[status] = len([job for job in self.jobs if job.status == status])
        self.summary["jobs"] = job_summary

        if not self.complete:
            return

        # Summarize physics and performance
        phys_summary = {
                "n_events_initial" : 0,
                "n_events_selected" : {},
                "sum_weights" : 0.0,
        }

        performance = {
                "time" : 0.0,
                "time_load" : 0.0,
                "time_syst" : 0.0,
                "time_taggers" : 0.0
        }

        self.outputs = {}

        completed_jobs = [job for job in self.jobs if job.status == "completed"]
        for job in completed_jobs:
            with open(job.summary_file, "r") as f_in:
                job_info = json.load(f_in)
         
            phys_summary["n_events_initial"] += job_info["n_events"]
            if not self.config["sample"]["is_data"]:
                phys_summary["sum_weights"] += job_info["sum_weights"]

            for syst_tag, n_events in job_info["n_events_selected"].items():
                if syst_tag not in phys_summary["n_events_selected"].keys():
                    phys_summary["n_events_selected"][syst_tag] = n_events
                else:
                    phys_summary["n_events_selected"][syst_tag] += n_events

            for syst_tag, output in job_info["outputs"].items():
                if syst_tag not in self.outputs.keys():
                    self.outputs[syst_tag] = []
                if job_info["config"]["remote_job"]:
                    output = job_info["config"]["output_dir"] + output
                if not job_info["n_events_selected"][syst_tag] > 0: # skip empty parquet files to avoid errors
                    continue
                else:
                    self.outputs[syst_tag].append(output)

            performance["time"] += job_info["time"]
            for portion in ["load", "syst", "taggers"]:
                performance["time_%s" % portion] += float(job_info["time"]) * job_info["time_frac_%s" % portion]

        if not self.config["sample"]["is_data"]:
            phys_summary["norm_factor"] = self.config["sample"]["norm_factor"]
            if phys_summary["sum_weights"] > 0.:
                phys_summary["scale1fb"] = (self.config["sample"]["norm_factor"] * 1000.) / phys_summary["sum_weights"]
            else:
                phys_summary["scale1fb"] = 0.
            self.lumi = self.config["sample"]["lumi"]
            self.scale1fb = phys_summary["scale1fb"]
        self.summary["physics"] = phys_summary
        self.summary["performance"] = performance

        self.merge_outputs()
        if not self.config["sample"]["is_data"]:
            self.apply_scale1fb_and_lumi()


    def merge_outputs(self):
        """

        """
        self.merged_outputs = {}
        for syst_tag, outputs in self.outputs.items():
            if not outputs:
                continue
            
            merged_output = self.output_dir + "/merged_%s.parquet" % (syst_tag) 
            self.merged_outputs[syst_tag] = merged_output
            merged_events = []

            for output in outputs:
                merged_events.append(awkward.from_parquet(output))

            logger.debug("[Task : merge_outputs] Task '%s' : merging %d outputs into file '%s'." % (self.name, len(outputs), merged_output))

            merged_events = awkward.concatenate(merged_events)
            awkward.to_parquet(merged_events, merged_output)


    def apply_scale1fb_and_lumi(self):
        """

        """
        for syst_tag, merged_output in self.merged_outputs.items():
            events = awkward.from_parquet(merged_output)

            for field in events.fields:
                # only apply scale1fb and lumi to central weight so variations are standalone
                if field == CENTRAL_WEIGHT:
                    logger.debug("[Task : apply_scale1fb_and_lumi] Task '%s' : scaling central weight branch '%s' in output file '%s' by scale1fb (%.9f) times lumi (%.2f)" % (self.name, field, merged_output, self.scale1fb, self.lumi))
                    awkward_utils.add_field(
                            events = events,
                            name = field,
                            data = events[field] * self.scale1fb * self.lumi,
                            overwrite = True
                    )

            awkward.to_parquet(events, merged_output)


    def add_process_id(self):
        """

        """
        if self.wrote_process_ids:
            return

        self.process_id = self.config["sample"]["process_id"]
        if self.process_id is None:
            return

        for syst_tag, merged_output in self.merged_outputs.items():
            events = awkward.from_parquet(merged_output)
            logger.debug("[Task : add_process_id] Task '%s' : adding field 'process_id' with value %d in output file '%s'" % (self.name, self.process_id, merged_output))
            awkward_utils.add_field(
                    events = events,
                    name = "process_id",
                    data = numpy.ones(len(events)) * self.process_id
            )
            awkward.to_parquet(events, merged_output)
        

        self.wrote_process_ids = True

    
    def add_year(self):
        """

        """
        if self.wrote_years:
            return

        self.year = int(self.config["sample"]["year"])

        for syst_tag, merged_output in self.merged_outputs.items():
            events = awkward.from_parquet(merged_output)
            logger.debug("[Task : add_process_id] Task '%s' : adding field 'year' with value %d in output file '%s'." % (self.name, self.year, merged_output))
            
            awkward_utils.add_field(
                    events = events,
                    name = "year",
                    data = numpy.ones(len(events)) * self.year
            )

            awkward.to_parquet(events, merged_output)


        self.wrote_years = True


def create_chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    chunks = []
    for i in range(0, len(lst), n):
        chunks.append(lst[i:i+n])
    return chunks

