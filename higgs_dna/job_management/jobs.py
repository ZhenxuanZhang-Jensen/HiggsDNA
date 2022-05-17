import os
import json
import subprocess

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils.metis_utils import do_cmd
from higgs_dna.utils.misc_utils import get_HiggsDNA_base, get_HiggsDNA_conda
from higgs_dna.job_management.constants import CONDOR_STATUS_FLAGS, HOST_PARAMS

class Job():
    """
    A ``Job`` instance will be created by its parent ``Task``.
    The base ``Job`` class is purely virtual, and it will be converted into either a ``LocalJob`` or a ``CondorJob`` by its parent ``JobsManager`` (which is also the parent of the ``Task``).

    A ``LocalJob`` has the following associated files under its ``dir``: 
        - config [<samplename>_<year>_config_job<n>.json] : same as ``config`` arg
        - python executable [<samplename>_<year>_executable_job<n>.py] : python script to run this job with ``python <executable_file>.py``
        - parquet file(s) [<output_job_<n>_<syst_var>.parquet] : actual physics content
        - summary [<samplename>_<year>_summary_job<n>.json] : FIXME what do we actually need here? Look in task.py

    A ``CondorJob`` has the following additional associated files:
        - condor executable [<samplename>_<year>_executable_job<n>.sh] : shell script to run the python executable on a condor node
        - condor submit [<samplename>_<year>_batch_submit_job<n>.txt] : file to submit the job to HTCondor with ``condor_submit <submit_file>.txt``
        - condor logs [<cluster_id>.<process_id>.{out,err,log}] : the ``.out`` and ``.err`` contain the ``stdout`` and ``stderr`` and the ``.log`` contains information printed by HTCondor about the job

    For a ``CondorJob``, all of these files will not always be in the same directory, depending on the host site and actual path of the ``dir``:
        - config, executables (.py and .sh), condor submit file, and logs will always go in ``dir``
        - summary and parquet files may be in a different directory, called ``output_dir``

    Due to the complication of a ``CondorJob``, every job has two directory attributes: ``dir`` and ``output_dir``:
        - ``dir == output_dir`` for every ``LocalJob``
        - may not be the case for a ``CondorJob``

    Making this distinction, which is not necessary for a ``LocalJob`` allows the logic for any ``Job`` in other parts of ``higgs_dna`` to remain the same.

    :param name: the name of the parent ``Task``, typically in the format <sample_name>_<year>
    :type name: str
    :param inputs: list of input files for this ``Job`` to run over
    :type inputs: list of str
    :param function: the ``higgs_dna`` function to be run by this job
    :type function: function
    :param config: dictionary that specifies the physics content (samples, tag sequence, systematics) of a given job
    :type config: dict
    :param dir: output directory of the parent ``Task``
    :type dir: str
    :param idx: integer to identify this job among the other jobs of the parent ``Task``
    :type idx: int
    """
    def __init__(self, name, inputs, function, config, dir, idx):
        self.name = name
        self.inputs = inputs
        self.function = function
        self.config = config
        self.idx = idx
        self.name_full = self.name + "_" + str(self.idx)
        self.n_attempts = 0
        self.force_retirement = False

        # If the ``dir`` for an analysis is /some/path/dir, the job will have an output directory of the format /some/path/dir/samplename_year/job_n
        self.task_dir = os.path.abspath(dir)
        self.dir = os.path.abspath(dir + "/job_%s/" % (str(self.idx)))
        os.system("mkdir -p %s" % self.dir)
        self.output_name = "output_job_%d" % self.idx
        
        self.output_dir = None # output_dir should be set in a derived LocalJob or CondorJob by a LocalManager or CondorManager, respectively
        self.log_output_dir = self.dir # this may be separate for lxplus jobs writing outputs to /eos, see https://batchdocs.web.cern.ch/troubleshooting/eos.html

        self.wrote_config = False
        self.wrote_python_executable = False
        self.wrote_condor_files = False
        self.cluster_id = None
        self.can_submit = True
        self.status = "waiting"
        self.processed = False # set to True by Task after recording its metadata (to avoid double counting jobs)


    def set_output_dir(self, base):
        """
        Set the output directory for this job from the base directory of the JobsManager.

        :param base: output directory path up to the "/sample_year/job_n" part (this will be added manually by Job)
        :type base: str
        """
        self.output_dir = os.path.abspath(base) + "/" + "/".join(os.path.abspath(self.dir).split("/")[-2:]) # tack on last 2 subdir
        os.system("mkdir -p %s" % self.output_dir)
        self.output_dir = os.path.abspath(self.output_dir)

    def set_log_output_dir(self, base):
        """
        Set the output directory for logs of this job from the base directory of the JobsManager.
        For why this might be necessary, see:
            - https://batchdocs.web.cern.ch/troubleshooting/eos.html

        :param base: output directory path up to the "/sample_year/job_n" part (this will be added manually by Job)
        :type base: str
        """ 
        self.log_output_dir = os.path.abspath(base) + "/" + "/".join(os.path.abspath(self.dir).split("/")[-2:]) # tack on last 2 subdir
        os.system("mkdir -p %s" % self.log_output_dir)
        self.log_output_dir = os.path.abspath(self.log_output_dir)

    def set_files(self):
        """
        Function to update all file paths associated with this job.
        In general, we can rely on config_file to also have executables, condor_submit, and logs with it
        and rely on summary_file to also have parquets with it.        
        """
        self.dir = os.path.abspath(self.dir)
        self.output_dir = os.path.abspath(self.output_dir)

        self.config_file = self.dir + "/%s_config_job%d.json" % (self.name, self.idx)
        self.python_executable_file = self.dir + "/%s_executable_job%d.py" % (self.name, self.idx)
        self.log_file = self.log_output_dir + "/%s_log_job%d.log" % (self.name, self.idx) # only actually logs here for LocalJob
        self.summary_file = self.output_dir + "/%s_summary_job%d.json" % (self.name, self.idx)
        self.condor_executable_file = self.log_output_dir + "/%s_executable_job%d.sh" % (self.name, self.idx)
        self.condor_submit_file = self.log_output_dir + "/%s_batch_submit_job%d.txt" % (self.name, self.idx)
        

    def write_config(self):
        """
        Write config file to json. 
        """
        self.config["name"] = self.name_full
        self.config["job_id"] = self.idx
        self.config["dir"] = self.dir
        self.config["output_dir"] = self.output_dir
        self.config["output_name"] = self.output_name
        self.config["summary_file"] = self.summary_file
        self.config["files"] = [file.name for file in self.inputs]

        if os.path.exists(self.config_file):
            logger.warning("[Job : write_config] Overwriting existing config file '%s'." % (self.config_file))
        with open(self.config_file, "w") as f_out:
            json.dump(self.config, f_out, sort_keys = True, indent = 4)
        self.wrote_config = True


    def query_status(self):
        """
        Update the job's status
        """
        if self.status in ["waiting", "failed", "completed", "retired"]:
            pass
        # Otherwise, the job may have previously been idle or running and has since started running/completed
        elif os.path.exists(self.summary_file): # a job writes its summary json as the very last step, so its presence implies successful completion
            self.status = "completed"
        else:
            self.monitor() 


    def write_condor_files(self):
        """

        """
        pass


    def write_python_executable(self):
        lines = []
        lines.append("import json")
        lines.append("import os")
        lines.append("from higgs_dna.utils.logger_utils import setup_logger")
        lines.append("from higgs_dna.analysis import run_analysis")
        lines.append("")
        lines.append("logger = setup_logger('DEBUG')")
        lines.append("config_file = '%s'" % self.config_file) # start with absolute path to config file
        lines.append("if not os.path.exists(config_file):") # in case this is a remote job, absolute path will not exist
        lines.append("    config_file = os.path.split(config_file)[-1]") # and we should have sent a copy that is present in the job's current dir
        lines.append("with open(config_file, 'r') as f_in:")
        lines.append("    config = json.load(f_in)")
        lines.append("")
        lines.append("run_analysis(config)") # FIXME: not compatible if another function is specified

        if os.path.exists(self.python_executable_file):
            logger.warning("[Job : write_python_executable] Overwriting existing python executable: %s" % (self.python_executable_file))

        with open(self.python_executable_file, "w") as f_out:
            for line in lines:
                f_out.write(line + "\n")
        self.wrote_python_executable = True
 

    def monitor(self):
        """

        """
        pass


    def kill(self):
        """

        """
        pass


    def submit(self, dry_run = False, reconfigure = False):
        """
        Write all associated job files:
            - python config
            - python executable
            - condor executable
            - condor submit
        if they have not already been written or if the file does not already exist (user may have manually deleted them),
        and then submit the job.

        :param dry_run: if True, don't actually submit the job
        :type dry_run: bool, defaults to False
        :param reconfigure: if True, force rewriting of associated job files
        :type reconfigure: bool, defaults to False
        :returns: whether the job was submitted (or is able to be submitted in the case of ``dry_run = True``)
        :rtype: bool
        """
        if reconfigure or (not self.wrote_config or not os.path.exists(self.config_file)):
            self.write_config()

        if reconfigure or (not self.wrote_python_executable or not os.path.exists(self.python_executable_file)):
            self.write_python_executable()

        if reconfigure or (not self.wrote_condor_files or not os.path.exists(self.condor_executable_file) or not os.path.exists(self.condor_submit_file)):
            self.write_condor_files()

        if self.status == "retired":
            return False

        if self.n_attempts >= 5 or self.force_retirement:
            logger.info("[Job : submit] Job '%s_%d' has been submitted %d times, retiring job. Jobs can be unretired with run_analysis.py through the `--unretire_jobs` option." % (self.name, self.idx, self.n_attempts))
            self.status = "retired"
            return False

        if dry_run:
            return True

        self.submit_to_batch()

        self.n_attempts += 1
        if self.n_attempts == 1 and isinstance(self, CondorJob):
            self.status = "idle"
            return True # don't waste time monitoring each individual job on the first try, they will be idle 
        self.monitor()
        logger.debug("[Job : submit] Job '%s_%d' has status '%s'." % (self.name, self.idx, self.status))
        return True


    def submit_to_batch(self):
        """

        """
        pass


class LocalJob(Job):
    """

    """
    def submit_to_batch(self):
        with open(self.log_file, "w") as f:
            self.p = subprocess.Popen(
                    "python %s" % (self.python_executable_file),
                    shell = True,
                    stdout = f,
                    stderr = f
            )


    def monitor(self):
        if self.p.poll() is not None:
            if os.path.exists(self.summary_file):
                self.status = "completed"
            else:
                self.status = "failed"
        else:
            self.status = "running"


    def kill(self):
        self.p.kill()
        self.status = "waiting"


class CondorJob(Job):
    """

    """
    REQUESTS = {
            "REQ_MEMORY" : 2048, # request 2GB of memory
            "REQ_DISK" : 5000, # request ~5GB of disk
            "REQ_NCPUS" : 1 # just 1 CPU
    }

    def write_condor_files(self):
        """
        Prepare condor executable and condor submit files based on 
        """
        # Set files for job according to whether the job has access to the host filesystem (remote_job == True means do not have access)
        if self.host_params["remote_job"]:
            self.job_config_file = os.path.split(self.config_file)[-1]
            self.job_python_executable_file = os.path.split(self.python_executable_file)[-1]
            self.job_summary_file = os.path.split(self.summary_file)[-1]
        else:
            self.job_config_file = self.config_file
            self.job_python_executable_file = self.python_executable_file
            self.job_summary_file = self.summary_file

        # Get proper templates according to host
        self.hdna_base = get_HiggsDNA_base()
        self.hdna_conda = get_HiggsDNA_conda()

        self.condor_exe_template = self.hdna_base + "/higgs_dna/job_management/condor/%s/exe_template.sh" % self.host
        self.condor_sub_template = self.hdna_base + "/higgs_dna/job_management/condor/%s/submit_template.txt" % self.host
        
        self.write_condor_executable_file()
        self.write_condor_submit_file()
        self.wrote_condor_files = True


    def update_file(self, old, new, replacement_map):
        """
        Reads in the lines from a file ``old``, takes each ``key`` and ``value`` from replacement map and replaces each instance of ``key`` with ``value``, and writes the resulting file out to ``new``.

        :param old: file to use as starting point
        :type old: str
        :param new: path to write resulting file to
        :type new: str
        :param replacement_map: dictionary of str : str specifying the replacements to make
        :type replacement_map: dict
        """

        with open(old, "r") as f_in:
            lines = f_in.readlines()

        lines = "".join(lines)
        for k,v in replacement_map.items():
            lines = lines.replace(k, str(v)) 

        with open(new, "w") as f_out:
            f_out.write(lines)


    def write_condor_executable_file(self):
        """
        Set placeholders in condor executable file to proper values for job.
        """
        replacement_map = {}
        replacement_map["PYTHON_FILE"] = self.job_python_executable_file
        replacement_map["SUMMARY_FILE"] = self.job_summary_file
        if self.host_params["needs_copy_proxy"]:
            replacement_map["GRID_PROXY"] = os.path.split(self.proxy)[-1]

        # update xrdcp placeholders for copying tar files into job
        if self.host_params["needs_tar"] and "xrd_redirector" in self.host_params.keys():
            replacement_map["XRD_CONDA_TARFILE"] = self.xrd_conda_tarfile
            replacement_map["XRD_ANALYSIS_TARFILE"] = self.xrd_analysis_tarfile

        # update gfal-copy placeholders
        if "gfal_redirector" in self.host_params.keys():
            to_replace, replace_with = self.host_params["gfal_redirector"]
            replacement_map["GFAL_BATCH_OUTPUT_DIR"] = self.output_dir.replace(to_replace, replace_with)    

        # if not a remote job, that means we are not sending a tar of the conda env and setting it up in the node, and we need to update the python path to point to our HiggsDNA version
        #if not self.host_params["remote_job"]:
        #    replacement_map["python"] = "%s/bin/python" % (self.hdna_conda)

        self.update_file(
                old = self.condor_exe_template,
                new = self.condor_executable_file,
                replacement_map = replacement_map
        )


    def write_condor_submit_file(self):
        """
        Set placeholders in condor submit file to proper values for job.
        """
        replacement_map = {}
        replacement_map["EXECUTABLE"] = self.condor_executable_file
        replacement_map["PYTHON_FILE"] = self.python_executable_file
        replacement_map["CONFIG_FILE"] = self.config_file
        replacement_map["HIGGS_DNA_BASE"] = self.hdna_base 

        # For hosts that need to copy proxy to home area, we pass that version
        if self.host_params["needs_copy_proxy"]:
            replacement_map["GRID_PROXY"] = self.hdna_base + "/" + os.path.split(self.proxy)[-1]
        else:
            replacement_map["GRID_PROXY"] = self.proxy # otherwise give the full path

        if "sites" in self.host_params.keys():
            replacement_map["TARGET_SITES"] = self.host_params["sites"]

        replacement_map["OUTPUT"] = self.log_output_dir + "/$(Cluster).$(Process).out" 
        replacement_map["ERROR"] = self.log_output_dir + "/$(Cluster).$(Process).err" 
        replacement_map["LOG"] = self.log_output_dir + "/$(Cluster).$(Process).log" 
        replacement_map["BATCH_NAME"] = self.name

        # memory, disk, cpu
        for k,v in self.REQUESTS.items():
            replacement_map[k] = v

        self.update_file(
                old = self.condor_sub_template,
                new = self.condor_submit_file,
                replacement_map = replacement_map
        )


    def submit_to_batch(self):
        if not self.can_submit:
            logger.info("[CondorJob : submit_to_batch] Job '%s_%d' with status '%s' can't be submitted to condor." % (self.name, self.idx, self.status)) 
            return

        result = do_cmd("condor_submit %s" % self.condor_submit_file)

        if "job(s) submitted to cluster" in result:
            success = True
            self.cluster_id = result.split("submitted to cluster ")[-1].split(".",1)[0].strip()
            logger.debug("[CondorJob : submit_to_batch] Submitted job '%s_%d' to cluster with id %s" % (self.name, self.idx, self.cluster_id))

        else:
            logger.exception("[CondorJob : submit_to_batch] Couldn't submit job to cluster because:\n----\n{0}\n----".format(result))
            raise RuntimeError()


    def monitor(self):
        if os.path.exists(self.summary_file):
            self.status = "completed"

        elif self.cluster_id is not None:
            cmd = "condor_q %s --json" % self.cluster_id
            result = do_cmd(cmd)

            if not result: # summary file does not exist and job is not found
                self.status = "failed"
            else:
                condor_status = json.loads(result)[0]["JobStatus"]
                self.status = CONDOR_STATUS_FLAGS[condor_status]
                if self.status == "submission_error":
                    self.can_submit = False


    def kill(self):
        self.monitor()
        if self.status == "idle" or self.status == "running" or self.status == "held":
            cmd = "condor_rm %s" % self.cluster_id
            do_cmd(cmd)
            self.status = "waiting"

        else:
            pass # if it's not idle, running, or held there's nothing to do
