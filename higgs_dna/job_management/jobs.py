import os
import json

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils.metis_utils import do_cmd

class Job():
    """

    """
    def __init__(self, name, inputs, function, config, output_dir, batch_output_dir, idx):
        self.name = name
        self.inputs = inputs
        self.function = function
        self.config = config
        self.config["remote_job"] = False

        self.idx = idx
        self.name_full = self.name + "_" + str(self.idx)
        self.output_dir = output_dir
        self.batch_output_dir = batch_output_dir
        self.batch_output_base_dir = "/".join(batch_output_dir.split("/")[:-2]) # remove last subdir to get the base dir
        self.job_batch_output_dir = self.batch_output_dir

        self.basepath = do_cmd("pwd").split("HiggsDNA")[0]

        self.summary_file = self.batch_output_dir + "/%s_summary_job%d.json" % (self.name, self.idx)
        self.config_file = self.output_dir + "/%s_config_job%d.json" % (self.name, self.idx)
        self.executable = self.output_dir + "/%s_executable_job%d.sh" % (self.name, self.idx)
        self.python_file = self.executable.replace(".sh", ".py")
        self.batch_submit_file = self.output_dir + "/%s_batch_submit_job%d.txt" % (self.name, self.idx)

        self.made_executable = False
        self.made_python_file = False
        self.made_batch_submit = False
        self.wrote_config = False

        self.cluster_id = None
        self.can_submit = True

        self.n_attempts = 0

        if os.path.exists(self.summary_file): # user may have ctrl+c-ed during the last run. This way we also prevent unintended overwriting of files.
            self.status = "completed"
        else:
            self.status = "waiting"
        

    def write_config(self):
        self.config["name"] = self.name_full
        self.config["job_id"] = self.idx
        self.config["output_dir"] = self.batch_output_dir
        self.config["job_batch_output_dir"] = self.job_batch_output_dir
        self.config["summary_file"] = self.summary_file
        self.config["files"] = [file.name for file in self.inputs]
        with open(self.config_file, "w") as f_out:
            json.dump(self.config, f_out, sort_keys = True, indent = 4)
        self.wrote_config = True


    def query_status(self):
        """

        """
        if self.status in ["waiting", "failed", "killed", "completed", "retired"]:
            pass
        # Otherwise, the job may have previously been idle or running and has since started running/completed
        elif os.path.exists(self.summary_file):
            self.status = "completed"

        else:
            self.monitor() 


    def create_executable(self):
        """

        """
        pass


    def create_python_file(self):
        lines = []
        lines.append("import json")
        lines.append("from higgs_dna.utils.logger_utils import setup_logger")
        lines.append("from higgs_dna.analysis import run_analysis")
        lines.append("")
        lines.append("logger = setup_logger('DEBUG')")
        lines.append("with open('%s', 'r') as f_in:" % self.config_file)
        lines.append("    config = json.load(f_in)")
        lines.append("")
        lines.append("run_analysis(config)")

        if os.path.exists(self.python_file):
            logger.warning("[LocalJob : create_executable] Overwriting existing executable: %s" % (self.executable))

        with open(self.python_file, "w") as f_out:
            for line in lines:
                f_out.write(line + "\n")

        self.made_python_file = True
 

    def monitor(self):
        """

        """
        pass

    def kill(self):
        """

        """
        pass

    def submit(self, dry_run = False):
        """

        """
        if not self.wrote_config:
            self.write_config()

        if not self.made_executable:
            self.create_executable()

        if self.status == "retired":
            return False

        if self.n_attempts >= 5:
            logger.info("[Job : submit] Job '%s_%d' has been submitted %d times, permanently retiring job." % (self.name, self.idx, self.n_attempts))
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


import subprocess
class LocalJob(Job):
    """

    """

    def create_executable(self):
        self.create_python_file()
        os.system("touch %s" % self.executable)
        self.made_executable = True


    def submit_to_batch(self):
        p = subprocess.Popen("python %s" % (self.python_file), shell = True)
        self.pid = p.pid

    def monitor(self):
        #p_status = self.p.poll()
        if not os.path.exists("/proc/%s" % self.pid):
            if os.path.exists(self.summary_file):
                self.status = "completed"
            else:
                self.status = "failed"
        else:
            self.status = "running"


    def kill(self):
        self.p.kill()
        self.status = "killed"


from higgs_dna.utils.misc_utils import expand_path
from higgs_dna.constants import CONDOR_EXE_TEMPLATE, CONDOR_SUB_TEMPLATE, CONDOR_STATUS_FLAGS, CONDOR_EXE_TEMPLATE_LXPLUS, CONDOR_SUB_TEMPLATE_LXPLUS, XRD_REDIRECTOR, GFAL_REDIRECTOR

class CondorJob(Job):
    """

    """

    def create_executable(self):
        self.create_condor_submit_file()
        self.create_python_file()
        
        if self.host == "lxplus":
            with open(expand_path(CONDOR_EXE_TEMPLATE_LXPLUS), "r") as f_in:
                lines = f_in.readlines()
            
            lines = "".join(lines)
            lines = lines.replace("OUTPUT_DIR", self.output_dir)
            lines = lines.replace("PYTHON_FILE", self.python_file)
            lines = lines.replace("HIGGS_DNA_PATH", do_cmd("echo $CONDA_PREFIX"))
            lines = lines.replace("HIGGS_DNA_BASE", do_cmd("pwd").split("HiggsDNA")[0] + "HiggsDNA/")
            lines = lines.replace("GRID_PROXY", self.proxy.split("/")[-1])

        else:
            self.job_config_file = self.config_file.split("/")[-1]
            self.job_summary_file = self.summary_file.split("/")[-1]
            self.job_python_file = self.python_file.split("/")[-1]

            if self.host in ["UCSD"]:
                self.job_batch_output_dir = self.batch_output_dir.replace("/hadoop/cms/", GFAL_REDIRECTOR[self.host])
            else:
                self.job_batch_output_dir = self.batch_output_dir

            # replace absolute path to config with relative path (since absolute path will no longer be relevant in job node)
            os.system("sed -i s@%s@%s@g %s" % (self.config_file, self.job_config_file, self.python_file)) 

            with open(expand_path(CONDOR_EXE_TEMPLATE), "r") as f_in:
                lines = f_in.readlines()

            lines = "".join(lines)
            lines = lines.replace("PYTHON_FILE", self.job_python_file)
            lines = lines.replace("CONFIG_FILE", self.job_config_file)
            lines = lines.replace("SUMMARY_FILE", self.job_summary_file)
            lines = lines.replace("BATCH_OUTPUT_DIR", self.job_batch_output_dir)
            lines = lines.replace("XRD_CONDA_TARFILE", self.batch_output_base_dir.replace("/hadoop/cms/", XRD_REDIRECTOR[self.host]) + "/" + self.conda_tarfile.split("/")[-1])
            lines = lines.replace("XRD_ANALYSIS_TARFILE", self.batch_output_base_dir.replace("/hadoop/cms/", XRD_REDIRECTOR[self.host]) + "/" + self.analysis_tarfile.split("/")[-1])

        with open(self.executable, "w") as f_out:
            f_out.write(lines)

        self.made_executable = True
        self.create_condor_submit_file()


    def create_condor_submit_file(self):
        if self.host == "lxplus":
            with open(expand_path(CONDOR_SUB_TEMPLATE_LXPLUS), "r") as f_in:
                lines = f_in.readlines()

        else:
            with open(expand_path(CONDOR_SUB_TEMPLATE), "r") as f_in:
                lines = f_in.readlines()

        lines = "".join(lines)
        lines = lines.replace("CONDA_TARFILE", self.conda_tarfile)
        lines = lines.replace("ANALYSIS_TARFILE", self.analysis_tarfile)
        lines = lines.replace("PYTHON_FILE", self.python_file)
        lines = lines.replace("CONFIG_FILE", self.config_file)
        lines = lines.replace("BASEPATH", self.basepath)
        lines = lines.replace("EXECUTABLE", self.executable)
        lines = lines.replace("OUTPUT", self.output_dir + "/$(Cluster).$(Process).out") 
        lines = lines.replace("ERROR", self.output_dir + "/$(Cluster).$(Process).err")
        lines = lines.replace("LOG", self.output_dir + "/$(Cluster).$(Process).log")
        lines = lines.replace("BATCH_NAME", self.name)
        if self.host == "lxplus":
            lines = lines.replace("GRID_PROXY", do_cmd("pwd").split("HiggsDNA")[0] + "HiggsDNA/" + self.proxy.split("/")[-1])
        else:
            lines = lines.replace("GRID_PROXY", self.proxy)


        with open(self.batch_submit_file, "w") as f_out:
            f_out.write(lines)

        self.made_batch_submit = True


    def submit_to_batch(self):
        if not self.can_submit:
            logger.info("[CondorJob : submit_to_batch] Job '%s_%d' with status '%s' can't be submitted to condor." % (self.name, self.idx, self.status)) 
            return

        result = do_cmd("condor_submit %s" % self.batch_submit_file)

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

            #logger.debug("[CondorJob : monitor] Job '%s_%d' was queried with command '%s', we found responses '%s' and obtained job status '%s'" % (self.name, self.idx, cmd, result, self.status))

    def kill(self):
        self.monitor()
        if self.status == "idle" or self.status == "running" or self.status == "held":
            cmd = "condor_rm %s" % self.cluster_id
            do_cmd(cmd)
            self.status = "killed"

        else:
            pass # if it's not idle, running, or held there's nothing to do
