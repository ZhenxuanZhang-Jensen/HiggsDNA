# If a site is supported for HTCondor submission, the default templates for submission files will be located at:
#   condor_submit   : higgs_dna/job_management/condor/SITE/submit_template.txt
#   condor_exe      : higgs_dna/job_management/condor/SITE/exe_template.sh  

# TODO : write instructions on how you could add a new host
HOST_PARAMS = {
        "UCSD" : {
            "xrd_redirector" : ("/ceph/cms/", "/ceph/cms/"),
            "copy_tar" : "cp",
            "gfal_redirector" : ("/ceph/cms/", "davs://redirector.t2.ucsd.edu:1095//"),
            "condor_base_path" : "/ceph/cms/store/user/USERNAME/HiggsDNA/DIRNAME/",
            "sites" : "T2_US_UCSD",
            "needs_tar" : True,
            "needs_copy_proxy" : False,
            "remote_job" : True # whether the job is unable to directly access the host site when running on condor. This means the job must write its files to relative paths, not absolute ones
        },
        "lxplus" : {
            "needs_tar" : True,
            "copy_tar" : "xrd",
            "xrd_redirector" : ("/eos/", "root://eosuser.cern.ch//eos/"),
            "condor_base_path" : "/eos/user/USERNAME_INITIAL/USERNAME/HiggsDNA/DIRNAME/",
            "needs_copy_proxy" : True,
            "remote_job" : False
        }
}

# HTCondor mappings (http://pages.cs.wisc.edu/~adesmet/status.html)
CONDOR_STATUS_FLAGS = {
    0 : "unexpanded", # appropriate behavior for this status not currently implemented
    1 : "idle",
    2 : "running",
    3 : "removed",
    4 : "completed",
    5 : "held",
    6 : "submission_error"
}

CONDA_TARFILE = "higgs-dna.tar.gz" # default conda tarfile area (so it doesn't get remade every time you run jobs. If you do install a new python package, simply remove the tarfile and it will be remade)
