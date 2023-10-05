# If a site is supported for HTCondor submission, the default templates for submission files will be located at:
#   condor_submit   : higgs_dna/job_management/condor/SITE/submit_template.txt
#   condor_exe      : higgs_dna/job_management/condor/SITE/exe_template.sh  

# TODO : write instructions on how you could add a new host
HOST_PARAMS = {
        "UCSD" : {
            "xrd_redirector" : ("/hadoop/cms/", "root://redirector.t2.ucsd.edu//"),
            "gfal_redirector" : ("/hadoop/cms/", "davs://redirector.t2.ucsd.edu:1094//"),
            "condor_base_path" : "/hadoop/cms/store/user/USERNAME/HiggsDNA/DIRNAME/",
            "sites" : "T2_US_UCSD",
            "needs_tar" : True,
            "needs_copy_proxy" : False,
            "remote_job" : True # whether the job is unable to directly access the host site when running on condor. This means the job must write its files to relative paths, not absolute ones
        },
        "lxplus" : {
            "needs_tar" : False,
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
