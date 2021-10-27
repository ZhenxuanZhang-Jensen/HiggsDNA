CENTRAL_WEIGHT = "weight_central" # name of the central weight branch
NOMINAL_TAG = "nominal" # name of the nominal events (relevant when we have multiple sets of events corresponding to systematics with independent collections)

CONDOR_EXE_TEMPLATE = "higgs_dna/job_management/condor/executable_template.sh"
CONDOR_SUB_TEMPLATE = "higgs_dna/job_management/condor/condor_submit_template.txt"
CONDOR_EXE_TEMPLATE_LXPLUS = "higgs_dna/job_management/condor/executable_template_lxplus.sh"
CONDOR_SUB_TEMPLATE_LXPLUS = "higgs_dna/job_management/condor/condor_submit_template_lxplus.txt"

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

LUMI = {
    "2016" : 35.9,
    "2017" : 41.5,
    "2018" : 59.8
}

XRD_REDIRECTOR = {
    "UCSD" : "root://redirector.t2.ucsd.edu//"
}

GFAL_REDIRECTOR = {
    "UCSD" : "davs://redirector.t2.ucsd.edu:1094//"
}
