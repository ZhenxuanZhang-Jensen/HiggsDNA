CENTRAL_WEIGHT = "weight_central" # name of the central weight branch
NOMINAL_TAG = "nominal" # name of the nominal events (relevant when we have multiple sets of events corresponding to systematics with independent collections)

#https://indico.cern.ch/event/560226/contributions/2277448/attachments/1324704/1988050/wgm_vfp_change_ebutz.pdf
#pre-VFP  runs: 273150-278800 lumi: 19480.4566773 /pb
#post-VFP runs: 278801-284044 lumi: 16755.0362868 /pb

LUMI = {
    "2016" : 35.9,
    "2016UL_preVFP" : 19.48, 
    "2016UL_postVFP" : 16.76, 
    "2017" : 41.5,
    "2018" : 59.8
}

# nanoAOD branches to always include
BRANCHES = {
    "data" : {
        "2016" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        ],
        "2016UL_postVFP" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        ],
        "2016UL_preVFP" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        ],
        "2017" : [
            "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        ],
        "2018" : [
            "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        ],
        "any" : ["event", "run", "luminosityBlock"]
    },
    "mc" : {
        "2016APV" : [],
        "2016" : [],
        "2016UL_preVFP" : [],
        "2016UL_postVFP" : [],
        "2017" : [],
        "2018" : [],
        "any" : ["event", "run", "luminosityBlock"]
    }
}
