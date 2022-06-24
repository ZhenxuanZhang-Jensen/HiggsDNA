CENTRAL_WEIGHT = "weight_central" # name of the central weight branch
NOMINAL_TAG = "nominal" # name of the nominal events (relevant when we have multiple sets of events corresponding to systematics with independent collections)

#https://indico.cern.ch/event/560226/contributions/2277448/attachments/1324704/1988050/wgm_vfp_change_ebutz.pdf
#pre-VFP  runs: 273150-278800 lumi: 19480.4566773 /pb
#post-VFP runs: 278801-284044 lumi: 16755.0362868 /pb

LUMI = {
    "2016" : 35.9,
    "2016UL_preVFP" : 19.48, # 2016 APV
    "2016UL_postVFP" : 16.76, # 2016
    "2017" : 41.5,
    "2018" : 59.8
}

GOLDEN_JSON = {
    "2016UL_preVFP" : "metadata/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2016UL_postVFP" : "metadata/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2017" : "metadata/golden_json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
    "2018" : "metadata/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
}

# nanoAOD branches to always include
BRANCHES = {
    "data" : {
        "2016" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
            "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
            "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55"
        ],
        "2016UL_postVFP" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
            "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
            "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55"
        ],
        "2016UL_preVFP" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
            "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
            "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55"
        ],
        "2017" : [
            "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
            "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"
        ],
        "2018" : [
            "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
            "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto"
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
