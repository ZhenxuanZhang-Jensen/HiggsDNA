# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2016_cfi.py#L25-L34 
PHOTON_ELECTRON_VETO_SF_2016 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 0.9973,
            "uncertainty" : 0.0024
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 0.9963,
            "uncertainty" : 0.0006
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.90],
            "value" : 0.9700,
            "uncertainty" : 0.0075
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 0.9909,
            "uncertainty" : 0.0018
        }
    ]
}

# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L26-L35
PHOTON_ELECTRON_VETO_SF_2017 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 0.9838,
            "uncertainty" : 0.0024
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 0.9913,
            "uncertainty" : 0.0009
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.90],
            "value" : 0.9777,
            "uncertainty" : 0.0180
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 0.9784,
            "uncertainty" : 0.0026
        }
    ]
}

# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2018_Legacy_cfi.py#L25-L34
PHOTON_ELECTRON_VETO_SF_2018 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 0.9819,
            "uncertainty" : 0.0020
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 0.9931,
            "uncertainty" : 0.0005
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.90],
            "value" : 0.9680,
            "uncertainty" : 0.0069
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 0.9785,
            "uncertainty" : 0.0018
        }
    ]
}

