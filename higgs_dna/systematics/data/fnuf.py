# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2016_cfi.py#L37-L46 
FNUF_2016 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.00044
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.00156
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.00003
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.00165
        }
    ]
}

# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L73-L82
FNUF_2017 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.00,
            "uncertainty" : 0.00062
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.94, 999.],
            "value" : 1.00,
            "uncertainty" : 0.00208
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.00005
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.00227
        }
    ]
}

# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2018_Legacy_cfi.py#L72-L81
FNUF_2018 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.0007
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.0022
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.00005
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.00251
        }
    ]
}
