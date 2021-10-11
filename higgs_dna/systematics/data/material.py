# Copied from flashgg
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2016_cfi.py#L369-L380
MATERIAL_2016 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {   
            "photon_eta" : [0.0, 1.0],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.000455
        },
        {   
            "photon_eta" : [0.0, 1.0],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.000233
        },
        {   
            "photon_eta" : [1.0, 1.5],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.002089
        },
        {   
            "photon_eta" : [1.0, 1.5],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.002089
        },
        {   
            "photon_eta" : [1.5, 999.],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.001090
        },
        {   
            "photon_eta" : [1.5, 999.],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.002377
        },
    ]
}

# Copied from flashgg
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L377-L399
MATERIAL_2017 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {   
            "photon_eta" : [0.0, 1.0],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.000455
        },
        {
            "photon_eta" : [0.0, 1.0],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.000233
        },
        {
            "photon_eta" : [1.0, 1.5],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.002089
        },
        {
            "photon_eta" : [1.0, 1.5],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.002089
        },
        {
            "photon_eta" : [1.5, 999.],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.001090
        },
        {
            "photon_eta" : [1.5, 999.],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.002377
        },
    ]
}

# Copied from flashgg:
# https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2018_Legacy_cfi.py#L360-L382
MATERIAL_2018 = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {   
            "photon_eta" : [0.0, 1.0],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.000455
        },
        {
            "photon_eta" : [0.0, 1.0],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.000233
        },
        {
            "photon_eta" : [1.0, 1.5],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.002089
        },
        {
            "photon_eta" : [1.0, 1.5],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.002089
        },
        {
            "photon_eta" : [1.5, 999.],
            "photon_r9" : [0.0, 0.94],
            "value" : 1.0,
            "uncertainty" : 0.001090
        },
        {
            "photon_eta" : [1.5, 999.],
            "photon_r9" : [0.94, 999.],
            "value" : 1.0,
            "uncertainty" : 0.002377
        },
    ]
}
