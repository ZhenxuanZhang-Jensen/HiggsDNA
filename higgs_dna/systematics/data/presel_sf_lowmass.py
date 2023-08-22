#DERIVED with hgg TnP tools for UL data 

#            2016 APV     2016         2017         2018
#EB low R9   0.999+-0.020 1.001+-0.020 1.007+-0.019 0.983+-0.025
#EB high R9  1.007+-0.023 1.006+-0.019 1.017+-0.016 0.998+-0.024
#EE high R9  1.015+-0.010 0.980+-0.015 1.037+-0.007  0.957+-0.015



PHOTON_PRESEL_SF_2016preVFP_LM = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 0.999,
            "uncertainty" : 0.020
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 1.007,
            "uncertainty" : 0.023
        },
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 1.015,
            "uncertainty" : 0.010
        }
    ]
}        

PHOTON_PRESEL_SF_2016postVFP_LM = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 1.001,
            "uncertainty" : 0.020
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 1.006,
            "uncertainty" : 0.019
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 0.980,
            "uncertainty" : 0.015
        }
    ]
}

PHOTON_PRESEL_SF_2017_LM = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 1.007,
            "uncertainty" : 0.019
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 1.017,
            "uncertainty" : 0.016
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 1.037,
            "uncertainty" : 0.007
        }
    ]
}

PHOTON_PRESEL_SF_2018_LM = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 0.983,
            "uncertainty" : 0.025
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 0.998,
            "uncertainty" : 0.024
        },  
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.90, 999.],
            "value" : 0.957,
            "uncertainty" : 0.015
        }
    ]
}


