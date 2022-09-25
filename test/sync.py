import os
import awkward
import random
import math

from higgs_dna.utils.logger_utils import setup_logger

COMMAND = 'python scripts/run_analysis.py --config "metadata/analysis/sync_syst.json" --merge_outputs --log-level "DEBUG" --fpo 1 --output_dir "sync"' 

os.system("rm -rf sync/")
os.system(COMMAND)

"""
The numbers below are based on the following validations with flashgg:
    - Diphoton preselection : https://indico.cern.ch/event/1071721/contributions/4551056/attachments/2320292/3950844/HiggsDNA_DiphotonPreselectionAndSystematics_30Sep2021.pdf
    - Systematics with independent collections ("VARIATIONS"):
        - Photon FNUF : https://indico.cern.ch/event/1071724/contributions/4580497/attachments/2332453/3975163/HiggsDNA_SystematicsValidation_20Oct2021.pdf
        - Photon Material : https://indico.cern.ch/event/1071724/contributions/4580497/attachments/2332453/3975163/HiggsDNA_SystematicsValidation_20Oct2021.pdf
        - Photon MC Smearings : to be presented
        - Photon Scales : TODO
    - Weight systematics ("WEIGHTS")
        - Electron veto SF : https://indico.cern.ch/event/1071724/contributions/4580497/attachments/2332453/3975163/HiggsDNA_SystematicsValidation_20Oct2021.pdf
        - Trigger SF : https://indico.cern.ch/event/1071724/contributions/4580497/attachments/2332453/3975163/HiggsDNA_SystematicsValidation_20Oct2021.pdf
"""

VARIATIONS = {
        "nominal" : {
            "n_events" : 238529,
            "mgg_mean" : 124.95029237436138,
            "mgg_std"  : 6.446602002050745
        },
        "fnuf_down" : {
            "n_events" : 238346,
            "mgg_mean" : 124.747417836477,
            "mgg_std" : 6.456442578415462
        },
        "fnuf_up" : {
            "n_events" : 238687,
            "mgg_mean" : 125.15355480666447,
            "mgg_std" : 6.442464282196147
        },  
        "material_down" : {
            "n_events" : 238424,
            "mgg_mean" : 124.84880904043195, 
            "mgg_std" : 6.455341426771809
        },  
        "material_up" : {
            "n_events" : 238635,
            "mgg_mean" : 125.05260697411336,
            "mgg_std" : 6.447698396205765
        },  
        "mc_smear_down" : {
            "n_events" : 238528,
            "mgg_mean" : 124.95002846091435,
            "mgg_std" : 6.448206670393567
        },  
        "mc_smear_up" : {
            "n_events" : 238537,
            "mgg_mean" : 124.95024373679699, 
            "mgg_std" : 6.446183469760313
        } 
}

WEIGHTS = {
        "electron_veto_sf_Diphoton_Photon" : {
            "central" : {
                "mean" : 0.9776013870430849,
                "std" : 0.0069613946750999404
            },
            "up" : {
                "mean" : 0.9805363027137162,
                "std" : 0.005814417121197521
            },
            "down" : {
                "mean" : 0.9746713187914258,
                "std" : 0.00886875436052018
            }
        },
        "trigger_sf" : {
            "central" : {
                "mean" : 0.9743179157460937,
                "std" : 0.037648196360608656
            },
            "up" : {
                "mean" : 0.9772434641070897,
                "std" : 0.03677957246910504
            },
            "down" : {
                "mean" : 0.9713970837927464,
                "std" : 0.03858247778515959
            }
        }
} 


### Perform sync ###
logger = setup_logger("DEBUG")

# 1. Variations
pass_variations = True
for variation, sync_info in VARIATIONS.items():
    events = awkward.from_parquet("sync/merged_%s.parquet" % variation)

    result = {
        "n_events" : len(events),
        "mgg_mean" : awkward.mean(events.Diphoton_mass),
        "mgg_std"  : awkward.std(events.Diphoton_mass)
    }

    pass_current_variation = True
    for x,y in result.items():
        if not math.isclose(y, sync_info[x], rel_tol = 1e-06):
            logger.warning("[sync.py] Variation '%s' has target of %.3g for field '%s' but sync exercise found value of %.3g." % (variation, sync_info[x], x, y))
            pass_current_variation = False
            pass_variations = False

    if pass_current_variation:
        logger.info("[sync.py] 'SUCCESS' : variation '%s' has expected number of events & mean/std dev of mgg distribution." % variation) 

if pass_variations:
    logger.info("[sync.py] 'SUCCESS' : All systematics with independent collections are properly synced.") 

# 2. Weights
events = awkward.from_parquet("sync/merged_nominal.parquet")
pass_weights = True
for weight, sync_info in WEIGHTS.items():
    pass_current_weight = True
    for variation, target in sync_info.items():
        result = {
            "mean" : awkward.mean(events["weight_%s_%s" % (weight, variation)]),
            "std" : awkward.std(events["weight_%s_%s" % (weight, variation)])
        }
        for x,y in result.items():
            if not math.isclose(y, target[x], rel_tol = 1e-06):
                logger.warning("[sync.py] Weight '%s', variation '%s' has target of %.3g for field '%s' but sync exercise found value of %.3g." % (weight, variation, target[x], x, y))
                pass_current_weight = False
                pass_weights = False

    if pass_current_weight:
        logger.info("[sync.py] 'SUCCESS' : weight '%s' has expected mean and std dev for the following variations '%s'." % (weight, ",".join(list(sync_info.keys()))))

if pass_weights:
    logger.info("[sync.py] 'SUCCESS' : All weight systematics are properly synced.")
