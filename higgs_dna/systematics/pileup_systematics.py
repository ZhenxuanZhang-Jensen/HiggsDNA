import awkward
import numpy

from correctionlib import _core

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils

PU_REWEIGHTING_FILE = {
    "2016UL_preVFP" : "jsonpog-integration/POG/LUM/2016preVFP_UL/puWeights.json",
    "2016UL_postVFP" : "jsonpog-integration/POG/LUM/2016postVFP_UL/puWeights.json",
    "2017" : "jsonpog-integration/POG/LUM/2017_UL/puWeights.json",
    "2018" : "jsonpog-integration/POG/LUM/2018_UL/puWeights.json"
}

PU_CAMPAIGN = {
    "2016UL_preVFP" : "Collisions16_UltraLegacy_goldenJSON",
    "2016UL_postVFP" : "Collisions16_UltraLegacy_goldenJSON",
    "2017" : "Collisions17_UltraLegacy_goldenJSON",
    "2018" : "Collisions18_UltraLegacy_goldenJSON"
}

def pu_reweight_sf(events, year, central_only):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/LUMI_puWeights_Run2_UL/LUMI_puWeights_2017_UL.html
    """

    required_fields = ["Pileup_nTrueInt"]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PU_REWEIGHTING_FILE[year]))

    num_true_interactions = awkward.to_numpy(events.Pileup_nTrueInt).astype(float)

    # Calculate SF and syst
    variations = {}

    variations["central"] = awkward.from_numpy(evaluator[PU_CAMPAIGN[year]].evalv(
            num_true_interactions,
            "nominal"
    ))

    if not central_only:
        for var in ["up", "down"]:
            variations[var] = awkward.from_numpy(evaluator[PU_CAMPAIGN[year]].evalv(
                    num_true_interactions,
                    var
            ))

    return variations
