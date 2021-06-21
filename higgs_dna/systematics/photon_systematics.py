import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import awkward_utils
from higgs_dna.systematics.utils import systematic_from_bins

PHOTON_PRESELECTION_SFs = {
    "variables" : ["photon_eta", "photon_r9"],
    "bins" : [
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.0, 0.85],
            "value" : 1.0057,
            "uncertainty" : 0.0010
        },
        {
            "photon_eta" : [0.0, 1.5],
            "photon_r9" : [0.85, 999.],
            "value" : 0.9988,
            "uncertainty" : 0.0009
        },
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.0, 0.9],
            "value" : 0.9443,
            "uncertainty" : 0.0072
        },
        {
            "photon_eta" : [1.5, 6.0],
            "photon_r9" : [0.9, 999.],
            "value" : 0.9947,
            "uncertainty" : 0.0051
        },

    ]
}

def photon_preselection_sf(events, central_only):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins = PHOTON_PRESELECTION_SFs,
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        central_only = central_only
    ) 
    
    return variations

def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations
