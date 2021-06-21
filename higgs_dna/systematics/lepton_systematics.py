import awkward
import numpy

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import awkward_utils
from higgs_dna.systematics.utils import systematic_from_bins

DUMMY_LEPTON_SFs = {
    "variables" : ["lepton_pt", "lepton_eta"],
    "bins" : [
        {
            "lepton_pt" : [25.0, 50.0],
            "lepton_eta" : [0.0, 1.5],
            "value" : 0.5,
            "uncertainty" : 0.25
        },
        {
            "lepton_pt" : [25.0, 50.0],
            "lepton_eta" : [1.5, 999],
            "value" : 1.5,
            "uncertainty" : 0.25
        },
        {
            "lepton_pt" : [50.0, 999999.],
            "lepton_eta" : [0.0, 999999.],
            "value" : 0.5,
            "uncertainty" : 0.25
        },
 
    ]
}

def dummy_lepton_sf(events, central_only, input_collection):
    """
    Dummy function illustrating a lepton scale factor.
    Should be deleted once real examples are implemented.
    """
    required_fields = [
            (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    leptons = events[input_collection]
    variations = systematic_from_bins(
            bins = DUMMY_LEPTON_SFs,
            variables = {
                "lepton_pt" : leptons.pt,
                "lepton_eta" : leptons.eta
            },
            central_only = central_only
    )

    return variations 
