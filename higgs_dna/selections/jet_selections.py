# TODO items
# 1. Make the selection more configurable
# 2. Should this (and selections more generally) be implemented as a class, where a Tagger will own multiple Selection instances?
#   - this would allow for TagSequence to more efficiently manage computation and avoid repeated computation (e.g. multiple taggers using the same jet selection)

import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_JETS = {
    "pt" : 25.0,
    "eta" : 2.4
}

def select_jets(jets, options, clean, name = "none", tagger = None):
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_JETS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name 

    standard_cuts = object_selections.select_objects(jets, options, clean, name, tagger)

    # TODO: jet ID

    all_cuts = standard_cuts

    if tagger is not None:
        tagger.register_cuts(
            names = ["all cuts"],
            results = [all_cuts],
            cut_type = name
        )

    return all_cuts 
