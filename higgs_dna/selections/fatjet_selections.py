import awkward

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_FATJETS = {
        "pt" : 150.,
        "eta" : 2.4
}

def select_fatjets(fatjets, options, clean, name = "none", tagger = None): 
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_FATJETS,
        new = options
    )

    standard_cuts = object_selections.select_objects(fatjets, options, clean, name, tagger)

    all_cuts = standard_cuts

    if tagger is not None:
        tagger.register_cuts(
            names = ["all cuts"],
            results = [all_cuts],
            cut_type = name
        )

    return all_cuts
