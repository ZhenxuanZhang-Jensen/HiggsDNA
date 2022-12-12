import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_JETS = {
    "pt" : 25.0,
    "eta" : 2.4,
    "looseID" : True,
    "PUid" : False,
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
    if options["looseID"]:
        id_cut = jets.jetId >= 1 # jetID stored bitwise for loose/tight/tightLepVeto
    else:
        id_cut = jets.pt > 0

    if options["PUid"]:
        pu_cut =  ( (jets.puId > 0) & (jets.pt < 50) ) | ( jets.pt > 50 )# PUjetID stored bitwise for loose/tight/tightLepVeto
    else:
        pu_cut = jets.pt > 0

    all_cuts = standard_cuts & id_cut & pu_cut

    if tagger is not None:
        tagger.register_cuts(
            names = ["std cuts", "id cut", "all cuts", "pu_cut"],
            results = [standard_cuts, id_cut, "pu_cut", all_cuts],
            cut_type = name
        )

    return all_cuts
