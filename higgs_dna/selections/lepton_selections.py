# TODO items
# 1. Make the selection more configurable
# 2. Should this (and selections more generally) be implemented as a class, where a Tagger will own multiple Selection instances?
#   - this would allow for TagSequence to more efficiently manage computation and avoid repeated computation (e.g. multiple taggers using the same jet selection)

import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_ELECTRONS = {
        "pt" : 10.0,
        "eta" : 2.4,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "WP90",
        "dr_photons" : 0.2
}

def select_electrons(electrons, options, clean, name = "none", tagger = None):
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_ELECTRONS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name

    standard_cuts = object_selections.select_objects(electrons, options, clean, name, tagger)

    if options["id"] == "WP90":
        id_cut = (electrons.mvaFall17V2Iso_WP90 == True) | ((electrons.mvaFall17V2noIso_WP90 == True) & (electrons.pfRelIso03_all < 0.3)) 
    elif not options["id"] or options["id"].lower() == "none":
        id_cut = electrons.pt > 0.
    else:
        logger.warning("[select_electrons] : Tagger '%s', id cut '%s' not recognized, not applying an ID cut." % (str(tagger), options["id"]))
        id_cut = electrons.pt > 0. 

    all_cuts = standard_cuts & id_cut

    if tagger is not None:
        tagger.register_cuts(
                names = ["standard object cuts", "id cut", "all cuts"],
                results = [standard_cuts, id_cut, all_cuts],
                cut_type = name
        )

    return all_cuts


DEFAULT_MUONS = {
        "pt" : 5.0,
        "eta" : 2.5,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "medium",       
        "pfRelIso03_all" : 0.3,
        "dr_photons" : 0.2,
        "global" : True
}

def select_muons(muons, options, clean, name = "none", tagger = None):
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_MUONS,
        new = options
    )
    
    tagger_name = "none" if tagger is None else tagger.name

    standard_cuts = object_selections.select_objects(muons, options, clean, name, tagger)

    if options["id"] == "medium":
        id_cut = muons.mediumId == True
    elif not options["id"] or options["id"].lower() == "none":
        id_cut = muons.pt > 0.
    else:
        logger.warning("[select_muons] : Tagger '%s', id cut '%s' not recognized, not applying an ID cut." % (str(tagger), options["id"]))
        id_cut = muons.pt > 0.

    if options["global"]:
        global_cut = muons.isGlobal == True
    else:
        global_cut = muons.pt > 0

    all_cuts = standard_cuts & id_cut & global_cut

    if tagger is not None:
        tagger.register_cuts(
                names = ["standard object cuts", "id cut", "global_muon cut", "all cuts"],
                results = [standard_cuts, id_cut, global_cut, all_cuts],
                cut_type = name
        )

    return all_cuts
