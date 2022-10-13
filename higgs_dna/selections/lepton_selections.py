import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_ELECTRONS = {
        "pt" : 10.0,
        "eta" : 2.47,
        "dxy" : 0.5,
        "dz" : 0.5,
        "id" : "WPL",
        "dr_photons" : 0.4,
        "dr_jets" : 0.4,
        "veto_transition" : True
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
    if options["id"] == "WPL":
        id_cut = (electrons.mvaFall17V2Iso_WPL == True) | ((electrons.mvaFall17V2Iso_WPL == True) & (electrons.pfRelIso03_all < 0.3))
    elif not options["id"] or options["id"].lower() == "none":
        id_cut = electrons.pt > 0.
    else:
        logger.warning("[select_electrons] : Tagger '%s', id cut '%s' not recognized, not applying an ID cut." % (str(tagger), options["id"]))
        id_cut = electrons.pt > 0. 
    if options["veto_transition"]:
        transition_cut = (abs(electrons.eta) < 1.37) | (abs(electrons.eta) > 1.52)
      #   transition_cut = (abs(electrons.eta) < 1.4442) | (abs(electrons.eta) > 1.566)


    all_cuts = standard_cuts & id_cut & transition_cut

    if tagger is not None:
        tagger.register_cuts(
                names = ["standard object cuts", "id cut", "ee-eb transition", "all cuts"],
                results = [standard_cuts, id_cut, transition_cut, all_cuts],
                cut_type = name
        )

    return all_cuts


DEFAULT_MUONS = {
        "pt" : 10.0,
        "eta" : 2.7,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "medium",       
        "pfRelIso03_all" : 0.3,
        "dr_photons" : 0.4,
        "dr_jets" : 0.4,
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

    if options["id"] == "loose":
        id_cut = muons.looseId == True
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
