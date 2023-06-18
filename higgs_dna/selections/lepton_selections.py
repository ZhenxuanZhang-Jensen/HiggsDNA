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
        "id" : "WPL",
        "non_iso": None,
        "iso": None,
        "dr_photons" : 0.2,
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
    if options["id"] == "WPL":
        id_cut = (electrons.mvaFall17V2Iso_WPL == True)
        logger.debug("[select_electrons] : id_cut =WPLiso True")
    elif options["id"] == "WPL_noniso":
        id_cut = (electrons.mvaFall17V2noIso_WPL == True) & (electrons.mvaFall17V2Iso_WPL == False) 
        # id_cut = (electrons.mvaFall17V2noIso_WPL == True)
        logger.debug("[select_electrons] : id_cut = WPLiso false WPL noniso true")
    if options["iso"] is not None:
        iso_cut = electrons.miniPFRelIso_all <= options["iso"]
    elif options["iso"] is None:
        iso_cut = electrons.pt > 0.
    if options["non_iso"] is not None:
        non_iso_cut = electrons.miniPFRelIso_all > options["non_iso"]
    elif options["non_iso"] is None:
        non_iso_cut = electrons.pt > 0.
    if options["id"] == "WP90":
        # id_cut = (electrons.mvaFall17V2Iso_WP90 == True) | ((electrons.mvaFall17V2noIso_WP90 == True) & (electrons.pfRelIso03_all < 0.3)) 
        id_cut = (electrons.mvaFall17V2Iso_WP90 == True) 
    elif options["id"] == "WP80":
        id_cut = (electrons.mvaFall17V2Iso_WP80 == True)
    if options["iso"] is not None:
        iso_cut = electrons.miniPFRelIso_all < options["iso"]
        
    if options["non_iso"] is not None:
        non_iso_cut = electrons.miniPFRelIso_all > options["non_iso"]
        
    elif not options["id"] or options["id"].lower() == "none":
        id_cut = electrons.pt > 0.
    
    # else:
    #     logger.warning("[select_electrons] : Tagger '%s', id cut '%s' not recognized, not applying an ID cut." % (str(tagger), options["id"]))
    #     id_cut = electrons.pt > 0. 
    if options["veto_transition"]:
        transition_cut = (abs(electrons.eta) < 1.4442) | (abs(electrons.eta) > 1.566)

    # if options["Z_mass_veto"]:
        
    all_cuts = standard_cuts & id_cut & transition_cut & iso_cut & non_iso_cut

    if tagger is not None:
        tagger.register_cuts(
                names = ["standard object cuts", "id cut", "ee-eb transition", "all cuts"],
                results = [standard_cuts, id_cut, transition_cut, all_cuts],
                cut_type = name
        )

    return all_cuts


DEFAULT_MUONS = {
        "pt" : 10.0,
        "eta" : 2.5,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "medium",  
        "non_iso": None,  
        "iso": None, 
        "non_pfRelIso04_all":None,  
        # "pfRelIso03_all" : 0.3,
        "pfRelIso04_all" : None,
        # "pfRelIso04_all" : 0.15,
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
    logger.debug("[select_muons] : options = %s" % str(options))
    tagger_name = "none" if tagger is None else tagger.name

    standard_cuts = object_selections.select_objects(muons, options, clean, name, tagger)
    if options["iso"] is not None:
        # iso_cut = muons.pfRelIso04_all < options["iso"]
        iso_cut = muons.miniPFRelIso_all < options["iso"]
    else:
        iso_cut = muons.pt > 0.
    if options["non_iso"] is not None:
        # non_iso_cut = muons.pfRelIso04_all > options["non_iso"]
        non_iso_cut = muons.miniPFRelIso_all > options["non_iso"]
    else:
        non_iso_cut = muons.pt > 0.
    if options["id"] == "medium":
        id_cut = muons.mediumId == True
    if options["id"] == "loose":
        id_cut = muons.looseId == True
    if options["id"] == "tight":
        id_cut = muons.tightId == True
    elif not options["id"] or options["id"].lower() == "none":
        id_cut = muons.pt > 0.
    else:
        logger.warning("[select_muons] : Tagger '%s', id cut '%s' not recognized, not applying an ID cut." % (str(tagger), options["id"]))
        id_cut = muons.pt > 0.

    if options["global"]:
        global_cut = muons.isGlobal == True
    else:
        global_cut = muons.pt > 0

    all_cuts = standard_cuts & id_cut & global_cut &iso_cut & non_iso_cut
    if tagger is not None:
        tagger.register_cuts(
                names = ["standard object cuts", "id cut", "global_muon cut", "all cuts"],
                results = [standard_cuts, id_cut, global_cut, all_cuts],
                cut_type = name
        )

    return all_cuts
