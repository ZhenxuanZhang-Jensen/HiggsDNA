import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_TAUS = {
    "pt" : 18.0,
    "eta" : 2.3,
    "dz" : 0.2,
    "deep_tau_vs_ele" : 1,
    "deep_tau_vs_mu" : 2,
    "deep_tau_vs_jet" : 8,
    "dr_photons" : 0.2,
    "dr_electrons" : 0.2,
    "dr_muons" : 0.2,
    "decay_mode" : True
}

def select_taus(taus, options, clean, name = "none", tagger = None):
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_TAUS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name

    standard_cuts = object_selections.select_objects(taus, options, clean, name, tagger)

    deep_tau_vs_ele_cut = taus.idDeepTau2017v2p1VSe >= options["deep_tau_vs_ele"]
    deep_tau_vs_mu_cut = taus.idDeepTau2017v2p1VSmu >= options["deep_tau_vs_mu"]
    deep_tau_vs_jet_cut = taus.idDeepTau2017v2p1VSjet >= options["deep_tau_vs_jet"]
 
    if options["decay_mode"]:
        decay_mode_cut = taus.idDecayModeNewDMs == True
    else:
        decay_mode_cut = taus.pt > 0

    all_cuts = standard_cuts & deep_tau_vs_ele_cut & deep_tau_vs_mu_cut & deep_tau_vs_jet_cut & decay_mode_cut

    if tagger is not None:
        tagger.register_cuts(
            names = ["standard object cuts", "id_vs_ele", "id_vs_mu", "id_vs_jet", "decay mode cut", "all cuts"],
            results = [standard_cuts, deep_tau_vs_ele_cut, deep_tau_vs_mu_cut, deep_tau_vs_jet_cut, decay_mode_cut, all_cuts],
            cut_type = name
        )

    return all_cuts


DEFAULT_ISO_TRACKS = {
    "pt" : 5.,
    "eta" : 2.5,
    "isPFCandidate" : True,
    "fromPV" : True,
    "non_lepton" : True
}

def select_iso_tracks(iso_tracks, options, clean, name = "none", tagger = None):
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_ISO_TRACKS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name

    standard_cuts = object_selections.select_objects(iso_tracks, options, clean, name, tagger)

    if options["isPFCandidate"]:
        pf_cand_cut = iso_tracks.isPFcand == True
    else:
        pf_cand_cut = iso_tracks.pt > 0

    if options["fromPV"]:
        pv_cut = iso_tracks.fromPV == True
    else:
        pv_cut = iso_tracks.pt > 0

    if options["non_lepton"]:
        non_lep_cut = (abs(iso_tracks.pdgId) != 11) & (abs(iso_tracks.pdgId) != 13) 
    else:
        non_lep_cut = iso_tracks.pt > 0

    all_cuts = standard_cuts & pf_cand_cut & pv_cut & non_lep_cut 

    if tagger is not None:
        tagger.register_cuts(
            names = ["standard object cuts", "pf cand cut", "pv cut", "non lep cut", "all"], 
            results = [standard_cuts, pf_cand_cut, pv_cut, non_lep_cut, all_cuts],
            cut_type = name
        )

    return all_cuts
  
