import awkward
import vector
import xgboost
import numpy

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.utils import awkward_utils, misc_utils

DEFAULT_OPTIONS = {
    "bdt_file" : "/home/users/hmei/myWorkspace/HiggsDNA/metadata/BDTs/20UL_30Mar2022_fixIsoTrk.xgb", # if running on condor, this file needs to be placed somewhere under higgs_dna/ so that it is included in the tar file. We probably want to think of a better long term solution for this.
    "bdt_features" : [
        "n_electrons", "n_muons", "n_taus", "n_isoTrks", "n_jets", "n_bjets",
        "MET_pt", "MET_gg_dPhi", "MET_ll_dPhi", "dPhi_MET_l", "lep12_dphi", "lep12_deta_bdt", "lep12_dr",
        "g1_ptmgg", "g1_eta_bdt", "g1_idmva", "g1_pixVeto", "g2_ptmgg", "g2_eta_bdt", "g2_idmva", "g2_pixVeto", "max_g_ptmgg", "min_g_ptmgg", "max_g_idmva", "min_g_idmva",
        "gg_ptmgg", "gg_eta", "gg_dR", "gg_dPhi", "gg_hel", "gg_tt_CS", "gg_tt_hel", "tt_hel",
        "tau_candidate_1_pt", "lep1_eta_bdt", "lep1_tightID", "tau_candidate_2_pt", "lep2_eta_bdt", "lep2_tightID", "max_lep_pt", "min_lep_pt",
        "Category", "jet_1_pt", "jet1_eta_bdt", "jet1_bTag", "jet_2_pt", "jet2_eta_bdt", "jet2_bTag", "max_bTag",
        "pt_tautau_SVFit", "eta_tautau_SVFit_bdt", "m_tautau_SVFit", "dR_tautau_SVFit", "dR_ggtautau_SVFit", "dPhi_tautau_SVFit", "dPhi_ggtautau_SVFit", "m_tautau_vis", "pt_tautau_vis", "eta_tautau_vis_bdt",
        "mX","m_llg_lead", "m_llg_subl"
    ],
    "bdt_cuts" : [0.9910, 0.972815]
}

class HHggTauTauNonResSRTagger(Tagger):
    """
    Signal region tagger for the non-resonant HH->ggTauTau analysis.
    """
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(HHggTauTauNonResSRTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, events):
        #####################################
        ### HH->ggTauTau Non-resonant SRs ###
        #####################################

        # Initialize BDT 
        bdt = xgboost.Booster()
        bdt.load_model(self.options["bdt_file"])

        # Convert events to proper format for xgb
        events_bdt = awkward.values_astype(events, numpy.float64)
        features_bdt = awkward.to_numpy(
                events_bdt[self.options["bdt_features"]]
        )
        features_bdt = xgboost.DMatrix(
                features_bdt.view((float, len(features_bdt.dtype.names)))
        )

        # Calculate BDT score
        events["bdt_score"] = bdt.predict(features_bdt)
        

        # Calculate SR cuts
        n_signal_regions = len(self.options["bdt_cuts"])
        sr_cuts = []
        for i in range(n_signal_regions):
            cut_sr = events.bdt_score >= self.options["bdt_cuts"][i]
            for j in range(len(sr_cuts)):
                cut_sr = cut_sr & ~(sr_cuts[j])

            # record which SR each event enters
            events["pass_sr_%d" % i] = cut_sr
            sr_cuts.append(cut_sr)


        # Calculate OR of all BDT cuts
        presel_cut = events.run < 0 # dummy all False
        for cut in sr_cuts:
            presel_cut = presel_cut | cut

        return presel_cut, events
