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
    "bdt_file" : "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/POGsel_04May2022.xgb", # if running on condor, this file needs to be placed somewhere under higgs_dna/ so that it is included in the tar file. We probably want to think of a better long term solution for this.
    "bdt_features" : [
        "n_electrons", "n_muons", "n_taus", "n_iso_tracks", "n_jets", "n_bjets",
        "MET_pt", "diphoton_met_dPhi", "MET_ll_dPhi", "lead_lepton_met_dphi", "ditau_dphi", "ditau_deta", "ditau_dR",
        ("LeadPhoton", "pt_mgg"), ("LeadPhoton", "eta"), ("LeadPhoton", "mvaID"), ("LeadPhoton", "pixelSeed"), ("SubleadPhoton", "pt_mgg"), ("SubleadPhoton", "eta"), ("SubleadPhoton", "mvaID"), ("SubleadPhoton", "pixelSeed"), ("Diphoton", "max_mvaID"), ("Diphoton", "min_mvaID"),
        ("Diphoton", "pt_mgg"), ("Diphoton", "eta"), ("Diphoton", "dR"), ("Diphoton", "dPhi"), ("Diphoton", "helicity"), "gg_tt_CS", "gg_tt_hel", "tt_hel",
        "lead_lepton_pt", "lead_lepton_eta", "sublead_lepton_pt", "sublead_lepton_eta",
        "category", "jet_1_pt", "jet_1_eta", "jet_1_btagDeepFlavB", "jet_2_pt", "jet_2_eta", "jet_2_btagDeepFlavB", "b_jet_1_btagDeepFlavB",
        "pt_tautau_SVFit", "eta_tautau_SVFit_bdt", "m_tautau_SVFit", "dR_tautau_SVFit", "dR_ggtautau_SVFit", "dPhi_tautau_SVFit", "dPhi_ggtautau_SVFit", "ditau_mass", "ditau_pt", "ditau_eta",
        "mX","dilep_leadpho_mass", "dilep_subleadpho_mass"
    ],
    "bdt_cuts" : [0.9898, 0.882222, 0.0]
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

        bdt_features = []
        for x in self.options["bdt_features"]:
            if isinstance(x, tuple):
                name_flat = "_".join(x)
                events_bdt[name_flat] = events_bdt[x]
                bdt_features.append(name_flat)
            else:
                bdt_features.append(x)
 
            if bdt_features[-1] in ["LeadPhoton_eta", "SubleadPhoton_eta", "ditau_deta", "ditau_eta", "jet_1_eta", "jet_2_eta", "lead_lepton_eta", "sublead_lepton_eta"]:
                events_bdt[bdt_features[-1]] = awkward.where(
                        events.Diphoton.eta < 0,
                        events_bdt[bdt_features[-1]] * -1,
                        events_bdt[bdt_features[-1]]
                )

        features_bdt = awkward.to_numpy(
                events_bdt[bdt_features]
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
