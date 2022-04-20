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
    "bdt_file" : "/home/users/iareed/HggAnalysisDev/MVAs/output/Flash_gg_sync.xgb",
    "bdt_features" : [
        ("Diphoton", "eta"), ("Diphoton", "pt_mgg"), ("Diphoton", "dR"),
        ("LeadPhoton", "mvaID"), ("SubleadPhoton", "mvaID"), ("LeadPhoton", "eta"), ("SubleadPhoton", "eta"), ("LeadPhoton", "pixelSeed"), ("SubleadPhoton", "pixelSeed"), ("LeadPhoton", "pt_mgg"), ("SubleadPhoton", "pt_mgg"),
        "MET_pt",
        "n_jets", "n_leptons", "n_electrons", "n_muons", "n_taus",
        "jet_1_pt", "jet_1_eta", "jet_1_btagDeepFlavB",
        "jet_2_pt", "jet_2_eta", "jet_2_btagDeepFlavB",
        "jet_3_pt", "jet_3_eta", "jet_3_btagDeepFlavB",
        "jet_4_pt", "jet_4_eta", "jet_4_btagDeepFlavB",
        "jet_5_pt", "jet_5_eta", "jet_5_btagDeepFlavB",
        "lepton_1_pt", "lepton_1_eta", "lepton_1_mass", "lepton_1_charge", "lepton_1_id",
        "lepton_2_pt", "lepton_2_eta", "lepton_2_mass", "lepton_2_charge", "lepton_2_id" 
    ],
    "bdt_cuts" : [0.9883, 0.969307]
}

class TTHHSRTagger(Tagger):
    """
    Signal region tagger for the non-resonant HH->ggTauTau analysis.
    """
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(TTHHSRTagger, self).__init__(name, options, is_data, year)

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
