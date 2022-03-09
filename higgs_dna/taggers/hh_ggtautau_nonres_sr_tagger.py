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
    "bdt_file" : "/path/to/some/bdt.xgb", # if running on condor, this file needs to be placed somewhere under higgs_dna/ so that it is included in the tar file. We probably want to think of a better long term solution for this.
    "bdt_features" : [
        ("Diphoton", "eta"),
        ["LeadPhoton", "eta"], ["LeadPhoton", "mvaID"],
        "ditau_pt", "ditau_eta", "ditau_phi", "ditau_mass", "ditau_dR",
    ],
    "bdt_cuts" : [0.99, 0.98]
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
            for j in range(sr_cuts):
                cut_sr = cut_sr & ~(sr_cuts[j])

            # record which SR each event enters
            events["pass_sr_%d" % i] = cut_sr
            sr_cuts.append(cut_sr)


        # Calculate OR of all BDT cuts
        presel_cut = events.run < 0 # dummy all False
        for cut in sr_cuts:
            presel_cut = presel_cut | cut

        return presel_cut, events
