import awkward
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.selections import lepton_selections, jet_selections
from higgs_dna.utils import awkward_utils, misc_utils

DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {
    "electrons" : {
        "pt" : 20.0,
        "dr_photons" : 0.2
    },
    "muons" : {
        "pt" : 20.0,
        "dr_photons" : 0.2
    },
    "jets" : {
        "pt" : 25.0,
        "eta" : 2.4,
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
    }
}   

class TTHPreselTagger(Tagger):
    """
    ttH Preselection tagger for tutorial
    """
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(TTHPreselTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )

    def calculate_selection(self, syst_tag, syst_events):
        # Electrons
        electron_cut = lepton_selections.select_electrons(
                electrons = syst_events.Electron,
                options = self.options["electrons"],
                clean = {
                    "photons" : {
                        "objects" : syst_events.Diphoton.Photon,
                        "min_dr" : self.options["electrons"]["dr_photons"]
                    }
                },
                name = "SelectedElectron",
                tagger = self
        )

        electrons = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedElectron",
                data = syst_events.Electron[electron_cut]
        )

        # Muons
        muon_cut = lepton_selections.select_muons(
                muons = syst_events.Muon,
                options = self.options["muons"],
                clean = {
                    "photons" : {
                        "objects" : syst_events.Diphoton.Photon,
                        "min_dr" : self.options["muons"]["dr_photons"]
                    }
                },
                name = "SelectedMuon",
                tagger = self
        )

        muons = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedMuon",
                data = syst_events.Muon[muon_cut]
        )

        # Jets
        jet_cut = jet_selections.select_jets(
                jets = syst_events.Jet,
                options = self.options["jets"],
                clean = {
                    "photons" : {
                        "objects" : syst_events.Diphoton.Photon,
                        "min_dr" : self.options["jets"]["dr_photons"]
                    },
                    "electrons" : {
                        "objects" : syst_events.SelectedElectron,
                        "min_dr" : self.options["jets"]["dr_electrons"]
                    },
                    "muons" : {
                        "objects" : syst_events.SelectedMuon,
                        "min_dr" : self.options["jets"]["dr_muons"]
                    }
                },
                name = "SelectedJet",
                tagger = self
        )

        jets = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedJet",
                data = syst_events.Jet[jet_cut]
        )

        bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis = 1, ascending = False)]

        # Preselection
        n_electrons = awkward.num(electrons)
        n_muons = awkward.num(muons)
        n_leptons = n_electrons + n_muons

        n_jets = awkward.num(jets)

        # Hadronic presel
        hadronic = (n_leptons == 0) & (n_jets >= 4)

        # Leptonic presel
        leptonic = (n_leptons >= 1) & (n_jets >= 2)

        presel_cut = hadronic | leptonic

        return presel_cut, syst_events
    
