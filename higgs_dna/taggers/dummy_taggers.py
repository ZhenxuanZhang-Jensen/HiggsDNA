import awkward
import vector

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.selections import lepton_selections, jet_selections
from higgs_dna.utils import awkward_utils

class TTHTagger(Tagger):
    """
    Dummy tth tagger
    """
    def calculate_selection(self, syst_tag, syst_events):
        # Require >= 1 lepton, >= 1 jets OR 0 leptons, >= 3 jets

        # Electrons
        electron_cut = lepton_selections.select_electrons(
                electrons = syst_events.Electron,
                diphotons = syst_events.Diphoton,
                options = { "pt" : 25. },
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
                diphotons = syst_events.Diphoton,
                options = { "pt" : 25. },
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
                diphotons = syst_events.Diphoton,
                muons = syst_events.SelectedMuon,
                electrons = syst_events.SelectedElectron,
                options = {},
                tagger = self
        )
                
        jets = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedJet",
                data = syst_events.Jet[jet_cut]
        )

        # Preselection
        n_leptons = awkward.num(electrons) + awkward.num(muons)
        n_jets = awkward.num(jets)

        dipho_pt_requirement = syst_events.Diphoton.pt > 100
        dipho_pt_cut = awkward.num(syst_events.Diphoton[dipho_pt_requirement]) >= 1

        tth_leptonic_presel = (n_leptons >= 1) & (n_jets >= 2)
        tth_hadronic_presel = (n_leptons == 0) & (n_jets >= 4)
        

        presel_cut = dipho_pt_cut & (tth_leptonic_presel | tth_hadronic_presel)
        self.register_cuts(
                names = ["diphoton pt", "tth leptonic", "tth hadronic", "tth preselection"],
                results = [dipho_pt_cut, tth_leptonic_presel, tth_hadronic_presel, presel_cut]
        )

        return presel_cut, syst_events

class THQTagger(Tagger):
    """
    Dummy thq tagger
    """
    def calculate_selection(self, syst_tag, syst_events):

        # Electrons
        electron_cut = lepton_selections.select_electrons(
                electrons = syst_events.Electron,
                diphotons = syst_events.Diphoton,
                options = { "pt" : 25. },
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
                diphotons = syst_events.Diphoton,
                options = { "pt" : 35. },
                tagger = self
        )
        muons = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedMuonTight",
                data = syst_events.Muon[muon_cut]
        )

        # Jets
        jet_cut = jet_selections.select_jets(
                jets = syst_events.Jet,
                diphotons = syst_events.Diphoton,
                muons = syst_events.SelectedMuonTight,
                electrons = syst_events.SelectedElectron,
                options = {},
                tagger = self
        )
                
        jets = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedJet",
                data = syst_events.Jet[jet_cut]
        )
 
        # Preselection
        n_leptons = awkward.num(electrons) + awkward.num(muons)
        n_jets = awkward.num(jets)
        
        thq_leptonic_presel = (n_leptons >= 1) & (n_jets >= 1)
        thq_hadronic_presel = (n_leptons == 0) & (n_jets >= 3)

        presel_cut = thq_leptonic_presel | thq_hadronic_presel
        self.register_cuts(
                names = ["thq_leptonic", "thq_hadronic", "thq_preselection"],
                results = [thq_leptonic_presel, thq_hadronic_presel, presel_cut]
        )

        return presel_cut, syst_events 
