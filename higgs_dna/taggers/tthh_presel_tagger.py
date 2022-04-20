import awkward
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.selections import object_selections, lepton_selections, jet_selections, tau_selections
from higgs_dna.utils import awkward_utils, misc_utils

DUMMY_VALUE = -999.
DEFAULT_OPTIONS = {
    "electrons" : {
        "pt" : 7.0,
        "eta" : 2.5,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "WP90",
        "dr_photons" : 0.2,
        "veto_transition" : True,
    },
    "muons" : {
        "pt" : 5.0,
        "eta" : 2.4,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : None,
        "pfRelIso03_all" : 0.3,
        "dr_photons" : 0.2
    },
    "taus" : {
        "pt" : 18.0,
        "eta" : 2.3,
        "dz" : 0.2,
        "deep_tau_vs_ele" : 1,
        "deep_tau_vs_mu" : 0,
        "deep_tau_vs_jet" : 7,
        "dr_photons" : 0.2,
        "dr_electrons" : 0.2,
        "dr_muons" : 0.2
    },
    "jets" : {
        "pt" : 25.0,
        "eta" : 2.4,
        "looseID" : True,
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
        "dr_taus" : 0.4,
        "dr_iso_tracks" : 0.4
    },
    "photon_mvaID" : -0.7
}

class TTHHPreselTagger(Tagger):
    """
    Preselection Tagger for the ttHH->ggXX analysis.
    """
    def __init__(self, name = "tthh_ggxx_tagger", options = {}, is_data = None, year = None):
        super(TTHHPreselTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS 
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, events):
        #################################
        ### HH->ggTauTau Preselection ###
        #################################

        ### Presel step 1 : select objects ###
        
        # Electrons
        electron_cut = lepton_selections.select_electrons(
                electrons = events.Electron,
                options = self.options["electrons"],
                clean = {
                    "photons" : {
                        "objects" : events.Diphoton.Photon,
                        "min_dr" : self.options["electrons"]["dr_photons"]
                    }
                },
                name = "SelectedElectron",
                tagger = self
        )

        electrons = awkward_utils.add_field(
                events = events,
                name = "SelectedElectron",
                data = events.Electron[electron_cut]
        )

        # Muons
        muon_cut = lepton_selections.select_muons(
                muons = events.Muon,
                options = self.options["muons"],
                clean = {
                    "photons" : {
                        "objects" : events.Diphoton.Photon,
                        "min_dr" : self.options["muons"]["dr_photons"]
                    }
                },
                name = "SelectedMuon",
                tagger = self
        )

        muons = awkward_utils.add_field(
                events = events,
                name = "SelectedMuon",
                data = events.Muon[muon_cut]
        )

        # Taus
        tau_cut = tau_selections.select_taus(
                taus = events.Tau,
                options = self.options["taus"],
                clean = {
                    "photons" : {
                        "objects" : events.Diphoton.Photon,
                        "min_dr" : self.options["taus"]["dr_photons"]
                    },
                    "electrons" : {
                        "objects" : events.SelectedElectron,
                        "min_dr" : self.options["taus"]["dr_electrons"]
                    },
                    "muons" : {
                        "objects" : events.SelectedMuon,
                        "min_dr" : self.options["taus"]["dr_muons"]
                    }
                },
                name = "AnalysisTau",
                tagger = self
        )

        taus = awkward_utils.add_field(
                events = events,
                name = "AnalysisTau",
                data = events.Tau[tau_cut]
        )

        # Jets
        jet_cut = jet_selections.select_jets(
                jets = events.Jet,
                options = self.options["jets"],
                clean = {
                    "photons" : {
                        "objects" : events.Diphoton.Photon,
                        "min_dr" : self.options["jets"]["dr_photons"]
                    },
                    "electrons" : {
                        "objects" : events.SelectedElectron,
                        "min_dr" : self.options["jets"]["dr_electrons"]
                    },
                    "muons" : {
                        "objects" : events.SelectedMuon,
                        "min_dr" : self.options["jets"]["dr_muons"]
                    },
                    "taus" : {
                        "objects" : events.AnalysisTau,
                        "min_dr" : self.options["jets"]["dr_taus"]
                    }
                },
                name = "SelectedJet",
                tagger = self
        )

        jets = awkward_utils.add_field(
                events = events,
                name = "SelectedJet",
                data = events.Jet[jet_cut]
        )

        awkward_utils.add_object_fields(
                events = events,
                name = "jet",
                objects = jets,
                n_objects = 6,
                dummy_value = DUMMY_VALUE
        )

        bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis = 1, ascending = False)]
        awkward_utils.add_object_fields(
                events = events,
                name = "b_jet",
                objects = bjets,
                n_objects = 4,
                fields = ["btagDeepFlavB"],
                dummy_value = DUMMY_VALUE
        )

        # Merge all taus, electrons, muons into one object (can be distinguished through "id" field)
        taus = awkward.with_field(taus, awkward.ones_like(taus.pt) * 15, "id")
        electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 11, "id")
        muons = awkward.with_field(muons, awkward.ones_like(muons.pt) * 13, "id") 

        # Leptons includes taus/eles/muons
        leptons = awkward.concatenate(
            [taus, electrons, muons],
            axis = 1
        )

        # Sort leptons by pT
        leptons = leptons[awkward.argsort(leptons.pt, ascending=False, axis=1)]
        leptons = awkward.Array(leptons, with_name = "Momentum4D")

        # Add leading 2 lepton fields to events array
        awkward_utils.add_object_fields(
            events = events,
            name = "lepton",
            objects = leptons,
            n_objects = 2,
            dummy_value = DUMMY_VALUE
        )

        # Add object multiplicity variables
        n_electrons = awkward.num(electrons)
        awkward_utils.add_field(events, "n_electrons", n_electrons)
        
        n_muons = awkward.num(muons)
        awkward_utils.add_field(events, "n_muons", n_muons)
        
        n_taus = awkward.num(taus)
        awkward_utils.add_field(events, "n_taus", n_taus)

        n_leptons = n_electrons + n_muons + n_taus
        awkward_utils.add_field(events, "n_leptons", n_leptons)
        
        n_jets = awkward.num(jets)
        awkward_utils.add_field(events, "n_jets", n_jets)

        # Add pho/dipho variables
        events[("Diphoton", "pt_mgg")] = events.Diphoton.pt / events.Diphoton.mass
        events[("LeadPhoton", "pt_mgg")] = events.LeadPhoton.pt / events.Diphoton.mass
        events[("SubleadPhoton", "pt_mgg")] = events.SubleadPhoton.pt / events.Diphoton.mass

        # Preselection cuts
        pho_id = (events.LeadPhoton.mvaID > self.options["photon_mvaID"]) & (events.SubleadPhoton.mvaID > self.options["photon_mvaID"])

        hadronic_cut = (n_leptons == 0) & (n_jets >= 4)
        semilep_cut  = (n_leptons == 1) & (n_jets >= 3)
        dilep_cut    = (n_leptons >= 2) & (n_jets >= 2)

        presel_cut = pho_id & (hadronic_cut | semilep_cut | dilep_cut) 
        self.register_cuts(
            names = ["photon ID MVA", "hadronic", "semileptonic", "dileptonic", "all"],
            results = [pho_id, hadronic_cut, semilep_cut, dilep_cut, presel_cut]
        )

        return presel_cut, events 
