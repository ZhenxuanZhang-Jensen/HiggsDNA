import awkward
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.selections import object_selections, lepton_selections, jet_selections, tau_selections, physics_utils
from higgs_dna.utils import awkward_utils, misc_utils

DUMMY_VALUE = -999.
DEFAULT_OPTIONS = {
    "electrons" : {
        "pt" : 10.0,
        "eta" : 2.5,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "WP90",
        "dr_photons" : 0.2,
        "veto_transition" : True,
    },
    "muons" : {
        "pt" : 15.0,
        "eta" : 2.4,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "medium",
        "pfRelIso03_all" : 0.3,
        "dr_photons" : 0.2
    },
    "taus" : {
        "pt" : 20.0,
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

        bjet_sort_idx = awkward.argsort(jets.btagDeepFlavB, axis = 1, ascending = False)
        bjets = jets[bjet_sort_idx]
        bjets["jet_idx"] = bjet_sort_idx
        awkward_utils.add_object_fields(
                events = events,
                name = "b_jet",
                objects = bjets,
                n_objects = 4,
                fields = ["btagDeepFlavB", "jet_idx"],
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

        # Add leading 3 lepton fields to events array
        awkward_utils.add_object_fields(
            events = events,
            name = "lepton",
            objects = leptons,
            n_objects = 3,
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
        events[("Diphoton", "dPhi")] = events.LeadPhoton.deltaphi(events.SubleadPhoton)
        events[("Diphoton", "helicity")] = physics_utils.abs_cos_theta_parentCM(events.LeadPhoton, events.SubleadPhoton)

        ### Ditau candidate reconstruction ###
        # 1. Create ditau candidates: all possible pairs of two leptons, with leptons = {taus, electrons, muons}
        ditau_pairs = awkward.combinations(leptons, 2, fields = ["LeadTauCand", "SubleadTauCand"])

        # 2. Keep only OS ditau pairs 
        os_cut = ditau_pairs.LeadTauCand.charge * ditau_pairs.SubleadTauCand.charge == -1
        ditau_pairs = ditau_pairs[os_cut]

        # 3. Assign priority to ditau pairs by lepton flavor
        # If there is more than one ditau candidate in an event, we first give preference by lepton flavor:
        # from highest to lowest priority : tau/tau, tau/mu, tau/ele, mu/ele, mu/mu, ele/ele 
        ditau_pairs = awkward.with_field(ditau_pairs, ditau_pairs.LeadTauCand.id * ditau_pairs.SubleadTauCand.id, "ProdID")
        ditau_pairs = awkward.with_field(ditau_pairs, awkward.ones_like(ditau_pairs.LeadTauCand.pt) * DUMMY_VALUE, "priority") 

        # tau_h / tau_h : 15 * 15 = 225
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 225, awkward.ones_like(ditau_pairs["priority"]) * 1, ditau_pairs["priority"])
        # tau_h / mu : 15 * 13 = 195
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 195, awkward.ones_like(ditau_pairs["priority"]) * 2, ditau_pairs["priority"])
        # tau_h / ele : 15 * 11 = 165
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 165, awkward.ones_like(ditau_pairs["priority"]) * 3, ditau_pairs["priority"])
        # mu / ele : 13 * 11 = 143
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 143, awkward.ones_like(ditau_pairs["priority"]) * 4, ditau_pairs["priority"])
        # mu / mu : 13 * 13 = 169
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 169, awkward.ones_like(ditau_pairs["priority"]) * 5, ditau_pairs["priority"])
        # ele / ele : 11 * 11 = 121
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 121, awkward.ones_like(ditau_pairs["priority"]) * 6, ditau_pairs["priority"])
        # tau / iso track : 15 * 1 = 15
        ditau_pairs["priority"] = awkward.where(ditau_pairs.ProdID == 15, awkward.ones_like(ditau_pairs["priority"]) * 7, ditau_pairs["priority"])

        # 4. Select only the highest priority di-tau candidate(s) in each event
        ditau_pairs = ditau_pairs[ditau_pairs.priority == awkward.min(abs(ditau_pairs.priority), axis = 1)]

        # 5. If still more than one ditau candidate in an event, take the one with m_vis closest to mH = 125 GeV
        ditau_pairs["ditau"] = ditau_pairs.LeadTauCand + ditau_pairs.SubleadTauCand
        ditau_pairs[("ditau", "dR")] = ditau_pairs.LeadTauCand.deltaR(ditau_pairs.SubleadTauCand)
        ditau_pairs[("ditau", "helicity")] = physics_utils.abs_cos_theta_parentCM(ditau_pairs.LeadTauCand, ditau_pairs.SubleadTauCand)
        # Make sure the pt/eta/phi/mass are explicitly calculated by vector so they will be saved in outputs
        ditau_pairs[("ditau", "pt")] = ditau_pairs.ditau.pt
        ditau_pairs[("ditau", "eta")] = ditau_pairs.ditau.eta
        ditau_pairs[("ditau", "phi")] = ditau_pairs.ditau.phi
        ditau_pairs[("ditau", "mass")] = ditau_pairs.ditau.mass
        awkward_utils.add_object_fields(
            events = events,
            name = "ditau",
            objects = ditau_pairs.ditau,
            n_objects = 1,
            dummy_value = DUMMY_VALUE
        )

        if awkward.any(awkward.num(ditau_pairs) >= 2):
            ditau_pairs = ditau_pairs[awkward.argsort(abs(ditau_pairs.ditau.mass - 125), axis = 1)] # if so, take the one with m_vis closest to mH
        ditau_pairs = awkward.firsts(ditau_pairs)

        # Preselection cuts
        pho_id = (events.LeadPhoton.mvaID > self.options["photon_mvaID"]) & (events.SubleadPhoton.mvaID > self.options["photon_mvaID"])

        hadronic_cut = (n_leptons == 0) & (n_jets >= 4)
        semilep_cut  = (n_leptons == 1) & (n_jets >= 2)
        dilep_cut    = (n_leptons >= 2) & (n_jets >= 0)

        presel_cut = pho_id & (hadronic_cut | semilep_cut | dilep_cut) 
        self.register_cuts(
            names = ["photon ID MVA", "hadronic", "semileptonic", "dileptonic", "all"],
            results = [pho_id, hadronic_cut, semilep_cut, dilep_cut, presel_cut]
        )

        return presel_cut, events 
