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
    },
    "photon_id" : 0.0,
    "btag_wp" : {
        "2016" : 0.3093,
        "2017" : 0.3040,
        "2018" : 0.2783
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

    def calculate_selection(self, events):
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

        bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis = 1, ascending = False)]
        bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]] 

        # Z-veto
        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        electrons = awkward.Array(electrons, with_name = "Momentum4D")
        muons = awkward.Array(muons, with_name = "Momentum4D")

        # Construct di-electron/di-muon pairs
        ee_pairs = awkward.combinations(
                electrons, # objects to make combinations out of
                2, # how many objects go in a combination
                fields = ["LeadLepton", "SubleadLepton"] # can access these as e.g. ee_pairs.LeadLepton.pt
        )
        mm_pairs = awkward.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
    
        # Concatenate these together
        z_cands = awkward.concatenate([ee_pairs, mm_pairs], axis = 1)
        z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton # these add as 4-vectors since we registered them as "Momentum4D" objects

        # Make Z candidate-level cuts
        os_cut = z_cands.LeadLepton.charge * z_cands.SubleadLepton.charge == -1
        mass_cut = (z_cands.ZCand.mass > 86.) & (z_cands.ZCand.mass < 96.)
        z_cut = os_cut & mass_cut
        z_cands = z_cands[z_cut] # OSSF lepton pairs with m_ll [86., 96.]

        # Make event level cut
        has_z_cand = awkward.num(z_cands) >= 1 # veto any event that has a Z candidate
        ee_event = awkward.num(electrons) >= 2
        mm_event = awkward.num(muons) >= 2
        z_veto = ~(has_z_cand & (ee_event | mm_event))

        # Preselection
        n_electrons = awkward.num(electrons)
        n_muons = awkward.num(muons)
        n_leptons = n_electrons + n_muons

        n_jets = awkward.num(jets)
        n_bjets = awkward.num(bjets)

        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (events.SubleadPhoton.mvaID > self.options["photon_id"]) 

        # Hadronic presel
        hadronic = (n_leptons == 0) & (n_jets >= 4) & (n_bjets >= 1) 

        # Leptonic presel
        leptonic = (n_leptons >= 1) & (n_jets >= 2)

        presel_cut = (hadronic | leptonic) & z_veto & photon_id_cut
        self.register_cuts(
            names = ["hadronic presel", "leptonic presel", "z veto", "photon ID cut", "all"],
            results = [hadronic, leptonic, z_veto, photon_id_cut, presel_cut]
        )

        return presel_cut, events
    
