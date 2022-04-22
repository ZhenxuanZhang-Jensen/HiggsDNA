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
    "iso_tracks" : {
        "pt" : 5.0,
        "eta" : 5.0,
        "dxy" : 0.2,
        "dz" : 0.1,
        "dr_photons" : 0.2,
        "dr_electrons" : 0.2,
        "dr_muons" : 0.2,
        "dr_taus" : 0.2
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
    "z_veto" : [80., 100.],
    "m_llg_veto_window" : 25,
    "photon_mvaID" : -0.7
}

class HHggTauTauPreselTagger(Tagger):
    """
    Preselection Tagger for the HH->ggTauTau analysis.
    """
    def __init__(self, name = "hh_ggtautau_tagger", options = {}, is_data = None, year = None):
        super(HHggTauTauPreselTagger, self).__init__(name, options, is_data, year)

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

        # IsoTrack
        # Add missing iso track fields
        events = awkward.with_field(
                base = events,
                what = awkward.zeros_like(events.IsoTrack.pt),
                where = ("IsoTrack", "mass")
        )
        events = awkward.with_field(
                base = events,
                what = awkward.nan_to_num(events.IsoTrack.pdgId / abs(events.IsoTrack.pdgId)),
                where = ("IsoTrack", "charge")
        )

        
        iso_track_cut = tau_selections.select_iso_tracks(
                iso_tracks = events.IsoTrack,
                options = self.options["iso_tracks"],
                clean = {
                    "photons" : {
                        "objects" : events.Diphoton.Photon,
                        "min_dr" : self.options["iso_tracks"]["dr_photons"]
                    },
                    "electrons" : {
                        "objects" : events.SelectedElectron,
                        "min_dr" : self.options["iso_tracks"]["dr_electrons"]
                    },
                    "muons" : {
                        "objects" : events.SelectedMuon,
                        "min_dr" : self.options["iso_tracks"]["dr_muons"]
                    },
                    "taus" : {
                        "objects" : events.AnalysisTau,
                        "min_dr" : self.options["iso_tracks"]["dr_taus"]
                    }
                },
                name = "SelectedIsoTrack",
                tagger = self
        )

        iso_tracks = awkward_utils.add_field(
                events = events,
                name = "SelectedIsoTrack",
                data = events.IsoTrack[iso_track_cut]
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
                    },
                    "iso_tracks" : {
                        "objects" : events.SelectedIsoTrack,
                        "min_dr" : self.options["jets"]["dr_iso_tracks"]
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
        awkward_utils.add_object_fields(
                events = events,
                name = "b_jet",
                objects = bjets,
                n_objects = 2,
                fields = ["btagDeepFlavB"],
                dummy_value = DUMMY_VALUE
        )

        # Add object fields to events array
        for objects, name in zip([electrons, muons, taus, jets], ["electron", "muon", "tau", "jet"]):
            awkward_utils.add_object_fields(
                    events = events,
                    name = name,
                    objects = objects,
                    n_objects = 2,
                    dummy_value = DUMMY_VALUE
            )

        awkward_utils.add_object_fields(
                events = events,
                name = "iso_track",
                objects = iso_tracks,
                n_objects = 1,
                dummy_value = DUMMY_VALUE,
        )

        n_electrons = awkward.num(electrons)
        awkward_utils.add_field(events, "n_electrons", n_electrons)
        
        n_muons = awkward.num(muons)
        awkward_utils.add_field(events, "n_muons", n_muons)
        
        n_leptons = n_electrons + n_muons
        awkward_utils.add_field(events, "n_leptons", n_leptons)
        
        n_taus = awkward.num(taus)
        awkward_utils.add_field(events, "n_taus", n_taus)
        
        n_iso_tracks = awkward.num(iso_tracks)
        awkward_utils.add_field(events, "n_iso_tracks", n_iso_tracks)

        n_jets = awkward.num(jets)
        awkward_utils.add_field(events, "n_jets", n_jets)


        ### Presel step 2: Z veto ###
        # 2.1 Register objects as 4 vectors for appropriate overloading of operators
        electrons = awkward.Array(electrons, with_name = "Momentum4D")
        muons = awkward.Array(muons, with_name = "Momentum4D")

        # 2.2 Create all SF dilep pairs
        ee_pairs = awkward.combinations(electrons, 2, fields = ["LeadLepton", "SubleadLepton"])
        mumu_pairs = awkward.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
        dilep_pairs = awkward.concatenate(
                [ee_pairs, mumu_pairs],
                axis = 1
        )

        # 2.3 Make z candidates
        dilep_pairs["z_candidate"] = dilep_pairs.LeadLepton + dilep_pairs.SubleadLepton

        # 2.4 Keep only OS cands in the specified mass range
        os_cut = dilep_pairs.LeadLepton.charge * dilep_pairs.SubleadLepton.charge == -1
        z_mass_cut = (dilep_pairs.z_candidate.mass >= self.options["z_veto"][0]) & (dilep_pairs.z_candidate.mass <= self.options["z_veto"][1])

        z_cut = os_cut & z_mass_cut
        dilep_pairs = dilep_pairs[z_cut]

        # 2.5 Event level cut
        has_z_cand = awkward.num(dilep_pairs) >= 1 # veto any event that has a Z candidate
        ee_event = awkward.num(electrons) >= 2
        mm_event = awkward.num(muons) >= 2
        z_veto = ~(has_z_cand & (ee_event | mm_event))

        ### Presel step 3: construct di-tau candidates and assign to category ###
        # 3.1 Merge all taus, electrons, muons, and iso tracks
        taus = awkward.with_field(taus, awkward.ones_like(taus.pt) * 15, "id")
        electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 11, "id")
        muons = awkward.with_field(muons, awkward.ones_like(muons.pt) * 13, "id")
        iso_tracks = awkward.with_field(iso_tracks, awkward.ones_like(iso_tracks.pt) * 1, "id") 
        
        
        tau_candidates = awkward.concatenate(
            [taus, electrons, muons, iso_tracks],
            axis = 1
        )
        tau_candidates = tau_candidates[awkward.argsort(tau_candidates.pt, ascending = False, axis = 1)] 
        tau_candidates = awkward.Array(tau_candidates, with_name = "Momentum4D")

        awkward_utils.add_object_fields(
                events = events,
                name = "tau_candidate",
                objects = tau_candidates,
                n_objects = 3,
                dummy_value = DUMMY_VALUE
        ) 

        # 3.2 Create ditau candidates: all possible pairs of two objects, with objects = {taus, electrons, muons, iso_tracks}
        ditau_pairs = awkward.combinations(tau_candidates, 2, fields = ["LeadTauCand", "SubleadTauCand"])

        # 3.3 Keep only OS ditau pairs 
        os_cut = ditau_pairs.LeadTauCand.charge * ditau_pairs.SubleadTauCand.charge == -1
        ditau_pairs = ditau_pairs[os_cut]

        # 3.4 Assign priority to ditau pairs by lepton flavor
        # If there is more than one ditau candidate in an event, we first give preference by lepton flavor:
        # from highest to lowest priority : tau/tau, tau/mu, tau/ele, mu/ele, mu/mu, ele/ele, tau/iso_track 
        ditau_pairs = awkward.with_field(ditau_pairs, ditau_pairs.LeadTauCand.id * ditau_pairs.SubleadTauCand.id, "ProdID")
        ditau_pairs = awkward.with_field(ditau_pairs, awkward.ones_like(ditau_pairs.LeadTauCand.pt) * -999, "priority")

        # NOTE from sam: when you add a field like "priority" and you want to modify its value, you MUST do events["priority"] = blah and NOT events.priority = blah. The latter is very bad and will not work, I am not sure why.
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

        # 3.5 Select only the highest priority di-tau candidate(s) in each event
        ditau_pairs = ditau_pairs[ditau_pairs.priority == awkward.min(abs(ditau_pairs.priority), axis = 1)]

        # 3.6 If still more than one ditau candidate in an event, take the one with m_vis closest to mH = 125 GeV
        ditau_pairs["ditau"] = ditau_pairs.LeadTauCand + ditau_pairs.SubleadTauCand
        ditau_pairs["ditau", "dR"] = ditau_pairs.LeadTauCand.deltaR(ditau_pairs.SubleadTauCand)

        if awkward.any(awkward.num(ditau_pairs) >= 2):
            ditau_pairs = ditau_pairs[awkward.argsort(abs(ditau_pairs.ditau.mass - 125), axis = 1)] # if so, take the one with m_vis closest to mH
        ditau_pairs = awkward.firsts(ditau_pairs)

        # Add ditau-related fields to array
        for field in ["pt", "eta", "phi", "mass", "charge", "id"]:
            if not field in ["charge", "id"]:
                awkward_utils.add_field(
                        events,
                        "ditau_%s" % field,
                        awkward.fill_none(getattr(ditau_pairs.ditau, field), DUMMY_VALUE)
                )
            awkward_utils.add_field(
                    events,
                    "ditau_lead_lepton_%s" % field,
                    awkward.fill_none(ditau_pairs.LeadTauCand[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                    events,
                    "ditau_sublead_lepton_%s" % field,
                    awkward.fill_none(ditau_pairs.SubleadTauCand[field], DUMMY_VALUE)
            )
        awkward_utils.add_field(
                events,
                "ditau_dR",
                awkward.fill_none(ditau_pairs.ditau.dR, DUMMY_VALUE)
        )

        # Now assign the selected tau candidate pair in each event to a category integer
        category_map = {
            225 : 3, # tau/tau
            195 : 1, # tau/mu
            165 : 2, # tau/ele
            143 : 6, # mu/ele
            169 : 4, # mu/mu
            121 : 5, # ele/ele
            15  : 7  # tau/iso_track
        }

        category = awkward.zeros_like(n_jets)

        for prod_id, category_number in category_map.items():
            category = awkward.where(ditau_pairs.ProdID == prod_id, awkward.ones_like(category) * category_number, category)

        # Now assign 1tau / 0 lep, 0 isotrack events to category 8
        category = awkward.fill_none(category, 0)
        category = awkward.where((category == 0) & (n_taus >= 1), awkward.ones_like(category) * 8, category)
        awkward_utils.add_field(events, "category", category) 

        # Events must fall into one of 8 categories
        category_cut = category > 0

        # Photon ID cut
        pho_id = (events.LeadPhoton.mvaID > self.options["photon_mvaID"]) & (events.SubleadPhoton.mvaID > self.options["photon_mvaID"])

        # Veto on m_llgamma to reject Z->eeg and Z->mmg events 
        dilep_lead_photon = ditau_pairs.ditau + events.LeadPhoton
        dilep_sublead_photon = ditau_pairs.ditau + events.SubleadPhoton

        awkward_utils.add_field(events, "dilep_leadpho_mass", dilep_lead_photon.mass) 
        awkward_utils.add_field(events, "dilep_subleadpho_mass", dilep_sublead_photon.mass) 

        m_llg_veto_lead = abs(dilep_lead_photon.mass - 91.18) < self.options["m_llg_veto_window"]
        m_llg_veto_sublead = abs(dilep_sublead_photon.mass - 91.18) < self.options["m_llg_veto_window"]
        # The m_llg_veto_* cuts will have `None` values for events which do not have a ditau pair. Here we replace `None` values with False (if there is not a ditau pair, there is no llg system, and it cannot have a mass in the veto range)
        m_llg_veto_lead = awkward.fill_none(m_llg_veto_lead, value = False)
        m_llg_veto_sublead = awkward.fill_none(m_llg_veto_sublead, value = False)

        # Veto event if there are at least 2 OSSF leptons and they have m_llg (for either lead or sublead photon) in the z mass window
        m_llg_veto = ~(((n_muons >= 2) | (n_electrons >= 2)) & (m_llg_veto_lead | m_llg_veto_sublead)) 

        presel_cut = category_cut & z_veto & pho_id & m_llg_veto
        self.register_cuts(
            names = ["category", "z_veto", "photon ID MVA", "m_llg veto", "all cuts"],
            results = [category_cut, z_veto, pho_id, m_llg_veto, presel_cut]
        )

        return presel_cut, events 
