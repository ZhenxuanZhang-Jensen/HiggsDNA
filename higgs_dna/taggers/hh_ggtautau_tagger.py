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
        "veto_transition" : True
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
    "iso_tracks" : {
        "pt" : 5.0,
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
    "photon_mvaID" : -0.7
}

class HHggTauTauTagger(Tagger):
    """
    Tagger for the non-resonant HH->ggTauTau analysis.
    """
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(HHggTauTauTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS 
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, syst_tag, syst_events):
        #################################
        ### HH->ggTauTau Preselection ###
        #################################

        ### Presel step 1 : select objects ###
        
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

        # Taus
        tau_cut = tau_selections.select_taus(
                taus = syst_events.Tau,
                options = self.options["taus"],
                clean = {
                    "photons" : {
                        "objects" : syst_events.Diphoton.Photon,
                        "min_dr" : self.options["taus"]["dr_photons"]
                    },
                    "electrons" : {
                        "objects" : syst_events.SelectedElectron,
                        "min_dr" : self.options["taus"]["dr_electrons"]
                    },
                    "muons" : {
                        "objects" : syst_events.SelectedMuon,
                        "min_dr" : self.options["taus"]["dr_muons"]
                    }
                },
                name = "AnalysisTau",
                tagger = self
        )

        taus = awkward_utils.add_field(
                events = syst_events,
                name = "AnalysisTau",
                data = syst_events.Tau[tau_cut]
        )

        # IsoTrack
        syst_events = awkward.with_field(
                base = syst_events,
                what = awkward.zeros_like(syst_events.IsoTrack.pt),
                where = ("IsoTrack", "mass")
        )
        syst_events = awkward.with_field(
                base = syst_events,
                what = awkward.nan_to_num(syst_events.IsoTrack.pdgId / abs(syst_events.IsoTrack.pdgId)),
                where = ("IsoTrack", "charge")
        )

        
        iso_track_cut = tau_selections.select_iso_tracks(
                iso_tracks = syst_events.IsoTrack,
                options = self.options["iso_tracks"],
                clean = {
                    "photons" : {
                        "objects" : syst_events.Diphoton.Photon,
                        "min_dr" : self.options["iso_tracks"]["dr_photons"]
                    },
                    "electrons" : {
                        "objects" : syst_events.SelectedElectron,
                        "min_dr" : self.options["iso_tracks"]["dr_electrons"]
                    },
                    "muons" : {
                        "objects" : syst_events.SelectedMuon,
                        "min_dr" : self.options["iso_tracks"]["dr_muons"]
                    },
                    "taus" : {
                        "objects" : syst_events.AnalysisTau,
                        "min_dr" : self.options["iso_tracks"]["dr_taus"]
                    }
                },
                name = "SelectedIsoTrack",
                tagger = self
        )

        iso_tracks = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedIsoTrack",
                data = syst_events.IsoTrack[iso_track_cut]
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
                    },
                    "taus" : {
                        "objects" : syst_events.AnalysisTau,
                        "min_dr" : self.options["jets"]["dr_taus"]
                    },
                    "iso_tracks" : {
                        "objects" : syst_events.SelectedIsoTrack,
                        "min_dr" : self.options["jets"]["dr_iso_tracks"]
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

        # Add object fields to events array
        for objects, name in zip([electrons, muons, taus, jets], ["electron", "muon", "tau", "jet"]):
            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = name,
                    objects = objects,
                    n_objects = 2,
                    dummy_value = DUMMY_VALUE
            )

        awkward_utils.add_object_fields(
                events = syst_events,
                name = "iso_track",
                objects = iso_tracks,
                n_objects = 1,
                dummy_value = DUMMY_VALUE,
        )

        n_electrons = awkward.num(electrons)
        awkward_utils.add_field(syst_events, "n_electrons", n_electrons)
        
        n_muons = awkward.num(muons)
        awkward_utils.add_field(syst_events, "n_muons", n_muons)
        
        n_leptons = n_electrons + n_muons
        awkward_utils.add_field(syst_events, "n_leptons", n_leptons)
        
        n_taus = awkward.num(taus)
        awkward_utils.add_field(syst_events, "n_taus", n_taus)
        
        n_iso_tracks = awkward.num(iso_tracks)
        awkward_utils.add_field(syst_events, "n_iso_tracks", n_iso_tracks)

        n_jets = awkward.num(jets)
        awkward_utils.add_field(syst_events, "n_jets", n_jets)


        ### Presel step 2: Z veto ###
        ee_pairs = awkward.combinations(electrons, 2, fields = ["LeadLepton", "SubleadLepton"])
        mumu_pairs = awkward.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
        dilep_pairs = awkward.concatenate(
                [ee_pairs, mumu_pairs],
                axis = 1
        )

        lead_lep_p4 = vector.awk({
            "pt" : dilep_pairs.LeadLepton.pt,
            "eta" : dilep_pairs.LeadLepton.eta,
            "phi" : dilep_pairs.LeadLepton.phi,
            "mass" : dilep_pairs.LeadLepton.mass
        })
        sublead_lep_p4 = vector.awk({
            "pt" : dilep_pairs.SubleadLepton.pt,
            "eta" : dilep_pairs.SubleadLepton.eta,
            "phi" : dilep_pairs.SubleadLepton.phi,
            "mass" : dilep_pairs.SubleadLepton.mass
        })
        z_candidates = lead_lep_p4 + sublead_lep_p4 

        os_cut = dilep_pairs["LeadLepton"].charge * dilep_pairs["SubleadLepton"].charge == -1
        z_mass_cut = (z_candidates.mass > self.options["z_veto"][0]) & (z_candidates.mass < self.options["z_veto"][1])

        z_veto = ~(os_cut & z_mass_cut) # z-veto on individual z candidates (in principle, can be more than 1 per event)
        z_veto = (awkward.num(z_candidates) == 0) | (awkward.any(z_veto, axis = 1)) # if any z candidate in an event fails the veto, the event is vetoed. If the event does not have any z candidates, we do not veto
        

        ### Presel step 3: construct di-tau candidates and assign to category ###
        taus = awkward.with_field(taus, awkward.ones_like(taus.pt) * 15, "id")
        electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 11, "id")
        muons = awkward.with_field(muons, awkward.ones_like(muons.pt) * 13, "id")
        iso_tracks = awkward.with_field(iso_tracks, awkward.ones_like(iso_tracks.pt) * 1, "id") 
        
        
        tau_candidates = awkward.concatenate(
            [taus, electrons, muons, iso_tracks],
            axis = 1
        )
        tau_candidates = tau_candidates[awkward.argsort(tau_candidates.pt, ascending = False, axis = 1)] 
        
        # Create ditau candidates: all possible pairs of two objects, with objects = {taus, electrons, muons, iso_tracks}
        tau_candidate_pairs = awkward.combinations(tau_candidates, 2, fields = ["LeadTauCand", "SubleadTauCand"])

        # Trim same-sign ditau candidates
        os_cut = tau_candidate_pairs.LeadTauCand.charge * tau_candidate_pairs.SubleadTauCand.charge == -1
        tau_candidate_pairs = tau_candidate_pairs[os_cut]

        # If there is more than one ditau candidate in an event, we first give preference by lepton flavor:
        # from highest to lowest priority : tau/tau, tau/mu, tau/ele, mu/ele, mu/mu, ele/ele, tau/iso_track 
        tau_candidate_pairs = awkward.with_field(tau_candidate_pairs, tau_candidate_pairs["LeadTauCand"].id * tau_candidate_pairs["SubleadTauCand"].id, "ProdID")
        tau_candidate_pairs = awkward.with_field(tau_candidate_pairs, awkward.ones_like(tau_candidate_pairs.LeadTauCand.pt) * -999, "priority")

        # NOTE from sam: when you add a field like "priority" and you want to modify its value, you MUST do events["priority"] = blah and NOT events.priority = blah. The latter is very bad and will not work, I am not sure why.
        # tau_h / tau_h : 15 * 15 = 225
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 225, awkward.ones_like(tau_candidate_pairs["priority"]) * 1, tau_candidate_pairs["priority"])
        # tau_h / mu : 15 * 13 = 195
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 195, awkward.ones_like(tau_candidate_pairs["priority"]) * 2, tau_candidate_pairs["priority"])
        # tau_h / ele : 15 * 11 = 165
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 165, awkward.ones_like(tau_candidate_pairs["priority"]) * 3, tau_candidate_pairs["priority"])
        # mu / ele : 13 * 11 = 143
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 143, awkward.ones_like(tau_candidate_pairs["priority"]) * 4, tau_candidate_pairs["priority"])
        # mu / mu : 13 * 13 = 169
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 169, awkward.ones_like(tau_candidate_pairs["priority"]) * 5, tau_candidate_pairs["priority"])
        # ele / ele : 11 * 11 = 121
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 121, awkward.ones_like(tau_candidate_pairs["priority"]) * 6, tau_candidate_pairs["priority"])
        # tau / iso track : 15 * 1 = 15
        tau_candidate_pairs["priority"] = awkward.where(tau_candidate_pairs.ProdID == 15, awkward.ones_like(tau_candidate_pairs["priority"]) * 7, tau_candidate_pairs["priority"])

        # Select only the highest priority di-tau candidate(s) in each event
        tau_candidate_pairs = tau_candidate_pairs[tau_candidate_pairs.priority == awkward.min(abs(tau_candidate_pairs.priority), axis = 1)]

        # If still more than one ditau candidate in an event, take the one with m_vis closest to mH = 125 GeV
        tau_candidates_lead_lepton_p4 = vector.awk({
            "pt" : tau_candidate_pairs.LeadTauCand.pt,
            "eta" : tau_candidate_pairs.LeadTauCand.eta,
            "phi" : tau_candidate_pairs.LeadTauCand.phi,
            "mass" : tau_candidate_pairs.LeadTauCand.mass
        })
        tau_candidates_sublead_lepton_p4 = vector.awk({
            "pt" : tau_candidate_pairs.SubleadTauCand.pt,
            "eta" : tau_candidate_pairs.SubleadTauCand.eta,
            "phi" : tau_candidate_pairs.SubleadTauCand.phi,
            "mass" : tau_candidate_pairs.SubleadTauCand.mass
        })

        ditau_candidates = tau_candidates_lead_lepton_p4 + tau_candidates_sublead_lepton_p4
        ditau_candidates["dR"] = tau_candidates_lead_lepton_p4.deltaR(tau_candidates_sublead_lepton_p4)

        tau_candidate_pairs["ditau"] = ditau_candidates

        tau_candidate_pairs = tau_candidate_pairs[awkward.argsort(abs(tau_candidate_pairs.ditau.mass - 125), axis = 1)]
        tau_candidate_pairs = awkward.firsts(tau_candidate_pairs)

        # Add ditau-related fields to array
        for field in ["pt", "eta", "phi", "mass", "charge", "id"]:
            if not field in ["charge", "id"]:
                awkward_utils.add_field(
                        syst_events,
                        "ditau_%s" % field,
                        awkward.fill_none(getattr(tau_candidate_pairs.ditau, field), DUMMY_VALUE)
                )
            awkward_utils.add_field(
                    syst_events,
                    "ditau_lead_lepton_%s" % field,
                    awkward.fill_none(tau_candidate_pairs.LeadTauCand[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                    syst_events,
                    "ditau_sublead_lepton_%s" % field,
                    awkward.fill_none(tau_candidate_pairs.SubleadTauCand[field], DUMMY_VALUE)
            )
        awkward_utils.add_field(
                syst_events,
                "ditau_dR",
                awkward.fill_none(tau_candidate_pairs.ditau.dR, DUMMY_VALUE)
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
            category = awkward.where(tau_candidate_pairs.ProdID == prod_id, awkward.ones_like(category) * category_number, category)

        # Now assign 1tau / 0 lep, 0 isotrack events to category 8
        category = awkward.fill_none(category, 0)
        category = awkward.where((category == 0) & (n_taus >= 1), awkward.ones_like(category) * 8, category)
        awkward_utils.add_field(syst_events, "category", category) 

        # Events must fall into one of 8 categories
        category_cut = category > 0

        # Photon ID cut
        pho_id = (syst_events.LeadPhoton.mvaID > self.options["photon_mvaID"]) & (syst_events.SubleadPhoton.mvaID > self.options["photon_mvaID"])

        presel_cut = category_cut & z_veto & pho_id 
        self.register_cuts(
            names = ["category", "z_veto", "photon ID MVA", "all cuts"],
            results = [category_cut, z_veto, pho_id, presel_cut]
        )

        return presel_cut, syst_events 
