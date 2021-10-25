import awkward
import vector

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
        "rel_iso" : 0.3,
        "dr_photons" : 0.2
    },
    "muons" : {
        "pt" : 5.0,
        "eta" : 2.5,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : None,
        "rel_iso" : 0.3,
        "dr_photons" : 0.2
    },
    "taus" : {
        "pt" : 18.0,
        "eta" : 2.3,
        "dz" : 0.2,
        "deeptau_vs_ele" : 1,
        "deeptau_vs_mu" : 2,
        "deeptau_vs_jet" : 4,
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
    def __init__(self, name, options = {}, sample = None):
        super(HHggTauTauTagger, self).__init__(name, options, sample)

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
                #fields = ["pt", "eta", "phi", "charge"]
        )

        # Preselection
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

        # Categories
        cat_1 = n_taus >= 2
        cat_2 = (n_taus >= 1) & (n_muons >= 1)
        cat_3 = (n_taus >= 1) & (n_electrons >= 1)
        cat_4 = (n_muons >= 1) & (n_electrons >= 1)
        cat_5 = (n_muons >= 2)
        cat_6 = (n_electrons >= 2)
        cat_7 = (n_taus >= 1) & (n_iso_tracks == 1)
        cat_8 = (n_taus >= 1)

        # Tau candidates
        syst_events = awkward.with_field(
                base = syst_events,
                what = awkward.zeros_like(syst_events.IsoTrack.pt),
                where = ("IsoTrack", "mass")
        )

        taus = awkward.with_field(taus, awkward.ones_like(taus.pt) * 15, "id")
        electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 11, "id")
        muons = awkward.with_field(muons, awkward.ones_like(muons.pt) * 13, "id")
        iso_tracks = awkward.with_field(iso_tracks, awkward.ones_like(iso_tracks.pt) * 1, "id")
 
        tau_cands = awkward.concatenate(
            [taus, electrons, muons, iso_tracks],
            axis = 1 # keep n_events constant, concatenate objects
        )

        tau_cands = tau_cands[awkward.argsort(tau_cands.pt, ascending=False, axis=1)]
        awkward_utils.add_object_fields(
            events = syst_events,
            name = "tau_candidate",
            objects = tau_cands,
            n_objects = 3,
            dummy_value = DUMMY_VALUE,
            fields = ["pt", "eta", "phi", "mass", "charge", "id"]
        )


        # Create dilep candidates
        lep1 = vector.arr({
            "pt" : syst_events["tau_candidate_1_pt"],
            "eta" : syst_events["tau_candidate_1_eta"],
            "phi" : syst_events["tau_candidate_1_phi"],
            "mass" : syst_events["tau_candidate_1_mass"]
        })
        lep2 = vector.arr({
            "pt" : syst_events["tau_candidate_2_pt"],
            "eta" : syst_events["tau_candidate_2_eta"],
            "phi" : syst_events["tau_candidate_2_phi"],
            "mass" : syst_events["tau_candidate_2_mass"]
        })
 
        dilep = lep1 + lep2
        
        m_vis = dilep.mass
        eta_vis = dilep.eta
        phi_vis = dilep.phi
        pt_vis = dilep.pt

        for field, array in zip(["mass", "eta", "phi", "pt"], [m_vis, eta_vis, phi_vis, pt_vis]):
            awkward_utils.add_field(syst_events, "ditau_%s_vis" % field, array)

        os_lep = syst_events["tau_candidate_1_charge"] * syst_events["tau_candidate_2_charge"] == -1 
        sf_lep = syst_events["tau_candidate_1_id"] == syst_events["tau_candidate_2_id"]

        z_veto = ~(os_lep & sf_lep & (m_vis > self.options["z_veto"][0]) & (m_vis < self.options["z_veto"][1]))

        # Photon ID cut
        pho_id = (syst_events.LeadPhoton.mvaID > self.options["photon_mvaID"]) & (syst_events.SubleadPhoton.mvaID > self.options["photon_mvaID"])

        for idx, cat in enumerate([cat_1,cat_2,cat_3,cat_4,cat_5,cat_6,cat_7,cat_8]):
            awkward_utils.add_field(
                syst_events,
                "cat_%d" % (idx+1),
                cat
            ) 

        presel_cut = (cat_1 | cat_2 | cat_3 | cat_4 | cat_5 | cat_6 | cat_7 | cat_8) & z_veto & pho_id
        self.register_cuts(
            names = ["2tau_0lep", "1tau_1mu", "1tau_1ele", "1mu_1ele", "0tau_2mu", "0tau_2ele", "1tau_1isotrack", "1tau", "z_veto", "photon ID MVA", "inclusive"],
            results = [cat_1, cat_2, cat_3, cat_4, cat_5, cat_6, cat_7, cat_8, z_veto, pho_id, presel_cut]
        )

        return presel_cut, syst_events 
