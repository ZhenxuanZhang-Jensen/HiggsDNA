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
        "pt" : 10.0,
        "dxy" : 0.2,
        "dz" : 0.1,
        "ch_iso" : 5.0,
        "ch_rel_iso" : 0.2,
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
    }
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

        # Preselection
        n_electrons = awkward.num(electrons)
        awkward_utils.add_field(syst_events, "n_electrons", n_electrons)
        
        n_muons = awkward.num(muons)
        awkward_utils.add_field(syst_events, "n_muons", n_muons)
        
        n_leptons = n_electrons + n_muons
        awkward_utils.add_field(syst_events, "n_leptons", n_leptons)
        
        n_taus = awkward.num(taus)
        awkward_utils.add_field(syst_events, "n_taus", n_taus)
        
        n_jets = awkward.num(jets)
        awkward_utils.add_field(syst_events, "n_jets", n_jets)

        tau1_lep0 = (n_taus == 1) & (n_leptons == 0)
        tau2_lep0 = (n_taus == 2) & (n_leptons == 0)
        tau1_lep1 = (n_taus == 1) & (n_leptons == 1)
        tau0_lep2 = (n_taus == 0) & (n_leptons == 2)

        presel_cut = tau1_lep0 | tau2_lep0 | tau1_lep1 | tau0_lep2
        self.register_cuts(
            names = ["1tau_0lep", "2tau_0lep", "1tau_1lep", "0tau_2lep", "inclusive"],
            results = [tau1_lep0, tau2_lep0, tau1_lep1, tau0_lep2, presel_cut]
        )

        return presel_cut, syst_events 
