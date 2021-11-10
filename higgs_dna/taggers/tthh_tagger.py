import awkward
import vector
import numpy

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.selections import object_selections, lepton_selections, jet_selections, tau_selections, fatjet_selections, gen_selections
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
        "id" : "medium",
        "pfRelIso03_all" : 0.3,
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
        "dr_photons" : 0.2,
        "dr_electrons" : 0.2,
        "dr_muons" : 0.2,
        "dr_taus" : 0.2
    },
    "fatjets" : {
        "pt" : 150.,
        "eta" : 2.4,
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
        "dr_taus" : 0.4,
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
    "z_veto" : [80., 100.]
}

class TTHHTagger(Tagger):
    """
    ttHH->ggXX tagger
    """

    def __init__(self, name, options = {}, is_data = None, year = None):
        super(TTHHTagger, self).__init__(name, options, is_data, year)

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

        # Fat jets
        fatjet_cut = fatjet_selections.select_fatjets(
                fatjets = syst_events.FatJet,
                options = self.options["fatjets"],
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
                name = "SelectedFatJet",
                tagger = self
        )

        fatjets = awkward_utils.add_field(
                events = syst_events,
                name = "SelectedFatJet",
                data = syst_events.FatJet[fatjet_cut]
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

        bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis = 1, ascending = False)]

        ### Z veto ###
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

        leptons = awkward.concatenate(
            [taus, electrons, muons],
            axis = 1
        )
        leptons = leptons[awkward.argsort(leptons.pt, ascending = False, axis = 1)]

        # Create ditau candidates: all possible pairs of two objects, with objects = {taus, electrons, muons, iso_tracks}
        tau_candidate_pairs = awkward.combinations(leptons, 2, fields = ["LeadTauCand", "SubleadTauCand"])

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

        if awkward.any(awkward.num(tau_candidate_pairs) >= 2): # are there any events still with more than one ditau candidate?
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

        for objects, name in zip([leptons, jets, fatjets], ["lepton", "jet", "fatjet"]):
            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = name,
                    objects = objects,
                    n_objects = 8 if name == "jet" else 4,
                    dummy_value = DUMMY_VALUE
            )

        awkward_utils.add_object_fields(
                events = syst_events,
                name = "bjet",
                objects = bjets,
                n_objects = 8,
                fields = ["btagDeepFlavB"],
                dummy_value = DUMMY_VALUE
        )

        

        boosted_ww_to_4q_cands = fatjets[awkward.argsort(fatjets.deepTagMD_H4qvsQCD, ascending = False, axis = 1)]
        boosted_h_to_bb_cands = fatjets[awkward.argsort(fatjets.deepTagMD_HbbvsQCD, ascending = False, axis = 1)]

        for objects, name in zip([boosted_ww_to_4q_cands, boosted_h_to_bb_cands], ["h_ww_4q_boosted_cand", "h_bb_boosted_cand"]):
            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = name,
                    objects = objects,
                    n_objects = 1,
                    dummy_value = DUMMY_VALUE
            )
        
        # Gen info
        if not self.is_data:
            h_to_gg = gen_selections.select_x_to_yz(syst_events.GenPart, 25, 22, 22)   
            h_to_bb = gen_selections.select_x_to_yz(syst_events.GenPart, 25, 5, 5)
            h_to_ww = gen_selections.select_x_to_yz(syst_events.GenPart, 25, 24, 24)
            h_to_tautau = gen_selections.select_x_to_yz(syst_events.GenPart, 25, 15, 15)

            h_to_xx = awkward.concatenate([h_to_gg, h_to_bb, h_to_ww, h_to_tautau], axis = 1)
            if awkward.any(awkward.num(h_to_xx) >= 2):
                h_to_xx = h_to_xx[awkward.argsort(h_to_xx.GenParent.pt, ascending = False, axis = 1)]

            t_to_bw = gen_selections.select_x_to_yz(syst_events.GenPart, 6, 5, 24)

            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = "gen_higgs",
                    objects = h_to_xx.GenParent,
                    n_objects = 2,
                    dummy_value = DUMMY_VALUE
            )

            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = "gen_higgs_lead_child",
                    objects = h_to_xx.LeadGenChild,
                    n_objects = 2,
                    dummy_value = DUMMY_VALUE
            )

            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = "gen_higgs_sublead_child",
                    objects = h_to_xx.SubleadGenChild,
                    n_objects = 2,
                    dummy_value = DUMMY_VALUE
            )


            awkward_utils.add_object_fields(
                    events = syst_events,
                    name = "gen_top",
                    objects = t_to_bw.GenParent,
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

        n_lep_tau = n_leptons + n_taus
        awkward_utils.add_field(syst_events, "n_lep_tau", n_lep_tau)

        n_jets = awkward.num(jets)
        awkward_utils.add_field(syst_events, "n_jets", n_jets)

        n_fatjets = awkward.num(fatjets)
        awkward_utils.add_field(syst_events, "n_fatjets", n_fatjets)

        boosted = (n_fatjets >= 1) 
        hadronic = (n_jets >= 2) & (n_leptons == 0) & (n_taus == 0)
        semi_lep = (n_jets >= 2) & (n_leptons == 1) & (n_taus == 0)
        semi_tau = (n_jets >= 2) & (n_leptons == 0) & (n_taus == 1)
        dilepton = (n_jets >= 1) & (n_lep_tau == 2)
        multilep = (n_jets >= 1) & (n_lep_tau >= 3)

        pho_idmva_cut = (syst_events.LeadPhoton.mvaID > -0.7) & (syst_events.SubleadPhoton.mvaID > -0.7)

        presel_cut = (boosted | hadronic | semi_lep | semi_tau | dilepton | multilep) & z_veto & pho_idmva_cut

        self.register_cuts(
            names = ["boosted", "hadronic", "semi_lep", "semi_tau", "dilepton", "multilepton", "z_veto", "pho_idmva_cut", "inclusive"],
            results = [boosted, hadronic, semi_lep, semi_tau, dilepton, multilep, z_veto, pho_idmva_cut, presel_cut]
        )

        return presel_cut, syst_events
    



