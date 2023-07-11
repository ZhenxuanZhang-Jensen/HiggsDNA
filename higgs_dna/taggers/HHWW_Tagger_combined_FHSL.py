import logging

import awkward
import vector
import numpy
from higgs_dna.selections import (fatjet_selections, jet_selections,
                                  lepton_selections,gen_selections)
from higgs_dna.taggers.tagger import NOMINAL_TAG, Tagger
from higgs_dna.utils import awkward_utils, misc_utils

vector.register_awkward()

def delta_R(objects1, objects2, max_dr):
    """
    Select objects from objects1 which are at least max_dr away from all objects in objects2.
    :param objects1: objects which are required to be at least max_dr away from all objects in objects2 
    :type objects1: awkward.highlevel.Array
    :param objects2: objects which are all objects in objects1 must be at leats max_dr away from
    :type objects2: awkward.highlevel.Array
    :param max_dr: minimum delta R between objects
    :type max_dr: float
    :return: boolean array of objects in objects1 which pass delta_R requirement
    :rtype: awkward.highlevel.Array
    """
    #obj1:jet obj2:muon     
    # if awkward.count(objects1) == 0 or awkward.count(objects2) == 0:
    #     return objects1.pt < 0. 
    if awkward.count(objects1) == 0: 
        return objects1.pt < 0. 
    if awkward.count(objects2) == 0:
        return objects1.pt > 0.

    if not isinstance(objects1, vector.Vector4D):
        objects1 = awkward.Array(objects1, with_name = "Momentum4D")
    if not isinstance(objects2, vector.Vector4D):
        objects2 = awkward.Array(objects2, with_name = "Momentum4D")

    obj1 = awkward.unflatten(objects1, counts = 1, axis = -1) # shape [n_events, n_obj, 1]
    obj2 = awkward.unflatten(objects2, counts = 1, axis = 0) # shape [n_events, 1, n_obj]

    dR = obj1.deltaR(obj2) # shape [n_events, n_obj1, n_obj2]

    selection = awkward.all(numpy.abs(dR) <= max_dr, axis = -1)
    return selection
    
logger = logging.getLogger(__name__)


DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {
    "electrons": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id": "WP80iso_WP90noniso",
    },
    "electrons_iso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id": "WP90",
    },
    "muons": {
        "pt" : 10.0, 
        "dr_photons": 0.4,
        "id" : "tight",
        "non_pfRelIso04_all":0.15,
        "non_iso": None,
        "global" : True
        },
    "muons_iso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id":"tight",
        "dr_photons": 0.4,
        "pfRelIso04_all":0.15,
        "iso" : None,
    },
    "jets": {
        "pt": 20.0, # attention this is the one exact same as old framework, make this 20 GeV(loose) for further analysis, we all know the higgs-like ak4 jets pt can be very small
        "eta": 2.4,
        "dr_photons": 0.4,
        "dr_electrons": 0.4,
        "dr_muons": 0.4,
        "dr_jets": 0.4,
    },
    "fatjets": {
        "pt": 200.0,
        "eta": 2.4,
        "Hqqqq_vsQCDTop": 0.4,
        "dr_photons": 0.8,
        "dr_electrons": 0.8,
        "dr_muons": 0.8
    },
    "fatjets_H": {
        "pt": 300,
        "eta": 2.4,
        # "Hqqqq_qqlv_vsQCDTop" :0.2,
        "Hqqqq_vsQCDTop" :0.4,
        "dr_photons": 0.8,
        "dr_electrons": 0.8,
        "dr_muons": 0.8
    },
    "photon_id": -0.7,
    "btag_wp": {
        "2016": 0.3093,
        "2017": 0.3040,
        "2018": 0.2783
    },
    "gen_info" : {
        "is_Signal" : False, #attention: change in HHWW_preselection.json
    }  
}


class HHWW_Preselection_FHSL(Tagger):
    """
    HHWW Preselection tagger for tutorial
    """

    def __init__(self, name, options={}, is_data=None, year=None, output_dir=None):
        super(HHWW_Preselection_FHSL, self).__init__(name, options, is_data, year, output_dir)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                original=DEFAULT_OPTIONS,
                new=options
            )

    def calculate_selection(self, events):
        # Gen selection
        # data will not select gen level infos 
        # need to comment when run bkgs
        # logger.debug("Is Signal: %s" %self.options["gen_info"]["is_Signal"])
        if not self.is_data and self.options["gen_info"]["is_Signal"]:    
            # fake_pho,prompt_pho = gen_selections.gen_Hww_4q(events)        
            gen_l1_p4, gen_q1_p4,gen_q2_p4 = gen_selections.gen_Hww_2q2l(events)        
        logger.debug("event fields: %s" %events.fields)
        original_electrons = awkward_utils.add_field(
            events=events,
            name="OriginalElectron",
            data=events.Electron)
        events["electron_pt"] = original_electrons.pt
        events["electron_eta"] = original_electrons.eta
        events["electron_phi"] = original_electrons.phi
        events["electron_mass"] = original_electrons.mass
        events["electron_dxy"] = original_electrons.dxy
        events["electron_dz"] = original_electrons.dz
        events["electron_mvaFall17V2Iso_WP90"] = original_electrons.mvaFall17V2Iso_WP90
        events["electron_mvaFall17V2noIso_WP90"] = original_electrons.mvaFall17V2noIso_WP90
        events["electron_mvaFall17V2Iso_WPL"] = original_electrons.mvaFall17V2Iso_WPL
        events["electron_mvaFall17V2noIso_WPL"] = original_electrons.mvaFall17V2noIso_WPL
        events["electron_mvaFall17V2Iso_WP80"] = original_electrons.mvaFall17V2Iso_WP80
        events["electron_mvaFall17V2noIso_WP80"] = original_electrons.mvaFall17V2noIso_WP80
        events["electron_pfRelIso03_all"] = original_electrons.pfRelIso03_all
        events["electron_pfRelIso03_chg"] = original_electrons.pfRelIso03_chg
        events["electron_hoe"] = original_electrons.hoe
        
        # Electrons
        electron_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons"]["dr_photons"]
                }
            },
            name="SelectedElectron",
            tagger=self
        )

        electrons = awkward_utils.add_field(
            events=events,
            name="SelectedElectron",
            data=events.Electron[electron_cut]
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "electron_noniso",
            objects = electrons[awkward.argsort(electrons.pt, axis=-1, ascending=False)],
            n_objects = 2,
            dummy_value = -999
        )
        electron_iso_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons_iso"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons"]["dr_photons"]
                }
            },
            name="SelectedElectron_iso",
            tagger=self
        )

        electrons_iso = awkward_utils.add_field(
            events=events,
            name="SelectedElectron_iso",
            data=events.Electron[electron_iso_cut]
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "electron_iso",
            objects = electrons_iso,
            n_objects = 1,
            dummy_value = -999
        )
        original_muons = awkward_utils.add_field(
            events=events,
            name="OriginalMuon",
            data=events.Muon)
        events["muon_pt"] = original_muons.pt
        events["muon_eta"] = original_muons.eta
        events["muon_phi"] = original_muons.phi
        events["muon_mass"] = original_muons.mass
        events["muon_dxy"] = original_muons.dxy
        events["muon_dz"] = original_muons.dz
        events["muon_pfRelIso04_all"] = original_muons.pfRelIso04_all
        events["muon_highPtId"] = original_muons.highPtId
        events["muon_looseId"] = original_muons.looseId
        events["muon_mediumId"] = original_muons.mediumId
        events["muon_tightId"] = original_muons.tightId
        events["muon_isGlobal"] = original_muons.isGlobal
        # Muons
        muon_cut = lepton_selections.select_muons(
            muons=events.Muon,
            options=self.options["muons"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["muons"]["dr_photons"]
                }
            },
            name="SelectedMuon",
            tagger=self
        )

        muons = awkward_utils.add_field(
            events=events,
            name="SelectedMuon",
            data=events.Muon[muon_cut]
        )
        

        awkward_utils.add_object_fields(
            events =events,
            name = "muon_noniso",
            objects = muons[awkward.argsort(muons.pt, ascending=False, axis=-1)],
            n_objects = 1,
            dummy_value = -999
        )
        # isolated muon where we applied pfRelIso03_all
        muon_iso_cut = lepton_selections.select_muons(
            muons=events.Muon,
            options=self.options["muons_iso"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["muons_iso"]["dr_photons"]
                }
            },
            name="SelectedMuon_iso",
            tagger=self
        )

        muons_iso = awkward_utils.add_field(
            events=events,
            name="SelectedMuon_iso",
            data=events.Muon[muon_iso_cut]
        )
        

        awkward_utils.add_object_fields(
            events =events,
            name = "muon_iso",
            objects = muons_iso,
            n_objects = 1,
            dummy_value = -999
        )

        e_4p = vector.obj(
            pt=events.electron_iso_pt,
            eta=events.electron_iso_eta,
            phi=events.electron_iso_phi,
            mass=events.electron_iso_mass
        )
        leadpho_4p = vector.obj(
            pt=events.LeadPhoton.pt,
            eta=events.LeadPhoton.eta,
            phi=events.LeadPhoton.phi,
            mass=events.LeadPhoton.mass
        )
        e_leadphoton=e_4p+leadpho_4p
        e_leadphoton.mass = numpy.nan_to_num(e_leadphoton.mass,nan=-999)
        Z_veto_cut = abs(e_leadphoton.mass-91.5)>5

        # Jets
        jet_cut = jet_selections.select_jets(
            jets=events.Jet,
            options=self.options["jets"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["jets"]["dr_photons"]
                },
                "electrons": {
                    "objects": events.SelectedElectron,
                    "min_dr": self.options["jets"]["dr_electrons"]
                },
                "muons": {
                    "objects": events.SelectedMuon,
                    "min_dr": self.options["jets"]["dr_muons"]
                }
            },
            name = "SelectedJet",
            tagger=self
        )
        jets = awkward_utils.add_field(
            events=events,
            name="SelectedJet",
            data=events.Jet[jet_cut]
        )

        # --------------------- if fatjet branches are not empty --------------------- #
        # Fat jets

        # add the H jet tagger for SL&FH channel((H3q+H4q+Hlvqq)/(H3q+H4q+Hlvqq+QCD+Top))
        # ----------------------------------- PN 4q ---------------------------------- #
        PN_sigs_4q = events.FatJet.inclParTMDV1_probHWqqWq0c+events.FatJet.inclParTMDV1_probHWqqWq1c+events.FatJet.inclParTMDV1_probHWqqWq2c+events.FatJet.inclParTMDV1_probHWqqWqq0c+events.FatJet.inclParTMDV1_probHWqqWqq1c+events.FatJet.inclParTMDV1_probHWqqWqq2c
        PN_bkgs_4q = events.FatJet.inclParTMDV1_probQCDb+events.FatJet.inclParTMDV1_probQCDbb+events.FatJet.inclParTMDV1_probQCDc+events.FatJet.inclParTMDV1_probQCDcc+events.FatJet.inclParTMDV1_probQCDothers+events.FatJet.inclParTMDV1_probTopbWq0c+events.FatJet.inclParTMDV1_probTopbWq1c+events.FatJet.inclParTMDV1_probTopbWqq0c+events.FatJet.inclParTMDV1_probTopbWqq1c
        # ------------------------------------ -- ------------------------------------ #
        # ---------------------------------- PN lvqq --------------------------------- #
        PN_sigs_lvqq = events.FatJet.inclParTMDV1_probHWqqWev0c+events.FatJet.inclParTMDV1_probHWqqWev1c+events.FatJet.inclParTMDV1_probHWqqWmv0c+events.FatJet.inclParTMDV1_probHWqqWmv1c+events.FatJet.inclParTMDV1_probHWqqWtauev0c+events.FatJet.inclParTMDV1_probHWqqWtauev1c+events.FatJet.inclParTMDV1_probHWqqWtauhv0c+events.FatJet.inclParTMDV1_probHWqqWtauhv1c+events.FatJet.inclParTMDV1_probHWqqWtaumv0c+events.FatJet.inclParTMDV1_probHWqqWtaumv1c
        PN_bkgs_lvqq = events.FatJet.inclParTMDV1_probQCDb+events.FatJet.inclParTMDV1_probQCDbb+events.FatJet.inclParTMDV1_probQCDc+events.FatJet.inclParTMDV1_probQCDcc+events.FatJet.inclParTMDV1_probQCDothers+events.FatJet.inclParTMDV1_probTopbWev+events.FatJet.inclParTMDV1_probTopbWmv+events.FatJet.inclParTMDV1_probTopbWtauev+events.FatJet.inclParTMDV1_probTopbWtauhv+events.FatJet.inclParTMDV1_probTopbWtaumv
        # ------------------------------------ -- ------------------------------------ #
        # -------------------------------- PN 4q lvqq -------------------------------- #
        PN_sigs_4q_lvqq = events.FatJet.inclParTMDV1_probHWqqWev0c+events.FatJet.inclParTMDV1_probHWqqWev1c+events.FatJet.inclParTMDV1_probHWqqWmv0c+events.FatJet.inclParTMDV1_probHWqqWmv1c+events.FatJet.inclParTMDV1_probHWqqWq0c+events.FatJet.inclParTMDV1_probHWqqWq1c+events.FatJet.inclParTMDV1_probHWqqWq2c+events.FatJet.inclParTMDV1_probHWqqWqq0c+events.FatJet.inclParTMDV1_probHWqqWqq1c+events.FatJet.inclParTMDV1_probHWqqWqq2c+events.FatJet.inclParTMDV1_probHWqqWtauev0c+events.FatJet.inclParTMDV1_probHWqqWtauev1c+events.FatJet.inclParTMDV1_probHWqqWtauhv0c+events.FatJet.inclParTMDV1_probHWqqWtauhv1c+events.FatJet.inclParTMDV1_probHWqqWtaumv0c+events.FatJet.inclParTMDV1_probHWqqWtaumv1c
        PN_bkgs_4q_lvqq = events.FatJet.inclParTMDV1_probQCDb+events.FatJet.inclParTMDV1_probQCDbb+events.FatJet.inclParTMDV1_probQCDc+events.FatJet.inclParTMDV1_probQCDcc+events.FatJet.inclParTMDV1_probQCDothers+events.FatJet.inclParTMDV1_probTopbWev+events.FatJet.inclParTMDV1_probTopbWmv+events.FatJet.inclParTMDV1_probTopbWq0c+events.FatJet.inclParTMDV1_probTopbWq1c+events.FatJet.inclParTMDV1_probTopbWqq0c+events.FatJet.inclParTMDV1_probTopbWqq1c+events.FatJet.inclParTMDV1_probTopbWtauev+events.FatJet.inclParTMDV1_probTopbWtauhv+events.FatJet.inclParTMDV1_probTopbWtaumv
        # ------------------------------------ -- ------------------------------------ #
        
        fatjet_tmp = events.FatJet
        fatjet_tmp['Hqqqq_qqlv_vsQCDTop'] = PN_sigs_4q_lvqq / (PN_bkgs_4q_lvqq + PN_sigs_4q_lvqq)
        fatjet_tmp['Hqqqq_vsQCDTop'] = PN_sigs_4q / (PN_bkgs_4q + PN_sigs_4q)
        fatjet_tmp['Hlvqq_vsQCDTop'] = PN_sigs_lvqq / (PN_bkgs_lvqq + PN_sigs_lvqq)
        
        events.FatJet = fatjet_tmp
        print("events.FatJet exsiting field:\n",dir(events.FatJet))
        print("events.FatJet tyep is: \n", type(events.FatJet))
        # ------------------------------------- - ------------------------------------ #
        fatjet_cut = fatjet_selections.select_fatjets(
            fatjets = events.FatJet,
            options = self.options["fatjets"],
            clean = {
                "photons" : {
                    "objects" : events.Diphoton.Photon,
                    "min_dr" : self.options["fatjets"]["dr_photons"]
                },
                "electrons" : {
                    "objects" : events.SelectedElectron_iso,
                    "min_dr" : self.options["fatjets"]["dr_electrons"]
                },
                "muons" : {
                    "objects" : events.SelectedMuon_iso,
                    "min_dr" : self.options["fatjets"]["dr_muons"]
                    },
                },
            name = "SelectedFatJet",
            tagger = self
        )
        fatjets = awkward_utils.add_field(
            events = events,
            name = "SelectedFatJet",
            data = events.FatJet[fatjet_cut]
        )   
        events['SlectedFatJet_pt']=fatjets.pt
        events['SlectedFatJet_eta']=fatjets.eta
        events['SlectedFatJet_phi']=fatjets.phi
        events['SlectedFatJet_mass']=fatjets.mass
        events['SlectedFatJet_Hqqqq_vsQCDTop']=fatjets.Hqqqq_vsQCDTop
        events['SlectedFatJet_particleNet_WvsQCD']=fatjets.particleNet_WvsQCD

        awkward_utils.add_object_fields(
        events=events,
        name="fatjet",
        objects=fatjets[awkward.argsort(fatjets.pt, ascending=False, axis=-1)],
        n_objects=3,
        dummy_value=-999
        )
        fatjet_H_cut = fatjet_selections.select_fatjets(
            fatjets = events.FatJet,
            options = self.options["fatjets_H"],
            clean = {
                "photons" : {
                    "objects" : events.Diphoton.Photon,
                    "min_dr" : self.options["fatjets_H"]["dr_photons"]
                },
                "electrons" : {
                    "objects" : events.SelectedElectron_iso,
                    "min_dr" : self.options["fatjets_H"]["dr_electrons"]
                },
                "muons" : {
                    "objects" : events.SelectedMuon_iso,
                    "min_dr" : self.options["fatjets_H"]["dr_muons"]
                    },
                },
            name = "SelectedFatJet_H",
            tagger = self
        )
        fatjets_H = awkward_utils.add_field(
            events = events,
            name = "SelectedFatJet_H",
            data = events.FatJet[fatjet_H_cut]
        )   

        awkward_utils.add_object_fields(
        events=events,
        name="fatjet_H",
        # objects=fatjets_H[awkward.argsort(fatjets_H.Hqqqq_qqlv_vsQCDTop, ascending=False, axis=-1)],
        objects=fatjets_H[awkward.argsort(fatjets_H.Hqqqq_vsQCDTop, ascending=False, axis=-1)],
        n_objects=1,
        dummy_value=-999
        ) 
        
        # produce the Hww3q tagger
        # Hww3qvsQCD = (fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq0c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq1c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq2c) / (fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq0c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq1c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq2c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDb+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDbb+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDc+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDcc+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDothers)



        # gen 4q deltaR with j1,j2,j3,j4
        # if not self.is_data and self.options["gen_info"]["is_Signal"]:    
            # print("debuggg")
            # gen_q1_p4,gen_q2_p4,gen_q3_p4,gen_q4_p4=gen_selections.gen_Hww_4q(events)
        # jet_p4 = vector.awk(
        #     {
        #         "pt" : jets["pt"],
        #         "eta" : jets["eta"],
        #         "phi" : jets["phi"],
        #         "mass" : jets["mass"]
        #     },
        #     with_name = "Momentum4D"
        # )

        # if not self.is_data and self.options["gen_info"]["is_Signal"]:    
        #     jets["deltaR_q1"] = jet_p4.deltaR(gen_q1_p4)
        #     jets["deltaR_q2"] = jet_p4.deltaR(gen_q2_p4)
        #     jets["deltaR_q3"] = jet_p4.deltaR(gen_q3_p4)
        #     jets["deltaR_q4"] = jet_p4.deltaR(gen_q4_p4)

        awkward_utils.add_object_fields(
            events=events,
            name="jet",
            objects=jets,
            n_objects=7,
            dummy_value=-999
        )
        electrons = awkward.Array(electrons, with_name="Momentum4D")
        electrons_iso = awkward.Array(electrons_iso, with_name="Momentum4D")
        muons = awkward.Array(muons, with_name="Momentum4D")
        muons_iso = awkward.Array(muons_iso, with_name="Momentum4D")
        # bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis=1, ascending=False)]
        # bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]]

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        # Preselection
        
        n_electrons = awkward.num(electrons)
        n_electrons_iso = awkward.num(electrons_iso)
        n_muons = awkward.num(muons)
        n_muons_iso = awkward.num(muons_iso)
        n_leptons = n_electrons + n_muons 
        n_leptons_iso = n_electrons_iso + n_muons_iso
        awkward_utils.add_field(events,"nGoodnonisoelectrons",n_electrons)
        awkward_utils.add_field(events,"nGoodisoelectrons",n_electrons_iso)
        awkward_utils.add_field(events,"nGoodnonisomuons",n_muons)
        awkward_utils.add_field(events,"nGoodisomuons",n_muons_iso)
        awkward_utils.add_field(events,"nGoodisoleptons",n_leptons_iso)
        awkward_utils.add_field(events,"nGoodnonisoleptons",n_leptons)
        # n_diphotons = awkward.num(events.Diphoton)
        # logger.debug(" the N_diphoton : %f" % (n_diphotons))
        n_jets = awkward.num(jets)
        awkward_utils.add_field(events,"nGoodAK4jets",n_jets)
        n_fatjets = awkward.num(fatjets)
        n_fatjets_H = awkward.num(fatjets_H)
        awkward_utils.add_field(events,"nGoodAK8jets",n_fatjets)
        # awkward_utils.add_field(events,"nGood_W_fatjets",n_fatjets_W)
        awkward_utils.add_field(events,"nGood_H_fatjets",n_fatjets_H)
        # n_bjets = awkward.num(bjets)

        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (
            events.SubleadPhoton.mvaID > self.options["photon_id"])

        # If have isolated lepton
        #attention: Semi leptonic channel                
        # add leptonic boosted category with (>=1 AK8 jets && WvsQCD > 0.94)
        # the the first fatjet with WvsQCD > 0.94
        selection_fatjet_WvsQCD_SL_cat0 = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)]) >= 1
        SL_boosted_cat = (n_leptons_iso == 1) & (n_fatjets >=1) & (selection_fatjet_WvsQCD_SL_cat0) # boosted 1 jet for SL channel with isolated lep
        # add Leptonic resolved category with (>=2 AK4 jets && WvsQCD < 0.5)
        # need the first fatjets with WvsQCD < 0.5
        selection_fatjet_WvsQCD_SL_cat1 = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)]) == 0
        SL_fullyresovled_cat = (n_leptons_iso == 1) & (n_jets >=2) & (selection_fatjet_WvsQCD_SL_cat1) # resolved 2 jets for SL channel with isolated lep


        # If no isolated lepton
        #attention: Fully leptonic channel

        # add boosted FH category with (>=1 Higgs jets && HvsQCD > 0.2)
        # the first Higgs jet with HvsQCD > 0.8
        selection_fatjet_HvsQCD_FH_cat0 = awkward.num(fatjets_H.Hqqqq_vsQCDTop[(fatjets_H.Hqqqq_vsQCDTop > 0.4)]) >= 1
        FH_boosted = (n_leptons_iso == 0) & (n_fatjets_H >=1) & (selection_fatjet_HvsQCD_FH_cat0) # boosted 1 jet for SL and FH channel wo isolated lep

        # add semi-boosted FH -1 category with (>=2 AK8 jets && WvsQCD > 0.5)
        # need the first two fatjets with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SB_2F = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)]) >= 2
        FH_SB_2Fatjet = (n_leptons_iso==0) & (n_fatjets >=2) & (selection_fatjet_WvsQCD_SB_2F) # 2 jets for FH

        # add semi-boosted FH -2 category with (==1 AK8 jets && WvsQCD > 0.9 && >=2 AK4 jets)
        # need the first fatjet with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SB_1F = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)]) >= 1
        FH_SB_1Fatjet = (n_leptons_iso==0) & (n_fatjets ==1) & (n_jets >=2) & (selection_fatjet_WvsQCD_SB_1F) # 1 jet for FH
        # add merged lepton boosted category with (==1 AK8 jets && WvsQCD < 0.9 && HvsQCD < 0.8 && btagDeepB < 0.4184 && pt > 300)
        selection_fatjet_fail_WvsQCD = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD < 0.94)]) >= 1
        selection_fatjet_fail_HvsQCD = awkward.num(fatjets.Hqqqq_vsQCDTop[(fatjets.Hqqqq_vsQCDTop < 0.4)]) >= 1
        # fatjet_b_veto=awkward.num(fatjets.btagDeepB[(fatjets.btagDeepB > 0.4184)]) == 0
        SL_Higgs_jet_pt=awkward.num(fatjets.pt[(fatjets.pt > 300)])
        SL_lep_merge_boosted_cat = ((n_leptons_iso==0)&(n_leptons>=0))&selection_fatjet_fail_WvsQCD & selection_fatjet_fail_HvsQCD & SL_Higgs_jet_pt
        # add resolved FH category with (>=4 AK4 jets )
        FH_fully_resovled_cat = (n_leptons_iso==0) & (n_jets>=4) # 4 jets for FH

        #attention: Semi leptonic channel
        # add muon merged boosted category with (>=1 AK8 jets && WvsQCD > 0.5, and deltaR(muon, fatjet) < 0.8))
        # selection_fatjet_WvsQCD = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.9)]) >= 1
        # muon deltaR with fatjet < 0.8
        # find the largest particleNet_WvsQCD fatjets
        # fatjets_leading_patricleNet_WvsQCD = fatjets[awkward.argmax(fatjets.particleNet_WvsQCD, axis=1, keepdims=True)]
        # muon_fatjet_deltaR = delta_R(muons, fatjets_leading_patricleNet_WvsQCD, 0.8)
        # selection_muon_fatjet_deltaR = awkward.num(muon_fatjet_deltaR[muon_fatjet_deltaR==True]) >= 1 # at least one muon deltaR with all fatjets < 0.8
        # one non-isolated muon, zero isolated lepton
        # SL_muon_merge_boosted_cat = (n_muons >=1) & (n_leptons_iso ==0) & (selection_muon_fatjet_deltaR) & (selection_fatjet_WvsQCD) # boosted 1 jet for SL channel with merge non-iso muon

        # add muon merged fully resolved category with (>=2 AK4 jets && WvsQCD > 0.5 && deltaR(lepton, jet) < 0.4)
        # deltaR(lepton, jet) < 0.4
        # lepton_jet_deltaR = delta_R(muons, jets, 0.4)
        # selection_lepton_jet_deltaR = awkward.num(lepton_jet_deltaR[lepton_jet_deltaR==True]) >= 1 # at least one lepton deltaR with all jets < 0.4
        # SL_muon_merge_full_resolved_cat = (n_muons >=1) & (n_leptons_iso ==0) & (n_jets >=2) & (selection_lepton_jet_deltaR)# resolved 2 jets for SL channel with merge non-iso muon


        # Hadronic presel
        # use priority to mark different category
        flatten_n_jets = awkward.num(jets.pt)
        category = awkward.zeros_like(flatten_n_jets)
        category = awkward.fill_none(category, 0)
        # add the priority for each category
        # 1: FH_fully_resovled_cat
        # 2: FH_SB_1Fatjet
        # 3: FH_SB_2Fatjet
        # 4: FH_boosted
        # 5: SL_muon_merge_full_resolved_cat
        # 6: SL_muon_merge_boosted_cat
        # 7: SL_fullyresovled_cat
        # 8: SL_boosted_cat
        # with no isolated lepton
        category = awkward.where(FH_fully_resovled_cat, awkward.ones_like(category)*1, category)
        category = awkward.where(FH_SB_1Fatjet, awkward.ones_like(category)*2, category)
        category = awkward.where(FH_SB_2Fatjet, awkward.ones_like(category)*3, category)
        category = awkward.where(FH_boosted, awkward.ones_like(category)*4, category)
        category = awkward.where(SL_lep_merge_boosted_cat, awkward.ones_like(category)*5, category)
        # category = awkward.where(SL_muon_merge_full_resolved_cat, awkward.ones_like(category)*5, category)

        # category = awkward.where(SL_muon_merge_boosted_cat, awkward.ones_like(category)*6, category)
        # with isolated lepton
        category = awkward.where(SL_fullyresovled_cat, awkward.ones_like(category)*6, category)
        category = awkward.where(SL_boosted_cat, awkward.ones_like(category)*7, category)
        category_cut = (category >= 0) # cut the events with category == 0
        awkward_utils.add_field(events, "category", category) 

        presel_cut = (photon_id_cut) & (category_cut) & (Z_veto_cut)

        self.register_cuts(
            names=["Photon id Selection","category_cut", "Z_veto_cut"],
            results=[photon_id_cut, category_cut, Z_veto_cut]
        )
        return presel_cut, events