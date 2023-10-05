import logging

import awkward
import vector
import numpy
from higgs_dna.selections import (fatjet_selections, jet_selections,
                                  lepton_selections,gen_selections)
from higgs_dna.taggers.tagger import NOMINAL_TAG, Tagger
from higgs_dna.utils import awkward_utils, misc_utils

vector.register_awkward()

logger = logging.getLogger(__name__)


DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {
    "electrons": {
        "pt": 10.0,
        "dr_photons": 0.4
    },
    "muons": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "pfRelIso04_all" : 0.15,
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
        "Hqqqq_qqlv_vsQCDTop": -999,
        "dr_photons": 0.8,
        "dr_electrons": 0.8,
        "dr_muons": 0.8
    },
    "fatjets_H": {
        "pt": 300,#300.0 fixed 0 
        "eta": 2.4,
        "Hqqqq_qqlv_vsQCDTop" :0.2,
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
        super(HHWW_Preselection_FHSL, self).__init__(name, options, is_data, year,output_dir)

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
        # if not self.is_data and self.options["gen_info"]["is_Signal"]:    
        #    gen_selections.gen_Hww_4q(events)        
        logger.debug("event fields: %s" %events.fields)

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
            name = "electron",
            objects = electrons,
            n_objects = 1,
            dummy_value = -999
        )
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
            name = "muon",
            objects = muons,
            n_objects = 1,
            dummy_value = -999
        )

        
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
                    "objects" : events.SelectedElectron,
                    "min_dr" : self.options["fatjets"]["dr_electrons"]
                },
                "muons" : {
                    "objects" : events.SelectedMuon,
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


        awkward_utils.add_object_fields(
        events=events,
        name="fatjet",
        objects=fatjets[awkward.argsort(fatjets.particleNet_WvsQCD, ascending=False, axis=-1)], # sort by particleNet_WvsQCD
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
                    "objects" : events.SelectedElectron,
                    "min_dr" : self.options["fatjets_H"]["dr_electrons"]
                },
                "muons" : {
                    "objects" : events.SelectedMuon,
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
        objects=fatjets_H[awkward.argsort(fatjets_H.Hqqqq_qqlv_vsQCDTop, ascending=False, axis=-1)],
        n_objects=3,
        dummy_value=-999
        ) 
        
        # produce the Hww3q tagger
        # Hww3qvsQCD = (fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq0c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq1c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq2c) / (fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq0c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq1c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probHWqqWqq2c + fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDb+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDbb+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDc+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDcc+fatjets[fatjet_W_cut].FatJet_inclParTMDV1_probQCDothers)



        # gen 4q deltaR with j1,j2,j3,j4
        # if not self.is_data and self.options["gen_info"]["is_Signal"]:    
        #     print("debuggg")
        #     gen_q1_p4,gen_q2_p4,gen_q3_p4,gen_q4_p4=gen_selections.gen_Hww_4q(events)
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
        muons = awkward.Array(muons, with_name="Momentum4D")
        # bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis=1, ascending=False)]
        # bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]]

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        # Preselection
        
        n_electrons = awkward.num(electrons)
        n_muons = awkward.num(muons)
        n_leptons = n_electrons + n_muons

        # n_diphotons = awkward.num(events.Diphoton)
        # logger.debug(" the N_diphoton : %f" % (n_diphotons))
        n_jets = awkward.num(jets)
        awkward_utils.add_field(events,"nGoodAK4jets",n_jets)
        n_fatjets = awkward.num(fatjets)
        n_fatjets_H = awkward.num(fatjets_H)
        awkward_utils.add_field(events,"nGoodAK4jets",n_jets)
        # awkward_utils.add_field(events,"nGood_W_fatjets",n_fatjets_W)
        awkward_utils.add_field(events,"nGood_H_fatjets",n_fatjets_H)
        # n_bjets = awkward.num(bjets)

        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (
            events.SubleadPhoton.mvaID > self.options["photon_id"])

        # ----- If isolated lepton ------------
        # add leptonic boosted category with (>=1 AK8 jets && WvsQCD > 0.5)
        # the the first fatjet with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SL_cat0 = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.5)]) >= 1
        SL_cat0 = (n_leptons == 1) & (n_fatjets >=1) & (selection_fatjet_WvsQCD_SL_cat0) # boosted 1 jet for SL channel with isolated lep
        # add Leptonic resolved category with (>=2 AK4 jets && WvsQCD < 0.5)
        # need the first fatjets with WvsQCD < 0.5
        selection_fatjet_WvsQCD_SL_cat1 = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.5)]) == 0
        SL_cat1 = (n_leptons == 1) & (n_jets >=2) & (selection_fatjet_WvsQCD_SL_cat1) # resolved 2 jets for SL channel with isolated lep
        
        # ----------  if no isolated lepton -------------

        # add boosted SL+FH category with (>=1 Higgs jets && HvsQCD > 0.2)
        # the first Higgs jet with HvsQCD > 0.2
        selection_fatjet_HvsQCD_SL_FH_cat0 = awkward.num(fatjets_H.Hqqqq_qqlv_vsQCDTop[(fatjets_H.Hqqqq_qqlv_vsQCDTop > 0.2)]) >= 1
        SL_FH_cat0 = (n_leptons == 0) & (n_fatjets_H >=1) & (selection_fatjet_HvsQCD_SL_FH_cat0) # boosted 1 jet for SL and FH channel wo isolated lep

        # add semi-boosted FH -1 category with (>=2 AK8 jets && WvsQCD > 0.5)
        # need the first two fatjets with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SB_2F = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.5)]) >= 2
        FH_cat0_SB_2F = (n_leptons==0) & (n_fatjets >=2) & (selection_fatjet_WvsQCD_SB_2F) # 2 jets for FH

        # add semi-boosted FH -2 category with (==1 AK8 jets && WvsQCD > 0.5 && >=2 AK4 jets)
        # need the first fatjet with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SB_1F = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.5)]) >= 1
        FH_cat0_SB_1F = (n_leptons==0) & (n_fatjets ==1) & (n_jets >=2) & (selection_fatjet_WvsQCD_SB_1F) # 1 jet for FH

        # add resolved FH category with (>=4 AK4 jets )
        FH_cat1 = (n_leptons==0) & (n_jets>=4) # 4 jets for FH

        # Hadronic presel
        # use priority to mark different category
        flatten_n_jets = awkward.num(jets.pt)
        category = awkward.zeros_like(flatten_n_jets)
        category = awkward.fill_none(category, 0)
        # add the priority for each category
        category = awkward.where(FH_cat1, awkward.ones_like(category)*6, category)
        category = awkward.where(FH_cat0_SB_1F, awkward.ones_like(category)*5, category)
        category = awkward.where(FH_cat0_SB_2F, awkward.ones_like(category)*4, category)
        category = awkward.where(SL_FH_cat0, awkward.ones_like(category)*3, category)
        category = awkward.where(SL_cat1, awkward.ones_like(category)*2, category)
        category = awkward.where(SL_cat0, awkward.ones_like(category)*1, category)

        category_cut = (category >= 0)
        awkward_utils.add_field(events, "category", category) 

        presel_cut = (photon_id_cut) & (category_cut)

        self.register_cuts(
            names=["Photon id Selection","category_cut"],
            results=[photon_id_cut, category_cut]
        )
        return presel_cut, events