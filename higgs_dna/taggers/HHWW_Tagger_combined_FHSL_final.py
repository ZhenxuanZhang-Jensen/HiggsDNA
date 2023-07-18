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
    "electrons_noiso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id": "WP80iso_WP90noniso",
    },
    "electrons_iso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id": "WP80",
    },
    "muons_noiso": {
        "pt" : 10.0, 
        "dr_photons": 0.4,
        "id" : "highptId",
        "non_pfRelIso04_all":0.15,
         },
    "muons_iso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id":"highptId",
        "dr_photons": 0.4,
        "pfRelIso04_all":0.15,
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
        # "Hqqqq_vsQCDTop": 0.4,
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
        # if not self.is_data and self.options["gen_info"]["is_Signal"]:    
        # fake_pho,prompt_pho = gen_selections.gen_Hww_4q(events)        
            # gen_l1_p4, gen_q1_p4,gen_q2_p4 = gen_selections.gen_Hww_2q2l(events)        
        logger.debug("event fields: %s" %events.fields)
        # logger.debug('After Diphoton selection')

        # Noiso-Electrons
        electron_noiso_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons_noiso"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons_noiso"]["dr_photons"]
                }
            },
            name="SelectedElectron_noiso",
            tagger=self
        )

        electrons_noiso = awkward_utils.add_field(
            events=events,
            name="SelectedElectron_noiso",
            data=events.Electron[electron_noiso_cut]
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "electron_noniso",
            objects = electrons_noiso[awkward.argsort(electrons_noiso.pt, axis=-1, ascending=False)],
            n_objects = 2,
            dummy_value = -999
        )
        # Iso-Electrons
        electron_iso_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons_iso"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons_noiso"]["dr_photons"]
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
        # Noiso-Muons
        muon_noiso_cut = lepton_selections.select_muons(
            muons=events.Muon,
            options=self.options["muons_noiso"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["muons_noiso"]["dr_photons"]
                }
            },
            name="SelectedMuon_noiso",
            tagger=self
        )

        muons_noiso = awkward_utils.add_field(
            events=events,
            name="SelectedMuon_noiso",
            data=events.Muon[muon_noiso_cut]
        )
        

        awkward_utils.add_object_fields(
            events =events,
            name = "muon_noniso",
            objects = muons_noiso[awkward.argsort(muons_noiso.pt, ascending=False, axis=-1)],
            n_objects = 1,
            dummy_value = -999
        )
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
                    "objects": events.SelectedElectron_iso,
                    "min_dr": self.options["jets"]["dr_electrons"]
                },
                "muons": {
                    "objects": events.SelectedMuon_iso,
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
        # ----------------------------------- PN 3q4q ---------------------------------- #
        PN_sigs_4q = events.FatJet.inclParTMDV1_probHWqqWq0c+events.FatJet.inclParTMDV1_probHWqqWq1c+events.FatJet.inclParTMDV1_probHWqqWq2c+events.FatJet.inclParTMDV1_probHWqqWqq0c+events.FatJet.inclParTMDV1_probHWqqWqq1c+events.FatJet.inclParTMDV1_probHWqqWqq2c
        PN_bkgs_4q = events.FatJet.inclParTMDV1_probQCDb+events.FatJet.inclParTMDV1_probQCDbb+events.FatJet.inclParTMDV1_probQCDc+events.FatJet.inclParTMDV1_probQCDcc+events.FatJet.inclParTMDV1_probQCDothers+events.FatJet.inclParTMDV1_probTopbWq0c+events.FatJet.inclParTMDV1_probTopbWq1c+events.FatJet.inclParTMDV1_probTopbWqq0c+events.FatJet.inclParTMDV1_probTopbWqq1c
        
        fatjet_tmp = events.FatJet
        fatjet_tmp['Hqqqq_vsQCDTop'] = PN_sigs_4q / (PN_bkgs_4q + PN_sigs_4q)        
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
        
       
        awkward_utils.add_object_fields(
            events=events,
            name="jet",
            objects=jets,
            n_objects=7,
            dummy_value=-999
        )
        electrons_noiso = awkward.Array(electrons_noiso, with_name="Momentum4D")
        electrons_iso = awkward.Array(electrons_iso, with_name="Momentum4D")
        muons_noiso = awkward.Array(muons_noiso, with_name="Momentum4D")
        muons_iso = awkward.Array(muons_iso, with_name="Momentum4D")
        lepton_noniso=awkward.concatenate([electrons_noiso,muons_noiso],axis=1)
        # bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis=1, ascending=False)]
        # bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]]

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        # Preselection
        
        n_electrons_noiso = awkward.num(electrons_noiso)
        n_electrons_iso = awkward.num(electrons_iso)
        n_muons_noiso = awkward.num(muons_noiso)
        n_muons_iso = awkward.num(muons_iso)
        n_leptons_noiso = n_electrons_noiso + n_muons_noiso 
        n_leptons_iso = n_electrons_iso + n_muons_iso
        awkward_utils.add_field(events,"nGoodnonisoelectrons",n_electrons_noiso)
        awkward_utils.add_field(events,"nGoodisoelectrons",n_electrons_iso)
        awkward_utils.add_field(events,"nGoodnonisomuons",n_muons_noiso)
        awkward_utils.add_field(events,"nGoodisomuons",n_muons_iso)
        awkward_utils.add_field(events,"nGoodisoleptons",n_leptons_iso)
        awkward_utils.add_field(events,"nGoodnonisoleptons",n_leptons_noiso)
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
        # ----------------------------------------------------------------------------------------------------#
        # If have isolated lepton
        #attention: Semi leptonic channel                
        # add leptonic boosted category with (>=1 AK8 jets && WvsQCD > 0.94)
        # the the first fatjet with WvsQCD > 0.94
        selection_fatjet_WvsQCD_SL_cat0 = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)]) >= 1
        SL_boosted_cat = (n_leptons_iso == 1) & (selection_fatjet_WvsQCD_SL_cat0) # boosted 1 jet for SL channel with isolated lep
        # add Leptonic resolved category with (>=2 AK4 jets && WvsQCD < 0.94)
        # need the first fatjets with WvsQCD < 0.5
        selection_fatjet_WvsQCD_SL_cat1 = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)]) == 0
        SL_fullyresovled_cat = (~SL_boosted_cat)&(n_leptons_iso == 1) & (n_jets >=2) & (selection_fatjet_WvsQCD_SL_cat1) # resolved 2 jets for SL channel with isolated lep
        # ----------------------------------------------------------------------------------------------------#
        # If no isolated lepton but >=1 non-isolated lepton with >=1 AK8 jets && AK8 jet with pt > 300 GeV && dR(lep, AK8) < 0.8
        #attention: Merged leptonic boosted channel
        selection_lepton_merged_SL_cat0 = awkward.num(fatjets.pt[(fatjets.pt > 300)]) >= 1
        
        # SL_merged_boosted_cat = (n_leptons_noiso == 1) & (n_leptons_iso==0) & (selection_lepton_merged_SL_cat0) & (awkward.num(dR_lep_fatjet < 0.8)==1) # merged boosted 1 jet for SL channel wo isolated lep
        dR_lep_fatjet=delta_R(fatjets,lepton_noniso,0.8)
        SL_merged_boosted_cat = (~SL_boosted_cat) & (~SL_fullyresovled_cat) & (n_leptons_noiso == 1) & (n_leptons_iso==0) & (selection_lepton_merged_SL_cat0) &(awkward.num((dR_lep_fatjet==True)>1))  # merged boosted 1 jet for SL channel wo isolated lep
        
        
        # fatjet_obj=awkward.fill_none(fatjets.mask[SL_merged_boosted_cat],dummy_fatjet,axis=0)
        # lepton_obj=awkward.fill_none(lepton_noniso.mask[SL_merged_boosted_cat],dummy_lepton,axis=0)
        # logger.debug(lepton_obj)
        # unflatten_fatjet_obj= awkward.unflatten(fatjet_obj[:,0],1)
        # fatjet_obj = awkward.Array(unflatten_fatjet_obj, with_name = "Momentum4D")
        # lepton_obj=lepton_obj[:,0]
        # unflatten_lepton_obj= awkward.unflatten(lepton_obj,1)
        # lepton_obj = awkward.Array(unflatten_lepton_obj, with_name = "Momentum4D")
        # dR_lep_fatjet=lepton_obj.deltaR(fatjet_obj)
        # SL_merged_boosted_cat = SL_merged_boosted_cat & (awkward.num(dR_lep_fatjet < 0.8)==1)
        # ----------------------------------------------------------------------------------------------------#
        # If no lepton
        #attention: Fully hadronic channel
        # add boosted FH category with (>=1 Higgs jets && HvsQCD > 0.4)
        # the first Higgs jet with HvsQCD > 0.4
        selection_fatjet_HvsQCD_FH_cat0 = awkward.num(fatjets_H.Hqqqq_vsQCDTop[(fatjets_H.Hqqqq_vsQCDTop > 0.4)]) >= 1
        FH_boosted = (~SL_boosted_cat) & (~SL_fullyresovled_cat) & (~SL_merged_boosted_cat) & (n_leptons_iso == 0) & (n_leptons_noiso == 0) & (n_fatjets_H >=1) & (selection_fatjet_HvsQCD_FH_cat0) # boosted 1 jet for SL and FH channel wo isolated lep
        # ----------------------------------------------------------------------------------------------------#
        # add semi-boosted FH -1 category with (>=2 AK8 jets && WvsQCD > 0.94)
        # need the first two fatjets with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SB_2F = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)&(fatjets.Hqqqq_vsQCDTop<=0.4)]) >= 2
        FH_2Wfatjet_cat = (~SL_boosted_cat) & (~SL_fullyresovled_cat) & (~SL_merged_boosted_cat) & (~FH_boosted) & (n_leptons_iso==0) & (n_leptons_noiso == 0) & (n_fatjets >=2) & (selection_fatjet_WvsQCD_SB_2F) # 2 jets for FH
        # ----------------------------------------------------------------------------------------------------#
        # add semi-boosted FH -2 category with (==1 AK8 jets && WvsQCD > 0.94 && >=2 AK4 jets)
        # need the first fatjet with WvsQCD > 0.5
        selection_fatjet_WvsQCD_SB_1F = awkward.num(fatjets.particleNet_WvsQCD[(fatjets.particleNet_WvsQCD > 0.94)&(fatjets.Hqqqq_vsQCDTop<=0.4)]) == 1
        FH_1Wfatjet_cat = (~SL_boosted_cat) & (~SL_fullyresovled_cat) & (~SL_merged_boosted_cat) & (~FH_boosted)&(~FH_2Wfatjet_cat)&(n_leptons_iso==0) & (n_leptons_noiso == 0) & (n_fatjets >=1) & (n_jets >=2) & (selection_fatjet_WvsQCD_SB_1F) # 1 jet for FH
        # ----------------------------------------------------------------------------------------------------#
        # add resolved FH category with (>=4 AK4 jets )
        FH_fully_resovled_cat =(~SL_boosted_cat) & (~SL_fullyresovled_cat) & (~SL_merged_boosted_cat) & (~FH_boosted) & (~FH_2Wfatjet_cat)&(~FH_1Wfatjet_cat)& (n_leptons_iso==0) & (n_leptons_noiso == 0) & (n_jets>=4)# 4 jets for FH

      
        flatten_n_jets = awkward.num(jets.pt)
        category = awkward.zeros_like(flatten_n_jets)
        category = awkward.fill_none(category, 0)
        # add the priority for each category
        # 0: no category
        # 1: SL_boosted_cat
        # 2: SL_fullyresovled_cat
        # 3: SL_merged_boosted_cat
        # 4: FH_boosted
        # 5: FH_2Wfatjet_cat
        # 6: FH_1Wfatjet_cat
        # 7: FH_fully_resovled_cat

        # no lepton
        category = awkward.where(SL_boosted_cat, awkward.ones_like(category)*1, category)
        category = awkward.where(SL_fullyresovled_cat, awkward.ones_like(category)*2, category)
        category = awkward.where(SL_merged_boosted_cat, awkward.ones_like(category)*3, category)
        category = awkward.where(FH_boosted, awkward.ones_like(category)*4, category)
        # with isolated lepton
        category = awkward.where(FH_2Wfatjet_cat, awkward.ones_like(category)*5, category)
        category = awkward.where(FH_1Wfatjet_cat, awkward.ones_like(category)*6, category)
        category = awkward.where(FH_fully_resovled_cat, awkward.ones_like(category)*7, category)
        category_cut = (category >= 0) # cut the events with category == 0
        awkward_utils.add_field(events, "category", category) 

        presel_cut = (photon_id_cut) & (category_cut) & (Z_veto_cut)

        self.register_cuts(
            names=["Photon id Selection","category_cut", "Z_veto_cut"],
            results=[photon_id_cut, category_cut, Z_veto_cut]
        )
        return presel_cut, events