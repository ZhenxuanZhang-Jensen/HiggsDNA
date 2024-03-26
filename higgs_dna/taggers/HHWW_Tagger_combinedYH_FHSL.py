import logging

import awkward
import vector
import numpy
from higgs_dna.selections import (fatjet_selections, jet_selections,
                                  lepton_selections,gen_selections)
from higgs_dna.taggers.tagger import NOMINAL_TAG, Tagger
from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.lepton_systematics import highptmuonsf
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
def delta_phi(obj1, obj2):
    """
    Select objects from objects1 which are at least min_dphi away from all objects in objects2.
    :param objects1: objects which are required to be at least min_dphi away from all objects in objects2 
    :type objects1: awkward.highlevel.Array
    :param objects2: objects which are all objects in objects1 must be at leats min_dphi away from
    :type objects2: awkward.highlevel.Array
    :param min_dphi: minimum delta R between objects
    :type min_dphi: float
    :return: boolean array of objects in objects1 which pass delta_R requirement
    :rtype: awkward.highlevel.Array
    """
    #obj1:fatjets obj2:puppiMET     
    fatjet_1_phi=awkward.fill_none(awkward.pad_none(obj1.phi,6),-999,axis=1)[:,0]
    fatjet_2_phi=awkward.fill_none(awkward.pad_none(obj1.phi,6),-999,axis=1)[:,1]
    fatjet_3_phi=awkward.fill_none(awkward.pad_none(obj1.phi,6),-999,axis=1)[:,2]
    fatjet_4_phi=awkward.fill_none(awkward.pad_none(obj1.phi,6),-999,axis=1)[:,3]
    fatjet_5_phi=awkward.fill_none(awkward.pad_none(obj1.phi,6),-999,axis=1)[:,4]
    fatjet_6_phi=awkward.fill_none(awkward.pad_none(obj1.phi,6),-999,axis=1)[:,5]

    event_dphi=awkward.concatenate([awkward.unflatten((fatjet_1_phi-awkward.flatten(obj2.phi)),counts=1),awkward.unflatten((fatjet_2_phi-awkward.flatten(obj2.phi)),counts=1),awkward.unflatten((fatjet_3_phi-awkward.flatten(obj2.phi)),counts=1),awkward.unflatten((fatjet_4_phi-awkward.flatten(obj2.phi)),counts=1),awkward.unflatten((fatjet_5_phi-awkward.flatten(obj2.phi)),counts=1),awkward.unflatten((fatjet_6_phi-awkward.flatten(obj2.phi)),counts=1)],axis=1)
    event_negative_fatjetphi=awkward.fill_none((event_dphi.mask[event_dphi<-numpy.pi]+2*numpy.pi),numpy.nan)
    event_positive_fatjetphi=awkward.fill_none((-event_dphi.mask[event_dphi>numpy.pi]+2*numpy.pi),numpy.nan)
    event_middphi=awkward.fill_none((event_dphi.mask[(event_dphi<=numpy.pi)&(event_dphi>=-numpy.pi)]),numpy.nan)

    fatjet_new_phi = numpy.where(~numpy.isnan(event_positive_fatjetphi), event_positive_fatjetphi, event_negative_fatjetphi)

    dphi=numpy.where(~numpy.isnan(fatjet_new_phi), fatjet_new_phi, event_middphi)
    # dphi=dphi[abs(dphi)<2*numpy.pi]
    return dphi
logger = logging.getLogger(__name__)


DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {
    "electrons_noiso": {
        "pt": 10,
        "dr_photons": 0.4,
        "id": "WP80iso_WP90noniso"
    },
    "electron_iso": {
        "pt": 10,
        "dr_photons": 0.4,
        "id": "WP80",
    },
    "electrons_all": {
        "pt": 10,
        "dr_photons": 0.4
    },
    "muons_all": {
        "pt": 15,
        "dr_photons": 0.4
    },
    "muons_noiso": {
        "Tunept" : 15.0, 
        "dr_photons": 0.4,
        "id" : "highptId",
        "global":True,
        "eta":2.4
         },
    "muons_iso": {
        "pt": 15,
        "dr_photons": 0.4,
        "global":True,
        "id":"tight",
        "eta":2.4,
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
        "dr_fatjets":0.8
    },
    "fatjets": {
        "pt": 100.0,
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
    # "photon_id": -0.9,
    "photon_id": -0.9,
    "btag_wp": {
        "2016": 0.3093,
        "2017": 0.3040,
        "2018": 0.2783
    },
    "gen_info" : {
        "is_Signal" : False, #attention: change in HHWW_preselection.json
    }  
}


class HHWW_Tagger_combinedYH_FHSL(Tagger):
    """
    HHWW Preselection tagger for tutorial
    """

    def __init__(self, name, options={}, is_data=None, year=None,output_dir=None):
        super(HHWW_Tagger_combinedYH_FHSL, self).__init__(name, options, is_data, year,output_dir)

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
        if not self.is_data and self.options["gen_info"]["is_Signal"]:    
            gen_obj_1, gen_obj_2, gen_obj_3, gen_obj_4, gen_lead_W_fromH, gen_sublead_W_fromH, gen_H_to_gg, gen_H_to_WW, gen_lead_g_fromH, gen_sublead_g_fromH = gen_selections.select_ww_to_qqlv_or_qqqq(events)
        # add object fields for gen objects
            gen_obj_1 = awkward_utils.add_field(
                events = events,
                name="gen_obj_1",
                data = awkward.unflatten(gen_obj_1,counts=1)
            )
            awkward_utils.add_object_fields(
                events = events,
                name = "gen_obj_1",
                objects = gen_obj_1,
                n_objects = 1,
                dummy_value = -999
            )
            gen_obj_2 = awkward_utils.add_field(
                events = events,
                name="gen_obj_2",
                data = awkward.unflatten(gen_obj_2,counts=1)
            )
            awkward_utils.add_object_fields(
                events = events,
                name = "gen_obj_2",
                objects = gen_obj_2,
                n_objects = 1,          
                dummy_value = -999                      
            )
            gen_obj_3 = awkward_utils.add_field(
                events = events,
                name="gen_obj_3",
                data = awkward.unflatten(gen_obj_3,counts=1)
            )                           
            awkward_utils.add_object_fields(                        
                events = events,                    
                name = "gen_obj_3",             
                objects = gen_obj_3,                
                n_objects = 1,                      
                dummy_value = -999              
            )                   
            gen_obj_4 = awkward_utils.add_field(
                events = events,
                name="gen_obj_4",
                data = awkward.unflatten(gen_obj_4,counts=1)
            )                           
            awkward_utils.add_object_fields(                        
                events = events,                    
                name = "gen_obj_4",             
                objects = gen_obj_4,                    
                n_objects = 1,                          
                dummy_value = -999
            )                  
            gen_lead_W_fromH = awkward_utils.add_field(
                events = events,
                name="gen_lead_W_fromH",
                data = awkward.unflatten(gen_lead_W_fromH,counts=1)
            )                           
            awkward_utils.add_object_fields(                        
                events = events,                    
                name = "gen_lead_W_fromH",             
                objects = gen_lead_W_fromH,                    
                n_objects = 1,                          
                dummy_value = -999
            )                  
            gen_sublead_W_fromH = awkward_utils.add_field(
                events = events,
                name="gen_sublead_W_fromH",
                data = awkward.unflatten(gen_sublead_W_fromH,counts=1)
            )                           
            awkward_utils.add_object_fields(                        
                events = events,                    
                name = "gen_sublead_W_fromH",             
                objects = gen_sublead_W_fromH,                    
                n_objects = 1,                          
                dummy_value = -999
            )                  
            gen_H_to_gg = awkward_utils.add_field(
                events = events,
                name="gen_H_to_gg",
                data = awkward.unflatten(gen_H_to_gg,counts=1)
            )                                                                                   
            awkward_utils.add_object_fields(                            
                events = events,                    
                name = "gen_H_to_gg",                   
                objects = gen_H_to_gg,          
                n_objects = 1,
                dummy_value = -999                          
            )    
            gen_H_to_gg = awkward_utils.add_field(
                events = events,
                name="gen_H_to_gg",
                data = gen_H_to_gg
            )                          
            awkward_utils.add_object_fields(                        
                events = events,                            
                name = "gen_H_to_gg",                       
                objects = gen_H_to_gg,                  
                n_objects = 1,
                dummy_value = -999
            )       
            gen_lead_g_fromH = awkward_utils.add_field(
                events = events,                        
                name="gen_lead_g_fromH",                    
                data = awkward.unflatten(gen_lead_g_fromH,counts=1)             
            )                       
            awkward_utils.add_object_fields(    
                events = events,
                name = "gen_lead_g_fromH",                  
                objects = gen_lead_g_fromH,                     
                n_objects = 1,                      
                dummy_value = -999              
            )                       
            gen_sublead_g_fromH = awkward_utils.add_field(
                events = events,                
                name="gen_sublead_g_fromH",                         
                data = awkward.unflatten(gen_sublead_g_fromH,counts=1)                  
            )                   
            awkward_utils.add_object_fields(            
                events = events,            
                name = "gen_sublead_g_fromH",
                objects = gen_sublead_g_fromH,
                n_objects = 1,                  
                dummy_value = -999  
            )
            gen_H_to_WW = awkward_utils.add_field(
                events = events,                
                name="gen_H_to_WW",                         
                data = gen_H_to_WW              
            )                   
            awkward_utils.add_object_fields(            
                events = events,            
                name = "gen_H_to_WW",
                objects = gen_H_to_WW,
                n_objects = 1,                  
                dummy_value = -999  
        )
            events= gen_selections.get_minmaxID(events)

        # add object fields for gen objects
        if not self.is_data and not self.options["gen_info"]["is_Signal"]:  
            events= gen_selections.get_minmaxID(events)

   
        logger.debug("event fields: %s" %events.fields)
        # logger.debug('After Diphoton selection')
        puppiMET = awkward.zip({'phi':events.PuppiMET_phi,'pt':events.PuppiMET_pt,'sumEt':events.PuppiMET_sumEt})
        puppiMET = awkward_utils.add_field(
            events=events,
            name="puppiMET",
            data=awkward.unflatten(puppiMET,counts=1)
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "PuppiMET",
            objects = puppiMET,
            n_objects = 1,
            dummy_value = -999
        )

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
        
        # Iso-Electrons
        electron_iso_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electron_iso"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons_noiso"]["dr_photons"]
                }
            },
            name="SelectedElectron_iso",
            tagger=self
        )

        electron_iso = awkward_utils.add_field(
            events=events,
            name="SelectedElectron_iso",
            data=events.Electron[electron_iso_cut]
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "electron_iso",
            objects = electron_iso,
            n_objects = 1,
            dummy_value = -999
        )
        
        # all electrons with pt > 200 for SL boosted 
        electrons_all_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons_all"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons_all"]["dr_photons"]
                }
            },
            name="SelectedElectrons_all",
            tagger=self
        )
        electrons_all = awkward_utils.add_field(
            events=events,
            name="SelectedElectrons_all",
            data=events.Electron[electrons_all_cut]
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "electrons_all",
            objects = electrons_all,
            n_objects = 3,
            dummy_value = -999
        )
        
        # Noiso-Muons
        muon_noiso_cut = lepton_selections.select_nonisomuons(
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
        muons_noiso['Tunept']= muons_noiso.pt * muons_noiso.tunepRelPt


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
        
        muons_all_cut = lepton_selections.select_muons(
            muons=events.Muon,
            options=self.options["muons_all"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["muons_all"]["dr_photons"]
                }
            },
            name="SelectedMuon_all",
            tagger=self
        )
        muons_all = awkward_utils.add_field(
            events=events,
            name="SelectedMuon_all",
            data=events.Muon[muons_all_cut]
        )
        awkward_utils.add_object_fields(
            events = events,
            name = "muons_all",
            objects = muons_all,
            n_objects = 3,
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



        # --------------------- if fatjet branches are not empty --------------------- #
        # Fat jets

        # add the H jet tagger for SL&FH channel((H3q+H4q+Hlvqq)/(H3q+H4q+Hlvqq+QCD+Top))
        # ----------------------------------- PN 3q4q ---------------------------------- #
        PN_sigs_4q = events.FatJet.inclParTMDV1_probHWqqWq0c+events.FatJet.inclParTMDV1_probHWqqWq1c+events.FatJet.inclParTMDV1_probHWqqWq2c+events.FatJet.inclParTMDV1_probHWqqWqq0c+events.FatJet.inclParTMDV1_probHWqqWqq1c+events.FatJet.inclParTMDV1_probHWqqWqq2c
        PN_bkgs_4q = events.FatJet.inclParTMDV1_probQCDb+events.FatJet.inclParTMDV1_probQCDbb+events.FatJet.inclParTMDV1_probQCDc+events.FatJet.inclParTMDV1_probQCDcc+events.FatJet.inclParTMDV1_probQCDothers+events.FatJet.inclParTMDV1_probTopbWq0c+events.FatJet.inclParTMDV1_probTopbWq1c+events.FatJet.inclParTMDV1_probTopbWqq0c+events.FatJet.inclParTMDV1_probTopbWqq1c
        
        fatjet_tmp = events.FatJet
        # fatjet_tmp['dphi_puppiMET']=(awkward.unflatten(events.PuppiMET_phi,counts=1)-fatjet_tmp.phi)
        fatjet_tmp['Hqqqq_vsQCDTop'] = PN_sigs_4q / (PN_bkgs_4q + PN_sigs_4q)  
        fatjet_tmp['WvsQCDMD']=(events.FatJet.particleNetMD_Xcc + events.FatJet.particleNetMD_Xqq)/(events.FatJet.particleNetMD_Xcc + events.FatJet.particleNetMD_Xqq + events.FatJet.particleNetMD_QCD)
        events.FatJet = fatjet_tmp
        print("events.FatJet exsiting field:\n",events.FatJet.fields)
        # ----------------------------------- Subjet ---------------------------------- #
        subjet = awkward_utils.add_field(
            events = events,
            name = "SelectedSubJet",
            data = events.SubJet
        )   
        # ------------------------------------- - ------------------------------------ #
        fatjet_cut = fatjet_selections.select_fatjets(
            fatjets = events.FatJet,
            subjets = events.SubJet,
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
        print("puppiMET", puppiMET)
        print("delta phi puppiMET", delta_phi(fatjets,puppiMET))
        try:
            fatjets['dphi_MET']=delta_phi(fatjets,puppiMET)[delta_phi(fatjets,puppiMET)>-10]
        except:
            fatjets['dphi_MET']= awkward.ones_like(fatjets.pt)*-999 # no fatjet passing the selection
        lepton_noniso=awkward.concatenate([electrons_noiso,muons_noiso],axis=1)

        awkward_utils.add_object_fields(
        events=events,
        name="fatjet",
        objects=fatjets[awkward.argsort(fatjets.pt, ascending=False, axis=-1)],
        n_objects=3,
        dummy_value=-999
        )
        fatjets=fatjets[awkward.argsort(fatjets.pt, ascending=False, axis=-1)]

       
        
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
                },
                "fatjets": {
                    "objects": events.SelectedFatJet,
                    "min_dr": self.options["jets"]["dr_fatjets"]
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
        awkward_utils.add_object_fields(
            events=events,
            name="jet",
            objects=jets,
            n_objects=7,
            dummy_value=-999
        )
        electrons_noiso = awkward.Array(electrons_noiso, with_name="Momentum4D")
        electron_iso = awkward.Array(electron_iso, with_name="Momentum4D")
        muons_noiso = awkward.Array(muons_noiso, with_name="Momentum4D")
        muons_iso = awkward.Array(muons_iso, with_name="Momentum4D")
        # bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis=1, ascending=False)]
        # bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]]

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        # Preselection
        
        n_electrons_noiso = awkward.num(electrons_noiso)
        n_electron_iso = awkward.num(electron_iso)
        n_muons_noiso = awkward.num(muons_noiso)
        n_muons_iso = awkward.num(muons_iso)
        n_leptons_noiso = n_electrons_noiso + n_muons_noiso 
        n_leptons_iso = n_electron_iso + n_muons_iso
        n_leptons_all_electrons = awkward.num(electrons_all) + n_muons_iso + n_muons_noiso 
        awkward_utils.add_field(events,"nGoodnonisoelectrons",n_electrons_noiso)
        # n_leptons_all_electrons
        awkward_utils.add_field(events,"nGoodallleptons",n_leptons_all_electrons)
        awkward_utils.add_field(events,"nGoodisoelectrons",n_electron_iso)
        awkward_utils.add_field(events,"nGoodnonisomuons",n_muons_noiso)
        awkward_utils.add_field(events,"nGoodisomuons",n_muons_iso)
        awkward_utils.add_field(events,"nGoodisoleptons",n_leptons_iso)
        awkward_utils.add_field(events,"nGoodnonisoleptons",n_leptons_noiso)
        n_jets = awkward.num(jets)
        awkward_utils.add_field(events,"nGoodAK4jets",n_jets)
        n_fatjets = awkward.num(fatjets)

        awkward_utils.add_field(events,"nGoodAK8jets",n_fatjets)

        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (
            events.SubleadPhoton.mvaID > self.options["photon_id"])
        # ----------------------------------------------------------------------------------------------------#
  

        # ----------------------------------------------------------------------------------------------------#
        # new YH category for PNN training
        # first category: 1 lepton(iso or noniso) + 1 Wfatjet
        selection_fatjet_WvsQCD = awkward.num(fatjets.WvsQCDMD[(fatjets.WvsQCDMD > 0.58)]) >= 1
        boosted_YH_SL_cat = (((n_leptons_iso >= 1) | (n_leptons_noiso >= 1)) & (n_fatjets >=1)) # boosted 1 jet for SL channel with isolated lep
        # boosted_YH_SL_cat_v2 = ((n_leptons_all_electrons==1) & (n_fatjets >=1) & (selection_fatjet_WvsQCD)) # boosted 1 jet for SL channel with isolated lep
        # second category: 0 lepton + 1 Wfatjet or 1 Higgs fatjet
        selection_fatjet_HvsQCD = awkward.num(fatjets.Hqqqq_vsQCDTop[(fatjets.Hqqqq_vsQCDTop > 0.2)]) >= 1
        boosted_YH_FH_cat = (n_leptons_iso == 0) & (n_leptons_noiso == 0) & (n_fatjets >=1) & ( (selection_fatjet_WvsQCD)|(selection_fatjet_HvsQCD)) # boosted 1 jet for SL channel with isolated lep
        # resolved case
        FH_fully_resovled_cat = (n_leptons_iso==0) & (n_leptons_noiso == 0) & (n_jets>=2)
        SL_fullyresovled_cat = (n_leptons_iso >= 1) & (n_jets >=1) # resolved 2 jets for SL channel with isolated lep
      
       
       
        # ----------------------------------------------------------------------------------------------------#
      
      
        flatten_n_jets = awkward.num(jets.pt)
        category = awkward.zeros_like(flatten_n_jets)
        category = awkward.fill_none(category, 0)
        
        # new YH category for PNN training
        # third category: SL_fullresovled_cat
        category = awkward.where(SL_fullyresovled_cat, awkward.ones_like(category)*3, category)
        # fourth category: FH_fully_resovled_cat
        category = awkward.where(FH_fully_resovled_cat, awkward.ones_like(category)*4, category)
        # first category: 1 lepton(iso or noniso) + 1 Wfatjet
        # category = awkward.where(boosted_YH_SL_cat_v2, awkward.ones_like(category)*1.5, category)
        category = awkward.where(boosted_YH_SL_cat, awkward.ones_like(category)*1, category)
        # second category: 0 lepton + 1 Wfatjet or 1 Higgs fatjet
        category = awkward.where(boosted_YH_FH_cat, awkward.ones_like(category)*2, category)
        # ----------------------------------------------------------------------------------------------------#
        category_cut = (category > 0) # cut the events with category == 0
        awkward_utils.add_field(events, "category", category) 
        # if not self.is_data and self.options["gen_info"]["is_Signal"]: 
        #     events=highptmuonsf(events)
        WP80=awkward.concatenate([awkward.unflatten(events.LeadPhoton.mvaID_WP80,counts=1),awkward.unflatten(events.SubleadPhoton.mvaID_WP80,counts=1)],axis=-1)
        WP90=awkward.concatenate([awkward.unflatten(events.LeadPhoton.mvaID_WP90,counts=1),awkward.unflatten(events.SubleadPhoton.mvaID_WP90,counts=1)],axis=-1)
        PhotonID=awkward.zip({"WP80":WP80,"WP90":WP90})
        events['is_passPhotonMVA80']=awkward.num(PhotonID[PhotonID.WP80==True])==2
        events['is_passPhotonMVA90']=awkward.num(PhotonID[PhotonID.WP90==True])==2
        
        presel_cut = (photon_id_cut) & (category_cut) & (Z_veto_cut)
        
        self.register_cuts(
                names=["Photon id Selection","category_cut", "Z_veto_cut"],
                results=[photon_id_cut, category_cut, Z_veto_cut])

        return presel_cut, events