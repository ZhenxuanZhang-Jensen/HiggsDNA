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
        "id": "WPL_noniso",
        "non_iso": 0.4 #This is for mini isolation
        # "non_iso": None
    },
    "electrons_iso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id": "WPL",
        "iso": 0.4 #This is for mini isolation
        # "iso": None
    },
    "muons": {
        "pt" : 10.0, 
        "dr_photons": 0.4,
        "id" : "tight",
        "non_pfRelIso04_all":0.15,
        "non_iso":0.4, #This is for mini isolation
        # "non_iso": None,
        "global" : True
        },
    "muons_iso": {
        "pt": 10.0,
        "dr_photons": 0.4,
        "id":"tight",
        "dr_photons": 0.4,
        "pfRelIso04_all":0.15,
        "iso" : 0.4, # This is for mini isolation
        # "iso" : None,
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
        # "Hqqqq_qqlv_vsQCDTop": -999,
        # "Hqqqq_vsQCDTop": -999,
        "dr_photons": 0.8,
        "dr_electrons": 0.8,
        "dr_muons": 0.8
    },
    "fatjets_H": {
        "pt": 300,#300.0 fixed 0 
        "eta": 2.4,
        # "Hqqqq_qqlv_vsQCDTop" :0.2,
        # "Hqqqq_vsQCDTop" :0.8,
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
            gen_l1_p4, gen_q1_p4,gen_q2_p4 = gen_selections.gen_Hww_2q2l(events)        
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

       
        fatjet_tmp = events.FatJet
        
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
                # "electrons" : {
                #     "objects" : events.SelectedElectron_iso,
                #     "min_dr" : self.options["fatjets"]["dr_electrons"]
                # },
                # "muons" : {
                #     "objects" : events.SelectedMuon_iso,
                #     "min_dr" : self.options["fatjets"]["dr_muons"]
                #     },
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
                # "electrons" : {
                #     "objects" : events.SelectedElectron_iso,
                #     "min_dr" : self.options["fatjets_H"]["dr_electrons"]
                # },
                # "muons" : {
                #     "objects" : events.SelectedMuon_iso,
                #     "min_dr" : self.options["fatjets_H"]["dr_muons"]
                #     },
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
        objects=fatjets_H[awkward.argsort(fatjets_H.pt, ascending=False, axis=-1)],
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
        electrons = awkward.Array(electrons, with_name="Momentum4D")
        electrons_iso = awkward.Array(electrons_iso, with_name="Momentum4D")
        muons = awkward.Array(muons, with_name="Momentum4D")
        muons_iso = awkward.Array(muons_iso, with_name="Momentum4D")
        
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


        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (events.SubleadPhoton.mvaID > self.options["photon_id"])
        category_p4 = (n_leptons == 1)
        category_p3 = (n_leptons == 1)
        category_p2 = (n_leptons == 1)
        category_p1 = (n_leptons == 1)
        flatten_n_jets = awkward.num(jets.pt)
        category = awkward.zeros_like(flatten_n_jets)
        category = awkward.fill_none(category, 0)
        category = awkward.where(category_p4, awkward.ones_like(category)*4, category)
        category = awkward.where(category_p3, awkward.ones_like(category)*3, category)
        category = awkward.where(category_p2, awkward.ones_like(category)*2, category)
        category = awkward.where(category_p1, awkward.ones_like(category)*1, category)
        awkward_utils.add_field(events, "category", category) 
        category_cut = category >= 0 # attention category equal to 0 mean don't pass any selection 

        category_cut = (category >= 0) # cut the events with category == 0
        awkward_utils.add_field(events, "category", category) 

        presel_cut = (photon_id_cut) & (category_cut) & (Z_veto_cut)

        self.register_cuts(
            names=["Photon id Selection","category_cut", "Z_veto_cut"],
            results=[photon_id_cut, category_cut, Z_veto_cut]
        )
        return presel_cut, events