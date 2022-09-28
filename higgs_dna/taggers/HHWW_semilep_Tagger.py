import logging

import awkward
import vector
import numpy
from higgs_dna.selections import (fatjet_selections, jet_selections,
                                  lepton_selections,gen_selections,object_selections)

from higgs_dna.taggers.tagger import NOMINAL_TAG, Tagger
from higgs_dna.taggers import diphoton_tagger
from higgs_dna.utils import awkward_utils, misc_utils

vector.register_awkward()

logger = logging.getLogger(__name__)

DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {
    "electrons": {
        "pt": 10.0,
        "dr_photons": 0.2
    },
    "muons": {
        "pt": 10.0,
        "eta":2.7,
        "dr_photons": 0.2
    },
    "jets": {
        "pt": 25.0, # attention this is the one exact same as old framework, make this 20 GeV(loose) for further analysis, we all know the higgs-like ak4 jets pt can be very small
        "eta": 2.5,
        "dr_photons": 0.4,
        "dr_electrons": 0.4,
        "dr_muons": 0.4,
        "dr_jets": 0.4,
    },
    "fatjets": {
        "pt": 200.0,
        "eta": 2.4,
        "dr_photons": 0.8,
        "dr_electrons": 0.8,
        "dr_muons": 0.8,
        "dr_jets":0.8
    },
    "photon_id": -0.9,
    "btag_wp": {
        "2016": 0.3093,
        "2017": 0.3040,
        "2018": 0.2783
    },
    "gen_info" : {
        "is_Signal" : True, #attention: change in HHWW_preselection.json
    }  
}


class HHWW_Preselection(Tagger):
    """
    HHWW Preselection tagger for tutorial
    """

    def __init__(self, name, options={}, is_data=None, year=None):
        super(HHWW_Preselection, self).__init__(name, options, is_data, year)

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
        logger.debug("Is Signal: %s" %self.options["gen_info"]["is_Signal"])
        if not self.is_data and self.options["gen_info"]["is_Signal"]:    
            gen_selections.gen_Hww_2q2l(events)
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
       # iso_cut = electrons["Electron_mvaFall17V2Iso_WPL"][electrons["Electron_mvaFall17V2Iso_WPL"]==True]
        #electrons=electrons[iso_cut]
        ######################################

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
        # --------------------- if fatjet branches are not empty --------------------- #
        # if len(events.FatJet.pt > 0 ):
        # Fat jets
        fatjet_cut = fatjet_selections.select_fatjets(
            fatjets = events.FatJet,
            options = self.options["fatjets"],
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
                    }
                },
            name = "SelectedFatJet",
            tagger = self
        )
        fatjets = awkward_utils.add_field(
            events = events,
            name = "SelectedFatJet",
            data = events.FatJet[fatjet_cut]
        )   

        fatjet_H_cut = (fatjets.deepTagMD_HbbvsQCD<0.6) & (fatjets.deepTagMD_H4qvsQCD>0.4) & (fatjets.pt>300)

        fatjets_H = awkward_utils.add_field(
            events = events,
            name = "SelectedFatJet_H",
            data = fatjets[fatjet_H_cut]
        )   


        awkward_utils.add_object_fields(
        events=events,
        name="fatjet_H",
        objects=fatjets[fatjet_H_cut][awkward.argsort(fatjets[fatjet_H_cut].deepTagMD_H4qvsQCD, ascending=False, axis=1)],
        n_objects=1,
        dummy_value=-999
        ) # apply the inverse bb cuts

        fatjet_W_cut = (fatjets.deepTagMD_HbbvsQCD<0.6) & (fatjets.deepTagMD_WvsQCD>0.4) & (fatjets.pt>200)

        fatjets_W = awkward_utils.add_field(
            events = events,
            name = "SelectedFatJet_W",
            data = fatjets[fatjet_W_cut]
        )   
        awkward_utils.add_object_fields(
        events=events,
        name="fatjet_W",
        objects=fatjets[fatjet_W_cut][awkward.argsort(fatjets[fatjet_W_cut].deepTagMD_WvsQCD, ascending=False, axis=1)],
        n_objects=1,
        dummy_value=-999
        ) # apply the inverse bb cuts

        # fatjets["Tau4_Tau2"]= (fatjets["tau4"]/fatjets["tau2"])
        # fatjets["Tau2_Tau1"]= (fatjets["tau2"]/fatjets["tau1"])
        # fatjets = fatjets[awkward.argsort(fatjets.deepTagMD_HbbvsQCD, ascending=True, axis=1)]

        
        # Jets
        jet_cut = jet_selections.select_jets(
            jets=events.Jet,
            options=self.options["jets"],
 #           clean = {},
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
                   }, #close for QCD samples
            name = "SelectedJet",
            tagger=self
        )
        jets = awkward_utils.add_field(
            events=events,
            name="SelectedJet",
            data=events.Jet[jet_cut]
        )
    

        # ---------------------------------------------------------------------------- #
        #              TODO try to add order with WinvM jets in output.parquet              #
        # ---------------------------------------------------------------------------- #
      
        gen_l1_p4,gen_q1_p4,gen_q2_p4=gen_selections.gen_Hww_2q2l(events)
        jet_p4 = vector.awk(
            {
                "pt" : jets["pt"],
                "eta" : jets["eta"],
                "phi" : jets["phi"],
                "mass" : jets["mass"]
            },
            with_name = "Momentum4D"
        )
        lead_photon_vec = vector.awk(
            {
                "pt" : events["LeadPhoton", "pt"],
                "eta" : events["LeadPhoton", "eta"],
                "phi" : events["LeadPhoton", "phi"],
                "mass" : events["LeadPhoton", "mass"] 
            },
            with_name = "Momentum4D"
        )
        sublead_photon_vec = vector.awk(
            {
                "pt" : events["SubleadPhoton", "pt"],
                "eta": events["SubleadPhoton", "eta"],
                "phi": events["SubleadPhoton", "phi"],
                "mass": events["SubleadPhoton", "mass"]
            },
            with_name = "Momentum4D"
        )

        jets["deltaR_q1"] = jet_p4.deltaR(gen_q1_p4)
        jets["deltaR_q2"] = jet_p4.deltaR(gen_q2_p4)

        jets["deltaR_pho1"] = jet_p4.deltaR(lead_photon_vec)
        jets["deltaR_pho2"] = jet_p4.deltaR(sublead_photon_vec)

        awkward_utils.add_object_fields(
            events=events,
            name="jet",
            objects=jets,
            n_objects=7,
            dummy_value=-999
        )
        # bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis=1, ascending=False)]
        # bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]]
 #       awkward_utils.add_fields(
  #          events=events,
   #         name="sigjet",
    #        data = selectedjet
       # )
        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        electrons = awkward.Array(electrons, with_name="Momentum4D")
        muons = awkward.Array(muons, with_name="Momentum4D")

        # Preselection
        n_electrons = awkward.num(electrons)
        n_muons = awkward.num(muons)
        n_photons = awkward.num(events.Photon)
        n_leptons = n_electrons + n_muons
        # n_diphotons = awkward.num(events.Diphoton)
        # logger.debug(" the N_diphoton : %f" % (n_diphotons))
        n_jets = awkward.num(jets)
        awkward_utils.add_field(events,"nGoodAK4jets",n_jets)
        n_fatjets = awkward.num(fatjets)
        n_fatjets_W = awkward.num(fatjets_W)
        n_fatjets_H = awkward.num(fatjets_H)
        awkward_utils.add_field(events,"nGoodAK4jets",n_jets)
        # n_bjets = awkward.num(bjets)

        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (
            events.SubleadPhoton.mvaID > self.options["photon_id"])
        # diphoton_pt_cut = (events.LeadPhoton.pt > self.options["photon_id"]) & (
        #     events.SubleadPhoton.pt > self.options["photon_id"])
        # photon_pixelseed_cut = (events.Photon.pixelSeed==0)

        # Hadronic presel
        # use priority to mark different category
        category_p3 = (n_jets>=2)
        category_p2 = (n_fatjets_W>=1) & (n_jets>=0)
        #category_p1 = (n_fatjets_H>=1)
        flatten_n_jets = awkward.num(jets.pt)
        category = awkward.zeros_like(flatten_n_jets)
        category = awkward.fill_none(category, 0)
        category = awkward.where(category_p3, awkward.ones_like(category)*3, category)
        category = awkward.where(category_p2, awkward.ones_like(category)*2, category)
       # category = awkward.where(category_p1, awkward.ones_like(category)*1, category)
        awkward_utils.add_field(events, "category", category) 
        category_cut = category >= 0 # attention category equal to 0 mean don't pass any selection 
        # OnlyFourJet_category = (n_leptons == 0) & (n_jets >= 4)
        Lepton_Selection = (n_leptons==1)&(n_jets >=2)

        # Semi-Leptonic presel
        # Semileptonic = (n_leptons == 1) & (n_jets >= 2)

        # Fully-Leptonic presel
        # FulllyLeptonic = (n_leptons >= 2)

        # presel_cut = (hadronic | Semileptonic | FulllyLeptonic)  & photon_id_cut
        # presel_FourJet_category = 
        presel_cut = (photon_id_cut) & (n_leptons==1) & (category_cut)

        self.register_cuts(
            names=["Photon Selection","Lepton Selection"],
            results=[photon_id_cut,Lepton_Selection]
        )


        return presel_cut, events 

