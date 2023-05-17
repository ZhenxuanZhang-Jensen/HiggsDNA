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

    def __init__(self, name, options={}, is_data=None, year=None):
        super(HHWW_Preselection_FHSL, self).__init__(name, options, is_data, year)

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
        
        
       
        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (
            events.SubleadPhoton.mvaID > self.options["photon_id"])

        presel_cut = (photon_id_cut)

        self.register_cuts(
            names=["Photon id Selection","category_cut"],
            results=[photon_id_cut]
        )
        return presel_cut, events