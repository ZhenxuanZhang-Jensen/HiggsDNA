# This script illustrates a full workflow example

import higgs_dna
import uproot
import json
import awkward
import copy
import pandas
import glob
import numpy
import os

from higgs_dna.analysis import AnalysisManager
from higgs_dna.taggers.dummy_taggers import THQTagger, TTHTagger
from higgs_dna.taggers.diphoton_tagger import DiphotonTagger
from higgs_dna.taggers.tag_sequence import TagSequence
from higgs_dna.systematics.systematics_producer import SystematicsProducer
from higgs_dna.systematics import photon_systematics, lepton_systematics, jet_systematics

from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger("DEBUG", "output/log.txt")

# 1. Define samples to run over
samples = {
        "sample_list" : ["ttH_M125"],
        "years" : ["2017", "2018"]
}

# Step 2: Define tag sequence (i.e. what selection you want to perform) 
# Step 2.1: Define individual taggers
tag_sequence = [
        {
            "module_name" : "higgs_dna.taggers.diphoton_tagger",
            "tagger" : "DiphotonTagger"
        },
        [
            { "module_name" : "higgs_dna.taggers.dummy_taggers", "tagger" : "TTHTagger" , "kwargs" : { "name" : "tth_dummy_tagger"}},
            { "module_name" : "higgs_dna.taggers.dummy_taggers", "tagger" : "THQTagger" , "kwargs" : { "name" : "thq_dummy_tagger"}}
        ]
]


"""
with open("metadata/example.json", "r") as f_in:
    options = json.load(f_in)

diphoton_tagger = DiphotonTagger(
        name = "diphoton_tagger",
        options = options["taggers"]["diphoton_tagger"]
)

tth_tagger = TTHTagger(
        name = "tth_dummy_tagger",
        options = {}
)

thq_tagger = THQTagger(
        name = "thq_dummy_tagger",
        options = {}
)

# Step 2.2: create TagSequence from Taggers 
tag_sequence = TagSequence(
        tag_list = [
            diphoton_tagger,
            [ # list within the list means these are enforced to be orthogonal
                tth_tagger, # first priority
                thq_tagger  # second priority
            ]
        ]
)
"""

# Step 3: define systematics (i.e. what corrections and systematic uncertainties you want to apply)
syst_options = {
        "weights" : {
            "dummy_theory_sf" : {
                "type" : "event", # event-level quantity
                "method" : "from_branch", # read the variations from branches in the input nanoAOD
                "branches" : { "central" : "genWeight", "up" : "genWeight", "down" : "genWeight" }, # these are the branches
                "modify_central_weight" : False # DO NOT multiply the central event weight by the central variation
            },
            "dummy_electron_sf" : {
                "type" : "object", # object-level quantity
                "method" : "from_function", # calculate this on-the-fly from a function
                "function" : { # use this function
                    "module_name" : "higgs_dna.systematics.lepton_systematics", # make sure this module is imported!
                    "name" : "dummy_lepton_sf" # name of actual function (i.e. higgs_dna.systematics.lepton_systematics.dummy_lepton_sf)
                },
                "modify_central_weight" : True, # DO multiply the central event weight by the central variation
                "input_collection" : "Electron", # Name of record to calculate the per-object weights on
                "target_collections" : [
                    "SelectedElectron" # Name of record to use when translating this to a per-event weight
                ]
            },
            "dummy_muon_sf" : {
                "type" : "object",
                "method" : "from_function",
                "function" : {
                    "module_name" : "higgs_dna.systematics.lepton_systematics", # make sure this module is imported!
                    "name" : "dummy_lepton_sf"
                },
                "modify_central_weight" : True,
                "input_collection" : "Muon",
                "target_collections" : [ # can have multiple target collections
                    "SelectedMuon",
                    "SelectedMuonTight"
                ],
                "modifies_taggers" :  {
                    "SelectedMuon" : ["tth_dummy_tagger"],
                    "SelectedMuonTight" : ["thq_dummy_tagger"]
                }
            },
            "dummy_ttH_theory_sf_v2" : {
                "type" : "event",
                "method" : "from_branch",
                "branches" : { "central" : "genWeight", "up" : "genWeight", "down" : "genWeight" },
                "modify_central_weight" : False,
                "samples" : ["ttH_M125"]
            },
            "diphoton_presel_sf" : {
                "type" : "object",
                "method" : "from_function",
                "function" : {
                    "module_name" : "higgs_dna.systematics.photon_systematics", # make sure this module is imported!
                    "name" : "photon_preselection_sf",
                },
                "modify_central_weight" : True,
                "input_collection" : "Photon",
                "target_collections" : [
                    ("Diphoton", "Photon")
                ]
            }
        },
        "independent_collections" : {
            "dummy_photon_pt_syst1" : {
                "method" : "from_branch",
                "branch_modified" : ("Photon", "pt"),
                "branches" : {
                    "up" : ("Photon", "pt"),
                    "down" : ("Photon", "pt")
                }
             },
            "dummy_photon_pt_syst2" : {
                "method" : "from_function",
                "branch_modified" : ("Photon", "pt"),
                "function" : {
                    "module_name" : "higgs_dna.systematics.photon_systematics", # make sure this module is imported!
                    "name" : "dummy_photon_pt_syst"
                }   
            },  
            "dummy_jes_syst" : {
                "method" : "from_function",
                "branch_modified" : ("Jet", "pt"),
                "function" : {
                    "module_name" : "higgs_dna.systematics.jet_systematics", # make sure this module is imported!
                    "name" : "dummy_jes_syst"
                }   
            }   
        }   
}       

# Add more dummy systematics to get closer to realistic number of systs
# Choose 30 additional weight systematics and 30 additional independent collections, based on t->Hq (H->gg) FCNC analysis
n_systs = 1
for i in range(n_systs):
    syst_options["weights"]["diphoton_presel_sf_%d" % i] = syst_options["weights"]["diphoton_presel_sf"]
    syst_options["independent_collections"]["dummy_photon_pt_syst2_%d" % i] = syst_options["independent_collections"]["dummy_photon_pt_syst2"]

# Step 4: Define content in output file
variables_of_interest = [("Diphoton", "mass"), ("Diphoton", "pt")]

# Step 5: Define job submission options
batch_system = "local"
n_files_per_job = 2

# Wrap all of this in an analysis manager class
my_analysis = AnalysisManager(
        name = "my_analysis",
        samples = samples,
        tag_sequence = tag_sequence,
        systematics = syst_options,
        variables_of_interest = variables_of_interest,
        batch_system = batch_system,
        n_files_per_job = n_files_per_job
)

my_analysis.run()
