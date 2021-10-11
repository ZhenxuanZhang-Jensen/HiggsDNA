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
        "sample_list" : ["Data"],
        "years" : ["2016"]
}

# Step 2: Define tag sequence (i.e. what selection you want to perform) 
# Step 2.1: Define individual taggers

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
            diphoton_tagger
        ]
)


# Step 3: define systematics (i.e. what corrections and systematic uncertainties you want to apply)
syst_options = {
        "weights" : {
        },
        "independent_collections" : {
        }   
}       

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

my_analysis.run(max_jobs = 8)
