### This script compares yields of the diphoton preselection using added branches to custom nanoAOD vs central nanoAOD on data and MC 

import higgs_dna
import uproot
import json
import awkward
import copy
import pandas
import glob
import numpy
import os

from higgs_dna.taggers.dummy_taggers import THQTagger, TTHTagger
from higgs_dna.taggers.diphoton_tagger import DiphotonTagger
from higgs_dna.taggers.tag_sequence import TagSequence
from higgs_dna.systematics.systematics_producer import SystematicsProducer
from higgs_dna.systematics import photon_systematics, lepton_systematics, jet_systematics

from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger("DEBUG", "output/log.txt")

# Currently, we need to load the options json before loading any files,
# as the options json tells us which branches we need to read.
# This could be improved in the future: maybe HiggsDNA automatically detects which branches should be read? 
with open("metadata/example.json", "r") as f_in:
    options = json.load(f_in)

###############################
### Step 1: pick your files ###
###############################

# Quick and easy way of doing it for now, intended to be replaced with a SampleManager class
#file_mc = "root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201021/test_nanoaod_IDX.root"
file_mc = "root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201021/test_nanoaod_IDX.root"
file_data = "root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau/EGamma_Run2018D-22Jan2019-v2_MINIAOD_v0.6_20201021/test_nanoaod_IDX.root"
idx_mc = numpy.arange(1,2)
idx_data = numpy.arange(1,2)

events_base_data = []
events_base_mc = []
for i in idx_data:
    file_idx_data = file_data.replace("IDX", str(i))
    print(file_idx_data)
    with uproot.open(file_idx_data) as f:
        tree = f["Events"]
        events_base_data.append(tree.arrays([b for b in options["branches"] if "gen" not in b], library = "ak", how = "zip"))

for i in idx_mc:
    file_idx_mc = file_mc.replace("IDX", str(i))
    print(file_idx_mc)
    with uproot.open(file_idx_mc) as f:
        tree = f["Events"]
        events_base_mc.append(tree.arrays(options["branches"], library = "ak", how = "zip"))


events_base_data = awkward.concatenate(events_base_data)
events_base_mc = awkward.concatenate(events_base_mc)
logger.debug("Running over %d events in data" % len(events_base_data))
logger.debug("Running over %d events in mc" % len(events_base_mc))

# Don't have SampleManager yet, so create a dummy central weight branch for now
events_base_data["weight_central"] = awkward.ones_like(events_base_data["run"])
events_base_mc["weight_central"] = events_base_mc["genWeight"]


options_central = copy.deepcopy(options)
options_central["taggers"]["diphoton_tagger"]["photons"]["use_central_nano"] = True

options_custom = copy.deepcopy(options)
options_custom["taggers"]["diphoton_tagger"]["photons"]["use_central_nano"] = False

# Central nanoAOD branches
# Data
diphoton_tagger_central_data = DiphotonTagger(
        name = "diphoton_tagger_central_data",
        options = options_central["taggers"]["diphoton_tagger"]
)
tag_sequence_central_data = TagSequence(
        tag_list = [
            diphoton_tagger_central_data,
        ]
)
events_base_data_central = awkward.copy(events_base_data)
events_data_central, x = tag_sequence_central_data.run(events_base_data_central)

# MC
diphoton_tagger_central_mc = DiphotonTagger(
        name = "diphoton_tagger_central_mc",
        options = options_central["taggers"]["diphoton_tagger"]
)
tag_sequence_central_mc = TagSequence(
        tag_list = [
            diphoton_tagger_central_mc,
        ]
)
events_base_mc_central = awkward.copy(events_base_mc)
events_mc_central, x = tag_sequence_central_mc.run(events_base_mc_central)

# Added branches in custom nanoAOD
# Data
diphoton_tagger_custom_data = DiphotonTagger(
        name = "diphoton_tagger_custom_data",
        options = options_custom["taggers"]["diphoton_tagger"]
)
tag_sequence_custom_data = TagSequence(
        tag_list = [
            diphoton_tagger_custom_data,
        ]
)
events_base_data_custom = awkward.copy(events_base_data)
events_data_custom, x = tag_sequence_custom_data.run(events_base_data_custom)

# MC
diphoton_tagger_custom_mc = DiphotonTagger(
        name = "diphoton_tagger_custom_mc",
        options = options_custom["taggers"]["diphoton_tagger"]
)
tag_sequence_custom_mc = TagSequence(
        tag_list = [
            diphoton_tagger_custom_mc,
        ]
)
events_base_mc_custom = awkward.copy(events_base_mc)
events_mc_custom, x = tag_sequence_custom_mc.run(events_base_mc_custom)


eff_data_central = float(len(events_data_central["nominal"])) / float(len(events_base_data_central))
eff_mc_central = float(len(events_mc_central["nominal"])) / float(len(events_base_mc_central))

eff_data_custom = float(len(events_data_custom["nominal"])) / float(len(events_base_data_custom))
eff_mc_custom = float(len(events_mc_custom["nominal"])) / float(len(events_base_mc_custom))

logger.info("Efficiency of diphoton preselection on data, central vs. custom: %.4f vs. %.4f, ratio: %.3f" % (eff_data_central, eff_data_custom, eff_data_central / eff_data_custom))
logger.info("Efficiency of diphoton preselection on mc, central vs. custom: %.4f vs. %.4f, ratio: %.3f" % (eff_mc_central, eff_mc_custom, eff_mc_central / eff_mc_custom))

