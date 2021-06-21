### This script is intended as an example to illustrate the tag sequence and systematics workflow. It will become outdated once we implement SamplesManager, JobsManager, etc.

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
file = "root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201021/test_nanoaod_IDX.root"
idx = numpy.arange(1,3)

events_base = []
for i in idx:
    file_idx = file.replace("IDX", str(i))
    print(file_idx)
    with uproot.open(file_idx) as f:
        tree = f["Events"]
        events_base.append(tree.arrays(options["branches"], library = "ak", how = "zip"))

events_base = awkward.concatenate(events_base) # TODO: look into awkward.lazy and/or awkward.iterate as alternatives
logger.debug("Running over %d events" % len(events_base))

# Notes:
#   - in the future, we should have a SamplesManager class
#   - user should be able to specify their samples and files like:
"""
sample_list = [
    "ttH_M125",
    "Run2_Data"
    ...
]

sample_manager = SampleManager(
    samples = sample_list,
    years = ["2016", "2017", "2018"]
)
"""
#   - each of the strings in sample_list point to entries in a comprehensive samples.json file:
"""
{
    "ttH_M125" : {
        "2018" : {
            "das_name" : "/ttHJetToGG_M125_TuneCP5_13TeV_amcatnloFXFX-madspin-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM",
            "xs" : 100.0
        },
        "2022" : {
            "das_name" : ...,
            "xs" : 110.0
        },
    ...
"""
#   - for computing scale1fb's, we should keep track of the sum of weights for every file we successfully run over and calculate scale1fb and scale weights appropriately at the very end
#       - this way, you can easily proceed with your analysis if 1 or 2 MC files are failing due to e.g. xrootd issues


# Don't have SampleManager yet, so create a dummy central weight branch for now
events_base["weight_central"] = events_base["genWeight"]

###########################################################
### Step 2: define the systematics you want to consider ###
###########################################################

# Create dummy systematics
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
                    "module" : lepton_systematics, # make sure this module is imported!
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
                    "module" : lepton_systematics,
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
            "diphoton_presel_sf" : {
                "type" : "object",
                "method" : "from_function",
                "function" : {
                    "module" : photon_systematics,
                    "name" : "photon_preselection_sf",
                },
                "modify_central_weight" : True,
                "input_collection" : "Photon",
                "target_collections" : [
                    ("Diphoton", "LeadPhoton"),
                    ("Diphoton", "SubleadPhoton")
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
                    "module" : photon_systematics,
                    "name" : "dummy_photon_pt_syst"
                }
            },
            "dummy_jes_syst" : {
                "method" : "from_function",
                "branch_modified" : ("Jet", "pt"),
                "function" : {
                    "module" : jet_systematics,
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

# Create SystematicsProducer instance to manage the calculation of systematics
syst_producer = SystematicsProducer(options = syst_options)

# NOTE:
# SystematicsProducer gets run twice: once before running taggers and once after
# In the first pass (before taggers) it computes all systematic variations which do not require the output of taggers (e.g. theory weights, photon/jet energy scale variations)
# In the second pass (after taggers) it computes remaining systematic variations which do require the output of taggers (e.g. diphoton preselection SF, lepton SFs for selected leptons in each tagger)
# In the second pass (after taggers) it computes remaining systematic variations which do require the output of taggers (e.g. diphoton preselection SF, lepton SFs for selected leptons in each tagger)
# As of now, I can only envision weight variations (and not independent collections) depending on the output of a tagger

#################################a
### Step 3: define TagSequence ###
##################################

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


# Run systematics producer and tag sequence

# Step 3.1: run systematics producer before any taggers: compute weight variations that don't require tagger outputs, compute independent collections
events = syst_producer.produce(events_base)

# Step 3.2: run the tag sequence
tag_sequence = TagSequence(
        tag_list = [
            diphoton_tagger,
            [ # list within the list means these are enforced to be orthogonal
                tth_tagger, # first priority
                thq_tagger  # second priority
            ]
        ]
)
selected_events, tag_idx_map = tag_sequence.run(events)

# Step 3.3: run systematics producer after tag sequence: compute weight variations that require tagger output
selected_events = syst_producer.apply_remaining_weight_systs(selected_events, tag_idx_map)
summary = syst_producer.summary

# Notes:
#   - Step 3 will change once we consider many files split over many jobs
#   - I envision the syntax to look something like this:
"""
samples = SamplesManager(
    samples = samples_list,
    years = ["2016", "2017", "2018", "2022"]
)

syst_producer = SystematicsProducer(options = syst_options)

tag_sequence = TagSequence(
    tag_list = ... # as before
)

jobs_manager = JobsManager(
    sample_manager = samples,
    systematics_producer = syst_producer,
    tag_sequence = tag_sequence,
    output_format = "parquet",
    output_variables = ["Diphoton_mass", "Diphoton_pt"],
    batch = "condor" # or "local", "dask", etc.
)

jobs_manager.run() # launch and monitor all jobs
    # for each job and corresponding set of files + metadata, we will do steps 3.1-3.3
    # we will also do step 4 (below), writing events to disk

jobs_manager.summarize() # write a json with summary info, including:
    # details of each job
    # metadata for each sample: e.g. various cut efficiencies for taggers, summary info for systematics, etc.
"""

#############################
### Step 4: write to disk ###
#############################

# Below is just temporary code. Will be handled by a JobsManager in the future, as described above.

weight_branches = list(summary["weights"].keys()) 

ext = "_toy_test" # extension to identify output files
os.system("mkdir -p output/")
for syst_tag, syst_events in selected_events.items():
    base_branches = ["genWeight", "weight_central", "tag_idx", "run", ("Diphoton", "mass"), ("Diphoton", "pt")]
    if syst_tag == "nominal":
        save_branches = base_branches + weight_branches
    else:
        save_branches = base_branches
        for branch in weight_branches:
            if "central" in branch:
                save_branches.append(branch)

    save_map = {}
    for branch in save_branches:
        if isinstance(branch, tuple):
            save_name = "_".join(branch)
        else:
            save_name = branch
        save_map[save_name] = syst_events[branch]

    syst_events = awkward.zip(save_map)

    awkward.to_parquet(syst_events, "output/%s%s.parquet" % (syst_tag, "_toy_test"))

# Now, explore output files
files = glob.glob("output/*%s.parquet" % ext)
for file in files:
    logger.info("\nLoading file %s" % file)
    array = awkward.from_parquet(file)
    df = awkward.to_pandas(array)
    print(df["Diphoton_mass"])
    print(df["Diphoton_pt"])
    print(df["tag_idx"])
    #if "nominal" in file:
    #    print(df["weight_diphoton_presel_sf_central"])
    #    print(df["weight_diphoton_presel_sf_up"])
    branches = df.columns
    logger.info("Branches: %s" % branches)
    tags = numpy.unique(df["tag_idx"])
    for tag in tags:
        logger.info("%d events for tag %d" % (len(df[df["tag_idx"] == tag]), tag))

array = awkward.from_parquet("output/nominal_toy_test.parquet")
df = awkward.to_pandas(array)

tth_events = df[df["tag_idx"] == 0]
thq_events = df[df["tag_idx"] == 1]

logger.info("tth events, mean value of gen weight: %.4f" % (tth_events["genWeight"].mean()))
logger.info("thq events, mean value of gen weight: %.4f" % (thq_events["genWeight"].mean()))
logger.info("tth events, mean value of central weight: %.4f" % (tth_events["weight_central"].mean()))
logger.info("thq events, mean value of central weight: %.4f" % (thq_events["weight_central"].mean()))
