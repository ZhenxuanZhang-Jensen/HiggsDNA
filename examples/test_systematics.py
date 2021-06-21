import uproot
import json
import awkward
import numpy
import time
import matplotlib.pyplot as plt

from higgs_dna.taggers.diphoton_tagger import DiphotonTagger
from higgs_dna.taggers.tag_sequence import TagSequence

from higgs_dna.utils.logger_utils import setup_logger

file = "root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201021/test_nanoaod_1.root"

with open("../higgs_dna/taggers/metadata/diphoton_default.json", "r") as f_in:
    options = json.load(f_in)

with uproot.open(file) as f:
    tree = f["Events"]
    events = tree.arrays(options["branches"], library = "ak", how = "zip")


### Explicit syst method ###
# Create dummy systematics for independent collection (IC)
photon_pt_up = events.Photon.pt + awkward.ones_like(events.Photon.pt)
photon_pt_down = events.Photon.pt - awkward.ones_like(events.Photon.pt)

logger = setup_logger("DEBUG")

x = [1, 2, 5, 10, 25, 50, 100]
y_explicit = []
y_implicit = []

for N_systs in x: 
    ICs = { "nominal" : events}
    for i in range(N_systs):
        ICs["photons_syst%d_up" % i] = awkward.with_field(events, photon_pt_up, ["Photon", "pt"])
        ICs["photons_syst%d_down" % i] = awkward.with_field(events, photon_pt_down, ["Photon", "pt"])


    start = time.time()
    for key, evts in ICs.items():
        start = time.time()
        pt_cut = evts.Photon.pt > 30
        elapsed_explicit = time.time() - start


    ### Implicit syst method ###
    events_implicit = awkward.copy(events)
    syst_array = numpy.ones(1 + (2*N_systs))
    events_implicit["Photon", "pt"] = events_implicit.Photon.pt * syst_array[None, None, :] 
    print(events_implicit.Photon.pt[0])

    # Compare performance
    start = time.time()
    pt_cut = events.Photon.pt > 30
    elapsed_single = time.time() - start

    start = time.time()
    pt_cut = events_implicit.Photon.pt > 30
    elapsed_implicit = time.time() - start

    start = time.time()
    for key, evts in ICs.items():
        pt_cut = evts.Photon.pt > 30
    elapsed_explicit = time.time() - start

    logger.info("For %d systs" % N_systs)
    logger.info("Time to perform loop through single collection: %.6f s" % elapsed_single)
    logger.info("Time to perform loop through ICs (explicit method): %.6f s (%.2fx slower)" % (elapsed_explicit, elapsed_explicit / elapsed_single))
    logger.info("Time to perform vectorized computation over ICs (implicit method): %.6f s (%.2fx slower)" % (elapsed_implicit, elapsed_implicit / elapsed_single))
    logger.info("Implicit method is a factor of %.2f faster\n\n" % (elapsed_explicit / elapsed_implicit))

    y_explicit.append(elapsed_explicit / elapsed_single)
    y_implicit.append(elapsed_implicit / elapsed_single)


fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(x, y_explicit, label = "explicit method", color = "blue")
ax1.plot(x, y_implicit, label = "implicit method", color = "red")

plt.xlabel("Number of photon systematics")
plt.ylabel("Time")
ax1.legend(loc = "upper left")
plt.savefig("syst_comparison.pdf")
