import uproot
import awkward
import numpy
import vector
import json
import numba
import time

from higgs_dna.utils import setup_logger

file = "root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201021/test_nanoaod_1.root"

with open("../higgs_dna/taggers/metadata/diphoton_default.json", "r") as f_in:
    options = json.load(f_in)

with uproot.open(file) as f:
    tree = f["Events"]
    events = tree.arrays(options["branches"], library = "ak", how = "zip")

logger = setup_logger("DEBUG")

####################################################################################
# Start with a "dummy" diphoton preselection: pt > 25, at least 2 photons in event #
####################################################################################

photons = events.Photon

pt_cut = photons.pt > 25.
photons = photons[pt_cut]

n_pho_cut = awkward.num(photons) >= 2
photons = photons[n_pho_cut]

##########################################################################################
# Construct lead/sublead photon vectors with `vector`. Nicely interfaces with `awkward`. #
##########################################################################################

pho_vecs = vector.awk({
    "pt" : photons.pt,
    "eta" : photons.eta,
    "phi" : photons.phi,
    "mass" : photons.mass
})

# Vector construction is very intuitive
lead_photons = vector.awk({
    "pt" : photons[:,0].pt, 
    "eta" : photons[:,0].eta,
    "phi" : photons[:,0].phi,
    "mass" : photons[:,0].mass
})

sublead_photons = vector.awk({
    "pt" : photons[:,1].pt,
    "eta" : photons[:,1].eta,
    "phi" : photons[:,1].phi,
    "mass" : photons[:,1].mass
})
 
# Access properties of vectors the same way we would as for awkward arrays
logger.debug("lead_photons.pt : %s" % lead_photons.pt)
logger.debug("lead_photons[0].pt : %s" % lead_photons[0].pt)

######################################################
# Construct diphoton pairs from lead/sublead photons #
######################################################

# Operators are overloaded in intuitive ways
diphotons = lead_photons + sublead_photons

# Can access all the typical quantities we need
for i in range(3):
    logger.debug("diphoton %d, mass: %.2f, pt: %.2f, rapidity: %.2f" % (i, diphotons[i].mass, diphotons[i].pt, diphotons[i].rapidity))

# Can add additional fields to `diphotons` object
diphotons["lead_photon"] = lead_photons
diphotons["sublead_photon"] = sublead_photons

logger.debug("diphotons[0].lead_photon.pt : %.2f" % diphotons[0].lead_photon.pt)

#################################################
# Calculate delta_R between lead/sublead photon #
#################################################

# Could start from the original awkward arrays of vectors for the photons
delta_R_v1 = lead_photons.deltaR(sublead_photons)

# Or, could access them through the diphotons object
delta_R_v2 = diphotons.lead_photon.deltaR(diphotons.sublead_photon)

# Can add this as a field to the diphotons object
diphotons["delta_R"] = delta_R_v1

for i in range(3):
    logger.debug("dR(lead pho, sublead pho) : %.3f (method 1)" % (delta_R_v1[i]))
    logger.debug("dR(lead pho, sublead pho) : %.3f (method 2)" % (delta_R_v2[i]))
    logger.debug("dR(lead pho, sublead pho) : %.3f (read from dipho object)" % (diphotons[i].delta_R))

######################################
# Pre-compiling functions with numba #
######################################

@numba.njit
def compute_mass(photons, n_photons):
    n_events = len(photons)
    out = numpy.empty(n_events, numpy.float64)
    for i in range(n_events):
        total = vector.obj(px=0.0, py=0.0, pz=0.0, E=0.0)
        for j in range(n_photons[i]):
            pho = vector.obj(
                pt = photons[i][j].pt,
                eta = photons[i][j].eta,
                phi = photons[i][j].phi,
                mass = photons[i][j].mass
            )
            total = total + pho
        out[i] = total.mass
    return out

def compute_mass_slow(photons, n_photons):
    n_events = len(photons)
    out = numpy.empty(n_events, numpy.float64)
    for i in range(n_events):
        total = vector.obj(px=0.0, py=0.0, pz=0.0, E=0.0)
        for j in range(n_photons[i]):
            pho = vector.obj(
                pt = photons[i][j].pt,
                eta = photons[i][j].eta,
                phi = photons[i][j].phi,
                mass = photons[i][j].mass
            )
            total = total + pho
        out[i] = total.mass
    return out
 

photons = events.Photon
n_pho_cut = awkward.num(photons) == 2
photons = photons[n_pho_cut]

logger.debug("Number of events: %d" % (len(photons)))

n_photons = numpy.array(awkward.num(photons), numpy.int64)

# Check time (including compilation)
start = time.time()
compute_mass(photons, n_photons)
elapsed = time.time() - start

logger.debug("Time to compute m_gg (including compilation): %.5f" % (elapsed))

# Check time (already compiled)
start = time.time()
compute_mass(photons, n_photons)
elapsed = time.time() - start

logger.debug("Time to compute m_gg (already compiled): %.5f" % (elapsed))

# Check time (no numba)
start = time.time()
compute_mass_slow(photons, n_photons)
elapsed = time.time() - start

logger.debug("Time to compute m_gg without numba: %.5f" % (elapsed))

# Check time (fully vectorized way)
lead_photons = vector.awk({
    "pt" : photons[:,0].pt,
    "eta" : photons[:,0].eta,
    "phi" : photons[:,0].phi,
    "mass" : photons[:,0].mass
})

sublead_photons = vector.awk({
    "pt" : photons[:,1].pt,
    "eta" : photons[:,1].eta,
    "phi" : photons[:,1].phi,
    "mass" : photons[:,1].mass
})

start = time.time()
mass = (lead_photons + sublead_photons).mass
elapsed = time.time() - start

logger.debug("Time to compute m_gg (fully vectorized): %.5f" % (elapsed))
