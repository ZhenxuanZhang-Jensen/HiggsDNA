import awkward
import vector
import numpy
import numba

vector.register_awkward()

from higgs_dna.utils import awkward_utils

def select_x_to_yy(gen_particles, x_pdgId, y_pdgId):
    """
    Find the x -> yy decay in the event.
    For now, assumes there is at most 1 x->yy, can be generalized to >1 in the future.
    """
    gen_x_from_y = gen_particles[(abs(gen_particles.pdgId) == x_pdgId) & (abs(gen_particles.pdgId[gen_particles.genPartIdxMother]) == y_pdgId)]
    gen_x = gen_particles[gen_x_from_y.genPartIdxMother]

    gen_x_from_y = gen_x_from_y[awkward.argsort(gen_x_from_y.pt, ascending = False, axis = 1)]
    gen_x_to_yy = awkward.combinations(gen_x_from_y, 2, fields = ["LeadGenY", "SubleadGenY"])

    gen_x_to_yy["GenX"] = gen_particles[gen_x_to_yy.LeadGenY.genPartIdxMother]
    
    lead_gen_y_p4 = vector.awk({
        "pt" : gen_x_to_yy.LeadGenY.pt,
        "eta" : gen_x_to_yy.LeadGenY.eta,
        "phi" : gen_x_to_yy.LeadGenY.phi,
        "mass" : gen_x_to_yy.LeadGenY.mass
    })
    sublead_gen_y_p4 = vector.awk({
        "pt" : gen_x_to_yy.SubleadGenY.pt,
        "eta" : gen_x_to_yy.SubleadGenY.eta,
        "phi" : gen_x_to_yy.SubleadGenY.phi,
        "mass" : gen_x_to_yy.SubleadGenY.mass
    })
    gen_x_to_yy["dR"] = lead_gen_y_p4.deltaR(sublead_gen_y_p4)
    gen_x_to_yy["n_GenX"] = awkward.num(gen_x_to_yy.GenX)

    return gen_x_to_yy

def select_x(gen_particles, pdgId, status_flags = None):
    x = gen_particles[(abs(gen_particles.pdgId) == pdgId)]
    if status_flags is not None:
        x = x[x.statusFlags == status_flags]
    x = x[awkward.argsort(x.pt, ascending = False, axis = 1)]
    return x
