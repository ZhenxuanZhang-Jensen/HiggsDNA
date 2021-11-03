import awkward
import vector
import numpy
import numba

vector.register_awkward()

from higgs_dna.utils import awkward_utils

def select_x_to_yz(gen_part, x_pdgId, y_pdgId, z_pdgId):
    """
    Return all x->yy decays, sorted by x_pt
    """


    gen_yz_from_x = gen_part[((abs(gen_part.pdgId) == y_pdgId) | (abs(gen_part.pdgId) == z_pdgId)) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == x_pdgId)]
    gen_yz_from_x = gen_yz_from_x[awkward.argsort(gen_yz_from_x.pt, ascending = False, axis = 1)]

    # Make all pairs of y's and z's that decay from x
    gen_child_pairs = awkward.combinations(gen_yz_from_x, 2, fields = ["LeadGenChild", "SubleadGenChild"])
    # Keep only the pairs where the y and the z decay from the same x
    gen_child_pairs = gen_child_pairs[gen_child_pairs.LeadGenChild.genPartIdxMother == gen_child_pairs.SubleadGenChild.genPartIdxMother] 
    
    # Grab the parent x
    gen_child_pairs["GenParent"] = gen_part[gen_child_pairs.LeadGenChild.genPartIdxMother]

    lead_gen_child_p4 = vector.awk({
        "pt" : gen_child_pairs.LeadGenChild.pt,
        "eta" : gen_child_pairs.LeadGenChild.eta,
        "phi" : gen_child_pairs.LeadGenChild.phi,
        "mass" : gen_child_pairs.LeadGenChild.mass
    })

    sublead_gen_child_p4 = vector.awk({
        "pt" : gen_child_pairs.SubleadGenChild.pt,
        "eta" : gen_child_pairs.SubleadGenChild.eta,
        "phi" : gen_child_pairs.SubleadGenChild.phi,
        "mass" : gen_child_pairs.SubleadGenChild.mass
    })

    gen_child_pairs[("GenParent", "dR")] = lead_gen_child_p4.deltaR(sublead_gen_child_p4) 

    if awkward.any(awkward.num(gen_child_pairs) >= 2):
        gen_child_pairs = gen_child_pairs[awkward.argsort(gen_child_pairs.GenParent.pt, ascending = False, axis = 1)] 

    return gen_child_pairs    
