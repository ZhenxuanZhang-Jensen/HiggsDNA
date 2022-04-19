import awkward
import vector
import numba
vector.register_awkward()

from higgs_dna.utils import awkward_utils

def select_x_to_yz(gen_part, x_pdgId, y_pdgId, z_pdgId):
    """
    Return all x->yy decays, sorted by x_pt.
    Values of `None` for any of `x_pdgId`, `y_pdgId`, or `z_pdgId` will result in no requirement on pdgId.
    For example, 
        - `x_pdgId = None`, `y_pdgId = 5`, `z_pdgId = 5` will select all X->bb decays
        - `x_pdgId = None`, `y_pdgId = 11`, `z_pdgId = None` will select all X->eY decays


    """

    if not isinstance(gen_part, vector.Vector4D):
        gen_part = awkward.Array(gen_part, with_name = "Momentum4D")

    if x_pdgId is None:
        x_pdgId = abs(gen_part.pdgId)
    if y_pdgId is None:
        y_pdgId = abs(gen_part.pdgId)
    if z_pdgId is None:
        z_pdgId = abs(gen_part.pdgId)

    gen_yz_from_x = gen_part[((abs(gen_part.pdgId) == y_pdgId) | (abs(gen_part.pdgId) == z_pdgId)) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == x_pdgId)]
    gen_yz_from_x = gen_yz_from_x[awkward.argsort(gen_yz_from_x.pt, ascending = False, axis = 1)]

    # Make all pairs of y's and z's that decay from x
    gen_child_pairs = awkward.combinations(gen_yz_from_x, 2, fields = ["LeadGenChild", "SubleadGenChild"])
    # Keep only the pairs where the y and the z decay from the same x
    gen_child_pairs = gen_child_pairs[gen_child_pairs.LeadGenChild.genPartIdxMother == gen_child_pairs.SubleadGenChild.genPartIdxMother] 
    
    # Grab the parent x
    gen_child_pairs["GenParent"] = gen_part[gen_child_pairs.LeadGenChild.genPartIdxMother]
    gen_child_pairs[("GenParent", "dR")] = gen_child_pairs.LeadGenChild.deltaR(gen_child_pairs.SubleadGenChild) 

    gen_child_pairs[("GenParent", "Child_1_Id")] = abs(gen_child_pairs.LeadGenChild.pdgId) 
    gen_child_pairs[("GenParent", "Child_2_Id")] = abs(gen_child_pairs.SubleadGenChild.pdgId) 

    if awkward.any(awkward.num(gen_child_pairs) >= 2):
        gen_child_pairs = gen_child_pairs[awkward.argsort(gen_child_pairs.GenParent.pt, ascending = False, axis = 1)] 

    return gen_child_pairs   
@numba.njit
def select_ww_to_qqqq(gen_part):
    """test for numba njit to improve the speed of loop in python"""
    #event loop
    W1_content = []
    W2_content = []
    H_content = []
    for i in range(len(gen_part)):
    # for i in range(1):
        list_quark_candidate = []
        #object loop
        for j in range(len(gen_part[i])):
            cut_Wplus = (abs(gen_part.pdgId[i][j]) <= 6 and (abs(gen_part.pdgId[i][gen_part.genPartIdxMother[i][j]]) == 24))
            if(cut_Wplus):
                # a = 1
                list_quark_candidate.append(gen_part[i][j])
        W1_candi = list_quark_candidate[0]+list_quark_candidate[1]
        W2_candi = list_quark_candidate[2]+list_quark_candidate[3]
        if(W1_candi.mass > W2_candi.mass):
            W1_candi = vector.obj(
                pt = (list_quark_candidate[0]+list_quark_candidate[1]).pt,
                eta = (list_quark_candidate[0]+list_quark_candidate[1]).eta,
                phi = (list_quark_candidate[0]+list_quark_candidate[1]).phi,
                mass = (list_quark_candidate[0]+list_quark_candidate[1]).mass
            )
            W2_candi = vector.obj(
                pt = (list_quark_candidate[2]+list_quark_candidate[3]).pt,
                eta = (list_quark_candidate[2]+list_quark_candidate[3]).eta,
                phi = (list_quark_candidate[2]+list_quark_candidate[3]).phi,
                mass = (list_quark_candidate[2]+list_quark_candidate[3]).mass
            )
        else:
            W2_candi = vector.obj(
                pt = (list_quark_candidate[0]+list_quark_candidate[1]).pt,
                eta = (list_quark_candidate[0]+list_quark_candidate[1]).eta,
                phi = (list_quark_candidate[0]+list_quark_candidate[1]).phi,
                mass = (list_quark_candidate[0]+list_quark_candidate[1]).mass
            )
            W1_candi = vector.obj(
                pt = (list_quark_candidate[2]+list_quark_candidate[3]).pt,
                eta = (list_quark_candidate[2]+list_quark_candidate[3]).eta,
                phi = (list_quark_candidate[2]+list_quark_candidate[3]).phi,
                mass = (list_quark_candidate[2]+list_quark_candidate[3]).mass
            )
        H_candi = vector.obj(px = 0., py = 0., pz = 0., E = 0.) # IMPORTANT NOTE: you need to initialize this to an empty vector first. Otherwise, you will get ZeroDivisionError exceptions for like 1 out of a million events (seemingly only with numba). 
        H_candi = H_candi + W1_candi + W2_candi
        W1_content.append([
            W1_candi.pt,
            W1_candi.eta,
            W1_candi.phi,
            W1_candi.mass,
        ])
        W2_content.append([
            W2_candi.pt,
            W2_candi.eta,
            W2_candi.phi,
            W2_candi.mass,
        ])
        H_content.append([
            H_candi.pt,
            H_candi.eta,
            H_candi.phi,
            H_candi.mass,
        ])

    return W1_content,W2_content,H_content