import awkward
import vector

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
