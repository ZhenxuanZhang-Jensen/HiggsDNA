import awkward
import vector
import numba
import numpy
vector.register_awkward()

from higgs_dna.utils import awkward_utils
import logging
logger = logging.getLogger(__name__)

def gen_Hww_4q(events):
    gen_part = awkward.Array(events.GenPart,with_name="Momentum4D")
    # -------------- gen level 4 signal quarks and 2 signal photons and the Higgs from gg -------------- #
    gen_qqqq = gen_part[(abs(gen_part.pdgId)<= 6) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24) ]
    gen_gg = gen_part[(abs(gen_part.pdgId) == 22) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 25) ]
    unflatten_gen_q1 = awkward.unflatten(gen_qqqq[:,0],1)
    unflatten_gen_q2 = awkward.unflatten(gen_qqqq[:,1],1)
    unflatten_gen_q3 = awkward.unflatten(gen_qqqq[:,2],1)
    unflatten_gen_q4 = awkward.unflatten(gen_qqqq[:,3],1)

    # the four quarks Momentum4D

    gen_q1_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q1["pt"],
            "eta" : unflatten_gen_q1["eta"],
            "phi" : unflatten_gen_q1["phi"],
            "mass" : unflatten_gen_q1["mass"]
        },
        with_name = "Momentum4D"
    )
    gen_q2_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q2["pt"],
            "eta" : unflatten_gen_q2["eta"],
            "phi" : unflatten_gen_q2["phi"],
            "mass" : unflatten_gen_q2["mass"]
        },
        with_name = "Momentum4D"
    )
    gen_q3_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q3["pt"],
            "eta" : unflatten_gen_q3["eta"],
            "phi" : unflatten_gen_q3["phi"],
            "mass" : unflatten_gen_q3["mass"]
        },
        with_name = "Momentum4D"
    )
    gen_q4_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q4["pt"],
            "eta" : unflatten_gen_q4["eta"],
            "phi" : unflatten_gen_q4["phi"],
            "mass" : unflatten_gen_q4["mass"]
        },
        with_name = "Momentum4D"
    )

    gen_q1 = awkward_utils.add_field(
    events = events,
    name = "gen_q1",
    data = unflatten_gen_q1
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GEN_q1",
    objects=gen_q1,
    n_objects=1,
    dummy_value=-999
    )
    gen_q2 = awkward_utils.add_field(
    events = events,
    name = "gen_q2",
    data = unflatten_gen_q2
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GEN_q2",
    objects=gen_q2,
    n_objects=1,
    dummy_value=-999
    )
    gen_q3 = awkward_utils.add_field(
    events = events,
    name = "gen_q3",
    data = unflatten_gen_q3
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GEN_q3",
    objects=gen_q3,
    n_objects=1,
    dummy_value=-999
    )
    gen_q4 = awkward_utils.add_field(
    events = events,
    name = "gen_q4",
    data = unflatten_gen_q4
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GEN_q4",
    objects=gen_q4,
    n_objects=1,
    dummy_value=-999
    )

    # ------------ gen level two W boson order with mass and the Higgs from WW ------------ #
    W1_candi,W2_candi,H_candi = select_ww_to_qqqq(gen_part)
    W1_candi4D = awkward.zip(
    {
    "pt": numpy.array(W1_candi)[:,0],
    "eta": numpy.array(W1_candi)[:,1],
    "phi": numpy.array(W1_candi)[:,2],
    "mass": numpy.array(W1_candi)[:,3],
    }, 
    with_name="Momentum4D")
    W2_candi4D = awkward.zip(
    {
    "pt": numpy.array(W2_candi)[:,0],
    "eta": numpy.array(W2_candi)[:,1],
    "phi": numpy.array(W2_candi)[:,2],
    "mass": numpy.array(W2_candi)[:,3],
    }, 
    with_name="Momentum4D")
    H_candi4D = awkward.zip(
    {
    "pt": numpy.array(H_candi)[:,0],
    "eta": numpy.array(H_candi)[:,1],
    "phi": numpy.array(H_candi)[:,2],
    "mass": numpy.array(H_candi)[:,3],
    }, 
    with_name="Momentum4D")
    unflatten_W1_candi4D = awkward.unflatten(W1_candi4D,1)
    unflatten_W2_candi4D = awkward.unflatten(W2_candi4D,1)
    unflatten_H_candi4D = awkward.unflatten(H_candi4D,1)
    genW1_qq = awkward_utils.add_field(
    events=events,
    name="genW1_qq",
    data=unflatten_W1_candi4D
    )
    genW2_qq = awkward_utils.add_field(
    events=events,
    name="genW2_qq",
    data=unflatten_W2_candi4D 
    )
    genHWW_qqqq = awkward_utils.add_field(
    events=events,
    name="genHWW_qqqq",
    data=unflatten_H_candi4D
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GENW1_qq",
    objects=genW1_qq,
    n_objects=1,
    dummy_value=-999
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GENW2_qq",
    objects=genW2_qq,
    n_objects=1,
    dummy_value=-999
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GENHWW_qqqq",
    objects=genHWW_qqqq,
    n_objects=1,
    dummy_value=-999
    )
    return gen_q1_p4, gen_q2_p4, gen_q3_p4, gen_q4_p4


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
def produce_genmatched(events,jets):
    padded_jets=jets
    jet_1_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,0] 
    jet_2_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,1] 
    jet_3_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,2] 
    jet_4_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,3] 
    jet_5_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,4] 
    jet_6_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,5] 
    jet_7_deltaR_q1 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q1"],7,axis=-1),999,axis=-1)[:,6]   
    jet_1_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,0] 
    jet_2_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,1] 
    jet_3_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,2] 
    jet_4_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,3] 
    jet_5_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,4] 
    jet_6_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,5] 
    jet_7_deltaR_q2 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q2"],7,axis=-1),999,axis=-1)[:,6]  
    jet_1_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,0] 
    jet_2_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,1] 
    jet_3_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,2] 
    jet_4_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,3] 
    jet_5_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,4] 
    jet_6_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,5] 
    jet_7_deltaR_q3 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q3"],7,axis=-1),999,axis=-1)[:,6]   
    jet_1_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,0] 
    jet_2_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,1] 
    jet_3_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,2] 
    jet_4_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,3] 
    jet_5_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,4] 
    jet_6_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,5] 
    jet_7_deltaR_q4 = awkward.fill_none(awkward.pad_none(padded_jets["deltaR_q4"],7,axis=-1),999,axis=-1)[:,6]
    genmatch = numpy.array([
    jet_1_deltaR_q1,jet_1_deltaR_q2,jet_1_deltaR_q3,jet_1_deltaR_q4,
    jet_2_deltaR_q1,jet_2_deltaR_q2,jet_2_deltaR_q3,jet_2_deltaR_q4,
    jet_3_deltaR_q1,jet_3_deltaR_q2,jet_3_deltaR_q3,jet_3_deltaR_q4,
    jet_4_deltaR_q1,jet_4_deltaR_q2,jet_4_deltaR_q3,jet_4_deltaR_q4,
    jet_5_deltaR_q1,jet_5_deltaR_q2,jet_5_deltaR_q3,jet_5_deltaR_q4,
    jet_6_deltaR_q1,jet_6_deltaR_q2,jet_6_deltaR_q3,jet_6_deltaR_q4,
    jet_7_deltaR_q1,jet_7_deltaR_q2,jet_7_deltaR_q3,jet_7_deltaR_q4])
    return genmatch, padded_jets, jet_1_deltaR_q1

#@numba.njit
def genmatched_fourjets(genmatch,jet_1_deltaR_q1): 
    genmatched_cut=[]
   # logger.debug("len(jet_1_deltaR_q1): %s" % (len(jet_1_deltaR_q1)))
    for i in range(len(jet_1_deltaR_q1)):
        jet_qqqqdR = numpy.array(genmatch[:,i])
        if len((jet_qqqqdR[jet_qqqqdR!=999]))>=16:
            gen_cut=[False,False,False,False,False,False,False]
            jet_qqqqdR=numpy.reshape(jet_qqqqdR,(int(len(jet_qqqqdR)/4),4))
            mindR=numpy.argmin(jet_qqqqdR)
            jet1_index=int(mindR/4)
            q1_index=mindR%4
            jet_qqqqdR[:,q1_index]=[999,999,999,999,999,999,999]
            jet_qqqqdR[jet1_index,:]=[999,999,999,999]
            mindR=numpy.argmin(jet_qqqqdR)
            jet2_index=int(mindR/4)
            q2_index=mindR%4
            jet_qqqqdR[:,q2_index]=[999,999,999,999,999,999,999]
            jet_qqqqdR[jet2_index,:]=[999,999,999,999]
            mindR=numpy.argmin(jet_qqqqdR)
            jet3_index=int(mindR/4)
            q3_index=mindR%4
            jet_qqqqdR[:,q3_index]=[999,999,999,999,999,999,999]
            jet_qqqqdR[jet3_index,:]=[999,999,999,999]
            mindR=numpy.argmin(jet_qqqqdR)
            jet4_index=int(mindR/4)
            q4_index=mindR%4
            jet_qqqqdR[:,q4_index]=[999,999,999,999,999,999,999]
            jet_qqqqdR[jet4_index,:]=[999,999,999,999]
            gen_cut[jet1_index]=True
            gen_cut[jet2_index]=True
            gen_cut[jet3_index]=True
            gen_cut[jet4_index]=True
        else:
            gen_cut=[False,False,False,False,False,False,False]
        genmatched_cut=numpy.append(genmatched_cut,gen_cut)
    genmatched_cut=numpy.reshape(genmatched_cut,(int(len(genmatched_cut)/7),7))
    genmatched_cut=genmatched_cut.astype(bool)
    unmatched_cut=~genmatched_cut
    genmatched_cut=genmatched_cut.astype(int)
    unmatched_cut=unmatched_cut.astype(int)
 #   logger.debug("genmatched_cut produced")

    return genmatched_cut, unmatched_cut
def gen_Hww_2q2l(events):
    gen_part = awkward.Array(events.GenPart,with_name="Momentum4D")
    # --------- gen level 2 signal quarks and 2 signal photons and 2 leptons and the Higgs from gg -------- #
    gen_qq = gen_part[(abs(gen_part.pdgId)<= 6) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24) ]
    gen_gg = gen_part[(abs(gen_part.pdgId) == 22) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 25) ]
    gen_lv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=11)&(abs(gen_part.pdgId)<=14))]  
    ngen_lv=awkward.num(gen_lv,axis=-1)
    gen_lvcut=(ngen_lv!=4)&(ngen_lv!=0)
    #no4lep_lv type:bool
    # -------There are some events which contains 4 leptons at genlevel, 2 e+e- and 2 leptons ---------#
    # -------Some events contain which contain 0 lepton at genlevel----------#
    logger.debug("origin num of events:%s, after gen level cut num of events:%s)")
    gen_lv = gen_lv.mask[gen_lvcut]
    gen_qq = gen_qq.mask[gen_lvcut]
    gen_gg = gen_gg.mask[gen_lvcut]
   # gen_part = gen_part.mask[gen_lvcut]
    #gen_part = gen_part[gen_lvcut]
    dummy = awkward.Array({
    'eta':[-999,-999],
     'phi':[-999,-999],
     'pt':[-999,-999],
    'mass':[-99,-999],
     'genPartIdxMother':[-999,-999],
     'pdgId':[-999,-999],
     'status':[-999,-999],
     'statusFlags':[-999,-999],
    }
    )
   # gen_part=awkward.fill_none(gen_part,dummy,axis=0)
    gen_qq = awkward.fill_none(gen_part,dummy,axis=0)
    gen_gg = awkward.fill_none(gen_gg,dummy,axis=0)
    gen_lv = awkward.fill_none(gen_lv,dummy,axis=0)
    unflatten_gen_q1 = awkward.unflatten(gen_qq[:,0],1)
    unflatten_gen_q2 = awkward.unflatten(gen_qq[:,1],1)
    unflatten_gen_l1 = awkward.unflatten(gen_lv[:,0],1)
    unflatten_gen_v1 = awkward.unflatten(gen_lv[:,1],1)
    # the four quarks Momentum4D
    gen_W1_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q1["pt"]+unflatten_gen_q2["pt"],
            "eta" : unflatten_gen_q1["eta"]+unflatten_gen_q2["eta"],
            "phi" : unflatten_gen_q1["phi"]+unflatten_gen_q2["phi"],
            "mass" : unflatten_gen_q1["mass"]+unflatten_gen_q2["mass"]
        },
        with_name = "Momentum4D"
    )    
    gen_l1_p4 = vector.awk(
        {
            "pt" : unflatten_gen_l1["pt"],
            "eta" : unflatten_gen_l1["eta"],
            "phi" : unflatten_gen_l1["phi"],
            "mass" : unflatten_gen_l1["mass"]
        },
        with_name = "Momentum4D"
    )       
    gen_q1_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q1["pt"],
            "eta" : unflatten_gen_q1["eta"],
            "phi" : unflatten_gen_q1["phi"],
            "mass" : unflatten_gen_q1["mass"]
        },
        with_name = "Momentum4D"
    )
    gen_q2_p4 = vector.awk(
        {
            "pt" : unflatten_gen_q2["pt"],
            "eta" : unflatten_gen_q2["eta"],
            "phi" : unflatten_gen_q2["phi"],
            "mass" : unflatten_gen_q2["mass"]
        },
        with_name = "Momentum4D"
    )
    
    gen_l1 = awkward_utils.add_field(
    events = events,
    name = "gen_l1",
    data = unflatten_gen_l1
    )
    #genW1_qq["deltaR_l1"] = 
    gen_l1["deltaR_Wqq"] = gen_l1_p4.deltaR(gen_W1_p4)
    awkward_utils.add_object_fields(
    events = events,
    name = "GEN_l1",
    objects = gen_l1,
    n_objects=1,
    dummy_value=-999,
    )
    gen_v1 = awkward_utils.add_field(
    events = events,
    name = "gen_v1",
    data = unflatten_gen_v1
    )
    awkward_utils.add_object_fields(
    events = events,
    name = "GEN_v1",
    objects = gen_v1,
    n_objects=1,
    dummy_value=-999,
    )
    gen_q1 = awkward_utils.add_field(
    events = events,
    name = "gen_q1",
    data = unflatten_gen_q1
    )
    gen_q1["deltaR_l1"] = gen_q1_p4.deltaR(gen_l1_p4)
    gen_q1["deltaR_q2"] = gen_q1_p4.deltaR(gen_q2_p4)
    awkward_utils.add_object_fields(
    events=events,
    name="GEN_q1",
    objects=gen_q1,
    n_objects=1,
    dummy_value=-999
    )
    gen_q2 = awkward_utils.add_field(
    events = events,
    name = "gen_q2",
    data = unflatten_gen_q1
    )
    gen_q2["deltaR_l1"] = gen_q2_p4.deltaR(gen_l1_p4)
    awkward_utils.add_object_fields(
    events=events,
    name="GEN_q2",
    objects=gen_q2,
    n_objects=1,
    dummy_value=-999
    )
    W1_candi,W2_candi,H_candi = select_ww_to_qqlv(gen_part)
    W1_candi4D = awkward.zip(
    {
    "pt": numpy.array(W1_candi)[:,0],
    "eta": numpy.array(W1_candi)[:,1],
    "phi": numpy.array(W1_candi)[:,2],
    "mass": numpy.array(W1_candi)[:,3],
    }, 
    with_name="Momentum4D")
    W2_candi4D = awkward.zip(
    {
    "pt": numpy.array(W2_candi)[:,0],
    "eta": numpy.array(W2_candi)[:,1],
    "phi": numpy.array(W2_candi)[:,2],
    "mass": numpy.array(W2_candi)[:,3],
    }, 
    with_name="Momentum4D")
    H_candi4D = awkward.zip(
    {
    "pt": numpy.array(H_candi)[:,0],
    "eta": numpy.array(H_candi)[:,1],
    "phi": numpy.array(H_candi)[:,2],
    "mass": numpy.array(H_candi)[:,3],
    }, 
    with_name="Momentum4D")
    unflatten_W1_candi4D = awkward.unflatten(W1_candi4D,1)
    unflatten_W2_candi4D = awkward.unflatten(W2_candi4D,1)
    unflatten_H_candi4D = awkward.unflatten(H_candi4D,1)
   # logger.debug("num(W1_candi)=%s"%len(W1_candi))
   # logger.debug("num(W2_candi)=%s"%len(W2_candi))
   # logger.debug("W1_candi%s"%W1_candi)
 #   W1_candi = awkward.pad_none(W1_candi,len(events),axis=0)

    logger.debug("num(W1_candi)=%s"%len(W1_candi4D))
    logger.debug("(W1_candi)=%s"% W1_candi4D)

    genW1_qq = awkward_utils.add_field(
    events=events,
    name="genW1_qq",
    data=unflatten_W1_candi4D
    )
    genW2_lv = awkward_utils.add_field(
    events=events,
    name="genW2_lv",
    data=unflatten_W2_candi4D 
    )
    genHWW_qqlv = awkward_utils.add_field(
    events=events,
    name="genHWW_qqlv",
    data=unflatten_H_candi4D
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GENW1_qq",
    objects=genW1_qq,
    n_objects=1,
    dummy_value=-999
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GENW2_lv",
    objects=genW2_lv,
    n_objects=1,
    dummy_value=-999
    )
    awkward_utils.add_object_fields(
    events=events,
    name="GENHWW_qqlv",
    objects=genHWW_qqlv,
    n_objects=1,
    dummy_value=-999
    )   
    return gen_l1_p4, gen_q1_p4,gen_q2_p4


@numba.njit
def select_ww_to_qqlv(gen_part):
    W1_content = []
    W2_content = []
    H_content = []
    for i in range(len(gen_part)):
    # for i in range(1):
        list_quark_candidate = []
        list_lepton_candidate = []
        #object loop
        for j in range(len(gen_part[i])):
            cut_W2 = (abs(gen_part.pdgId[i][j] <= 14) and abs(gen_part.pdgId[i][j]>=11) and (abs(gen_part.pdgId[i][gen_part.genPartIdxMother[i][j]])==24))
            cut_W1 = (abs(gen_part.pdgId[i][j]) <= 6 )and (abs(gen_part.pdgId[i][gen_part.genPartIdxMother[i][j]]) == 24)
            if(cut_W1):
                # a = 1
                list_quark_candidate.append(gen_part[i][j])
            if(cut_W2):
                list_lepton_candidate.append(gen_part[i][j])

        if len(list_lepton_candidate)!=2:
            W1_candi = vector.obj(
                pt = -999,
                eta = -999,
                phi = -999,
                mass = -999
            ) 
        if len(list_quark_candidate)!=2:   
            W1_candi = vector.obj(
                pt = -999,
                eta = -999,
                phi = -999,
                mass = -999
            )    

        else:
            W1_candi = list_quark_candidate[0]+list_quark_candidate[1]
            W1_candi = vector.obj(
                    pt = (list_quark_candidate[0]+list_quark_candidate[1]).pt,
                    eta = (list_quark_candidate[0]+list_quark_candidate[1]).eta,
                    phi = (list_quark_candidate[0]+list_quark_candidate[1]).phi,
                    mass = (list_quark_candidate[0]+list_quark_candidate[1]).mass
                )

            W2_candi = list_lepton_candidate[0]+list_lepton_candidate[1]
            W2_candi = vector.obj(
                    pt = (list_lepton_candidate[0]+list_lepton_candidate[1]).pt,
                    eta = (list_lepton_candidate[0]+list_lepton_candidate[1]).eta,
                    phi = (list_lepton_candidate[0]+list_lepton_candidate[1]).phi,
                    mass = (list_lepton_candidate[0]+list_lepton_candidate[1]).mass
                )
            H_candi = vector.obj(px = 0., py = 0., pz = 0., E = 0.)
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
    return W1_content, W2_content, H_content








