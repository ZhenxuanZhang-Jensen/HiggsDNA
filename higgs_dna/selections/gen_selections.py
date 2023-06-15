import awkward
import vector
import numba
import numpy
import logging
logger = logging.getLogger(__name__)

vector.register_awkward()

from higgs_dna.utils import awkward_utils
def select_t3_4(events):
    gen_part = awkward.Array(events.GenPart,with_name="Momentum4D")
    gen_qqqq = gen_part[(abs(gen_part.pdgId)<= 6) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24) ]
    gen_gg = gen_part[(abs(gen_part.pdgId) == 22) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 25) ]
    return None

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
def select_fake_photon(events,diphotons):
    gen_part = awkward.Array(diphotons.GenPart,with_name="Momentum4D")
    # gen_Part = awkward_utils.add_field(
    #     events = diphotons,
    #     name = "genpart",
    #     data=gen_part

    # )
    # awkward_utils.add_object_fields(
    # events=diphotons,
    # name="GenPart",
    # objects=gen_Part,
    # n_objects=1,
    # dummy_value=-999
    # )
    leadidx = awkward.unflatten(diphotons.LeadPhoton.genPartIdx,counts=awkward.to_list(awkward.full_like(diphotons.LeadPhoton.genPartIdx,1)),axis=-1)
    gen_matched_leadphoton = (abs(gen_part.pdgId[leadidx])==22) & (diphotons.LeadPhoton.genPartFlav == 1)
    gen_fake_leadphoton = ~gen_matched_leadphoton 
    is_leadphoton = awkward.flatten(gen_matched_leadphoton)
    fake_leadphoton = awkward.flatten(gen_fake_leadphoton)
    diphotons[("Diphoton", "is_LeadPhoton")]=is_leadphoton
    diphotons[("LeadPhoton", "is_LeadPhoton")]=is_leadphoton
    diphotons[("Diphoton","GenPart_pdgId")]=gen_part.pdgId
    diphotons[("Diphoton","GenPart_eta")]=gen_part.eta
    diphotons[("Diphoton","GenPart_pt")]=gen_part.pt
    diphotons[("Diphoton","GenPart_phi")]=gen_part.phi
    diphotons[("Diphoton","GenPart_mass")]=gen_part.mass
    diphotons[("Diphoton","GenPart_IdxMother")]=gen_part.genPartIdxMother
    diphotons[("Diphoton","GenPart_status")]=gen_part.status
    diphotons[("Diphoton","GenPart_statusFlags")]=gen_part.statusFlags
    test=awkward.where(diphotons.LeadPhoton.is_LeadPhoton == True, diphotons.LeadPhoton, -999)

    # real_leadphoton = awkward.zip(
    # {
    # "pt": awkward.unflatten(awkward.where(is_leadphoton == True, diphotons, -999),1).pt,
    # "eta": awkward.unflatten(awkward.where(is_leadphoton == True, diphotons, -999),1).eta,
    # "phi": awkward.unflatten(awkward.where(is_leadphoton == True, diphotons, -999),1).phi,
    # "mass": awkward.unflatten(awkward.where(is_leadphoton == True, diphotons, -999),1).mass
    # }), 
    prompt_leadpho = awkward_utils.add_field(
    events = diphotons,
    name = "real_leadphoton",
    data = awkward.unflatten(test,1)
    )
    awkward_utils.add_object_fields(
    events=diphotons,
    name="Prompt_leadphoton",
    objects=prompt_leadpho,
    n_objects=1,
    dummy_value=-999
    )
    fake_leadpho = awkward_utils.add_field(
    events = diphotons,
    name = "Fake_leadphoton",
    data = awkward.unflatten(awkward.where(fake_leadphoton == True, diphotons, -999),1)
    )
    awkward_utils.add_object_fields(
    events=diphotons,
    name="Fake_leadphoton",
    objects=fake_leadpho,
    n_objects=1,
    dummy_value=-999
    )
    subleadidx = awkward.unflatten(diphotons.LeadPhoton.genPartIdx,counts=awkward.to_list(awkward.full_like(diphotons.SubleadPhoton.genPartIdx,1)),axis=-1)
    gen_matched_subleadphoton = (abs(gen_part.pdgId[subleadidx])==22) & (diphotons.SubleadPhoton.genPartFlav == 1)
    gen_fake_subleadphoton = ~gen_matched_subleadphoton 
    is_subphoton = awkward.flatten(gen_matched_subleadphoton)
    fake_subleadphoton = awkward.flatten(gen_fake_subleadphoton)
    diphotons[("Diphoton", "is_SubleadPhoton")]=is_subphoton
    prompt_subleadpho = awkward_utils.add_field(
    events = diphotons,
    name = "real_subleadphoton",
    data = awkward.unflatten(awkward.where(is_subphoton == True, diphotons, -999),1)
    )
    awkward_utils.add_object_fields(
    events=diphotons,
    name="Prompt_subleadphoton",
    objects=prompt_subleadpho,
    n_objects=1,
    )
    print(fake_leadpho)
    fake_subleadpho = awkward_utils.add_field(
    events = diphotons,
    name = "Fake_subleadphoton",
    data = awkward.unflatten(awkward.where(fake_subleadphoton == True, diphotons, -999),1)
    )
    awkward_utils.add_object_fields(
    events=diphotons,
    name="Fake_subleadphoton",
    objects=fake_subleadpho,
    n_objects=1,
    )
    fake_pho = awkward_utils.add_field(
    events = diphotons,
    name = "Fake_photon",
    data = awkward.concatenate([fake_leadpho,fake_subleadpho],axis=1)
    )
    # awkward_utils.add_object_fields(
    # events=events,
    # name="Fake_photon",
    # objects=fake_pho,
    # n_objects=1,
    # )
    prompt_pho = awkward_utils.add_field(
    events = diphotons,
    name = "Prompt_photon",
    data = awkward.concatenate([prompt_leadpho,prompt_subleadpho],axis=1)
    )
    # awkward_utils.add_object_fields(
    # events=events,
    # name="Prompt_photon",
    # objects=prompt_pho,
    # n_objects=1,
    # )
    return fake_pho,prompt_pho
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
def gen_Hww_2q2l(events):
    gen_part = awkward.Array(events.GenPart,with_name="Momentum4D")
    # --------- gen level 2 signal quarks and 2 signal photons and 2 leptons and the Higgs from gg -------- #
    gen_qq = gen_part[(abs(gen_part.pdgId)<= 6) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24) ]
    gen_gg = gen_part[(abs(gen_part.pdgId) == 22) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 25) ]
    # gen_lv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=11)&(abs(gen_part.pdgId)<=14))]
    gen_lv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=11)&(abs(gen_part.pdgId)<=14))]
    dummy = awkward.full_like(gen_gg[0],99)
    ngen_lv=awkward.num(gen_lv,axis=-1)
    gen_lvcut=(ngen_lv==2)

   # gen_lvcut=(ngen_lv!=4)&(ngen_lv!=0)
    #no4lep_lv type:bool
    # -------There are some events which contains 4 leptons at genlevel, 2 e+e- and 2 leptons ---------#
    # -------Some events contain which contain 0 lepton at genlevel----------#
    logger.debug("origin num of events:%s, after gen level cut >=2 lepton num of events:%s, ==2 lepton num of events:%s)"%(str(len(events)),str(len(gen_lv[ngen_lv>=2])),str(len(gen_lv[ngen_lv==2]))))
    gen_lv = gen_lv.mask[gen_lvcut]
    gen_qq = gen_qq.mask[gen_lvcut]
    gen_gg = gen_gg.mask[gen_lvcut]

   # gen_part=awkward.fill_none(gen_part,dummy,axis=0)
    gen_qq = awkward.fill_none(gen_qq,dummy,axis=0)
    gen_gg = awkward.fill_none(gen_gg,dummy,axis=0)
    gen_lv = awkward.fill_none(gen_lv,dummy,axis=0)
    unflatten_gen_q1 = awkward.unflatten(gen_qq[:,0],1)
    unflatten_gen_q2 = awkward.unflatten(gen_qq[:,1],1)
    unflatten_gen_l1 = awkward.unflatten(gen_lv[:,0],1)
    unflatten_gen_v1 = awkward.unflatten(gen_lv[:,1],1)
    unflatten_gen_pho1 = awkward.unflatten(gen_gg[:,0],1)
    unflatten_gen_pho2 = awkward.unflatten(gen_gg[:,1],1)
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
    gen_pho1_p4 = vector.awk(
        {
            "pt" : unflatten_gen_pho1["pt"],
            "eta" : unflatten_gen_pho1["eta"],
            "phi" : unflatten_gen_pho1["phi"],
            "mass" : unflatten_gen_pho1["mass"]
        },
        with_name = "Momentum4D"
    )
    gen_pho2_p4 = vector.awk(
        {
            "pt" : unflatten_gen_pho2["pt"],
            "eta" : unflatten_gen_pho2["eta"],
            "phi" : unflatten_gen_pho2["phi"],
            "mass" : unflatten_gen_pho2["mass"]
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
    logger.debug("import select_ww_qqlv W1,W2,Hcandi")
    Wqq_candi,Wll_candi,H_candi = select_ww_to_qqlv(gen_qq,gen_lv)

    W1_candi4D = awkward.zip(
    {
    "pt": numpy.array(Wqq_candi)[:,0],
    "eta": numpy.array(Wqq_candi)[:,1],
    "phi": numpy.array(Wqq_candi)[:,2],
    "mass": numpy.array(Wqq_candi)[:,3],
    }, 
    with_name="Momentum4D")
    W2_candi4D = awkward.zip(
    {
    "pt": numpy.array(Wll_candi)[:,0],
    "eta": numpy.array(Wll_candi)[:,1],
    "phi": numpy.array(Wll_candi)[:,2],
    "mass": numpy.array(Wll_candi)[:,3],
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
def select_ww_to_qqlv(gen_qq,gen_lv):    
    W1_content = []
    W2_content = []
    H_content = []
    for i in range(len(gen_lv)):
        list_quark_candidate = []
        list_lepton_candidate = []
        for j in range(len(gen_lv[i])):
            list_quark_candidate.append(gen_qq[i][j])
            list_lepton_candidate.append(gen_lv[i][j])
        W1_candi = list_quark_candidate[0]+list_quark_candidate[1]

        W2_candi = list_lepton_candidate[0]+list_lepton_candidate[1]

        W1_candi = vector.obj(
                pt = (list_quark_candidate[0]+list_quark_candidate[1]).pt,
                eta = (list_quark_candidate[0]+list_quark_candidate[1]).eta,
                phi = (list_quark_candidate[0]+list_quark_candidate[1]).phi,
                mass = (list_quark_candidate[0]+list_quark_candidate[1]).mass)
        W2_candi = vector.obj(
                pt = (list_lepton_candidate[0]+list_lepton_candidate[1]).pt,
                eta = (list_lepton_candidate[0]+list_lepton_candidate[1]).eta,
                phi = (list_lepton_candidate[0]+list_lepton_candidate[1]).phi,
                mass = (list_lepton_candidate[0]+list_lepton_candidate[1]).mass)         
        H_candi = vector.obj(px = 0., py = 0., pz = 0., E = 0.)
        #zero=numpy.array([0.,0.,0.,0.])
        #zero=zero.astype(numpy.float32)
        #H_candi = vector.obj(pt = zero[0], eta =zero[1], phi = zero[2], mass = zero[3])
        ##H_candi = vector.obj(pt = 0., eta = 0., phi = 0., E = 0.)

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
    return W1_content, W2_content,H_content
def gen_3categories(events):
    gen_part = awkward.Array(events.GenPart,with_name="Momentum4D")
    gen_lv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=11)&(abs(gen_part.pdgId)<=16))]
    print(type(gen_lv))
    nlep = awkward.num(gen_lv, axis = -1)
    #nofourlep_cut = (nlep!=4)
    gen_part=gen_part[nlep!=4]
    events = events[nlep!=4]
    # remove four lepton events (~200 events)
    # --------- gen level 2 signal quarks and 2 signal photons and 2 leptons and the Higgs from gg -------- #
    gen_qq = gen_part[(abs(gen_part.pdgId)<= 6) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24) ]
    gen_gg = gen_part[(abs(gen_part.pdgId) == 22) & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 25) ]
    gen_lv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=11)&(abs(gen_part.pdgId)<=16))]
    gen_μv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=13)&(abs(gen_part.pdgId)<=14))]
    gen_ev = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=11)&(abs(gen_part.pdgId)<=12))]
    gen_τv = gen_part[(abs(gen_part.pdgId[gen_part.genPartIdxMother]) == 24)&((abs(gen_part.pdgId)>=15)&(abs(gen_part.pdgId)<=16))]

    # record the initial events num of three categories: muon, electron, tau
    n_muon_event = len(gen_μv[awkward.num(gen_μv)==2])
    n_electron_event = len(gen_ev[awkward.num(gen_ev)==2])
    n_tau_event = len(gen_τv[awkward.num(gen_τv)==2])
    e_events = events[awkward.num(gen_ev,axis=-1)!=0]
    μ_events = events[awkward.num(gen_μv,axis=-1)!=0]
    τ_events = events[awkward.num(gen_τv,axis=-1)!=0]
    logger.debug("At the beginning, there are %s muon lepton candidates, %s electron lepton candidates, %s tau lepton candidates events"%(n_muon_event, n_electron_event,n_tau_event))    
    return events, n_muon_event, n_electron_event, n_tau_event, e_events, μ_events, τ_events