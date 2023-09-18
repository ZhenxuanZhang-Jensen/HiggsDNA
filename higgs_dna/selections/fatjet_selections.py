import awkward

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_FATJETS = {
        "pt" : 150.,
        "eta" : 2.4,
        # "Hqqqq_qqlv_vsQCDTop": -999
        "Hqqqq_vsQCDTop": -999 
}

def select_fatjets(fatjets, subjets, options, clean, name = "none", tagger = None): 
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_FATJETS,
        new = options
    )
    standard_cuts = object_selections.select_objects(fatjets, options, clean, name, tagger)
    # can apply some additional cut
    # add the H jet tagger for SL&FH channel((H3q+H4q+Hlvqq)/(H3q+H4q+Hlvqq+QCD+Top))
    H_jet_cut = fatjets.Hqqqq_vsQCDTop > options["Hqqqq_vsQCDTop"]
    # H_jet_cut = fatjets.Hqqqq_qqlv_vsQCDTop > options["Hqqqq_qqlv_vsQCDTop"]
    # print("Hqqqq_qqlv_vsQCDTop cut value is : ",options["Hqqqq_qqlv_vsQCDTop"])
    print("Hqqqq_vsQCDTop cut value is : ",options["Hqqqq_vsQCDTop"])
    H_jet_cut=awkward.fill_none(awkward.pad_none(H_jet_cut,1,axis=1),False,axis=-1)
    subjets_cut = ((awkward.fill_none(awkward.pad_none(subjets.pt,1,axis=1),-999,axis=-1)[awkward.fill_none(awkward.pad_none(fatjets.subJetIdx1,1,axis=1),0,axis=-1)]>20)==True)&((awkward.fill_none(awkward.pad_none(subjets.pt,1,axis=1),-999,axis=-1)[awkward.fill_none(awkward.pad_none(fatjets.subJetIdx2,1,axis=1),0,axis=-1)]>20)==True)
    
    standard_cuts=awkward.fill_none(awkward.pad_none(standard_cuts,1,axis=1),False,axis=-1)
    all_cuts = standard_cuts & (H_jet_cut) & (subjets_cut)
    # all_cuts = standard_cuts & (H_jet_cut)

    if tagger is not None:
        tagger.register_cuts(
            names = ["standard_cuts", "H_jet_cut","SubJet_cut", "all_cuts"],
            results = [standard_cuts, H_jet_cut,subjets_cut, all_cuts],
            cut_type = name
        )

    return all_cuts