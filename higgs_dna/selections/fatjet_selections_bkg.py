import awkward

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_FATJETS = {
        "pt" : 150.,
        "eta" : 2.4,
        # "Hqqqq_qqlv_vsQCDTop": -999
        "Hqqqq_vsQCDTop": -999 
}

def select_fatjets(fatjets, options, clean, name = "none", tagger = None): 
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_FATJETS,
        new = options
    )
    print("option is !!!!!!!!!!!!!", options)
    standard_cuts = object_selections.select_objects(fatjets, options, clean, name, tagger)
    # can apply some additional cut
    # add the H jet tagger for SL&FH channel((H3q+H4q+Hlvqq)/(H3q+H4q+Hlvqq+QCD+Top))
    # H_jet_cut = fatjets.Hqqqq_vsQCDTop > options["Hqqqq_vsQCDTop"]
    # H_jet_cut = fatjets.Hqqqq_qqlv_vsQCDTop > options["Hqqqq_qqlv_vsQCDTop"]
    # print("Hqqqq_qqlv_vsQCDTop cut value is : ",options["Hqqqq_qqlv_vsQCDTop"])
    # print("Hqqqq_vsQCDTop cut value is : ",options["Hqqqq_vsQCDTop"])
    # if options["tau2_tau1"] != 0:
    #     tau2_tau1_cut = ((fatjets.tau2 / fatjets.tau1) < options["tau2_tau1"])
    # else:
    #     tau2_tau1_cut = fatjets.pt > 0
    # print("tau2_tau1", options["tau2_tau1"])

    all_cuts = standard_cuts# & (H_jet_cut)

    if tagger is not None:
        tagger.register_cuts(
            names = ["standard_cuts", "all_cuts"],
            results = [standard_cuts, all_cuts],
            cut_type = name
        )

    return all_cuts
