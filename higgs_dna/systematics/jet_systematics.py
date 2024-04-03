import awkward
import numpy

from correctionlib import _core

from higgs_dna.utils import misc_utils, awkward_utils

###################################
### b-tag continuous reshape SF ###
###################################

BTAG_RESHAPE_SF_FILE = {
    "2016" : "jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json", 
    "2016UL_preVFP" : "jsonpog-integration/POG/BTV/2016preVFP_UL/btagging.json", 
    "2016UL_postVFP" : "jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json", 
    "2017" : "jsonpog-integration/POG/BTV/2017_UL/btagging.json",
    "2018" : "jsonpog-integration/POG/BTV/2018_UL/btagging.json"
}

DEEPJET_RESHAPE_SF = {
    "2016" : "deepJet_shape", 
    "2016UL_preVFP" : "deepJet_shape",
    "2016UL_postVFP" : "deepJet_shape",
    "2017" : "deepJet_shape",
    "2018" : "deepJet_shape"
}

DEEPJET_VARIATIONS = { 
    "up_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm (4)
    "up_lf" : [5],
    "up_hfstats1" : [5],
    "up_hfstats2" : [5],
    "up_cferr1" : [4],
    "up_cferr2" : [4],
    "up_hf" : [0],
    "up_lfstats1" : [0],
    "up_lfstats2" : [0],
    "down_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm(4)
    "down_lf" : [5],
    "down_hfstats1" : [5],
    "down_hfstats2" : [5],
    "down_cferr1" : [4],
    "down_cferr2" : [4],
    "down_hf" : [0],
    "down_lfstats1" : [0],
    "down_lfstats2" : [0],
}

def btag_deepjet_reshape_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_bjets_Run2_UL/
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/btvExample.py

    Note: application of SFs should not change the overall normalization of a sample (before any b-tagging selection) and each sample should be adjusted by an overall weight derived in a phase space with no requirements on b-jets such that the normalization is unchanged. TODO: link BTV TWiki that describes this.
    """
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "hadronFlavour"), (input_collection, "btagDeepFlavB") 
    ]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
   
    jets = events[input_collection]
    jets["flavor"] = jets.hadronFlavour

    # Flatten jets then convert to numpy for compatibility with correctionlib
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.49999 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        99999999.
    )
    jet_disc = awkward.to_numpy(jets_flattened.btagDeepFlavB)        

    variations_list = ["central"]
    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}

    central_sf = evaluator[DEEPJET_RESHAPE_SF[year]].evalv(
            "central",
            jet_flavor,
            jet_abs_eta,
            jet_pt,
            jet_disc
    )

    variations["central"] = awkward.unflatten(central_sf, n_jets)

    for var in variations_list:
        if var == "central":
            continue
        applicable_flavors = DEEPJET_VARIATIONS[var] # the up/down variations are only applicable to specific flavors of jet
        var_sf = central_sf 
        for f in applicable_flavors:
            var_sf = numpy.where(
                jet_flavor == f,
                evaluator[DEEPJET_RESHAPE_SF[year]].evalv(
                    var,
                    numpy.ones_like(jet_flavor) * f,
                    jet_abs_eta,
                    jet_pt,
                    jet_disc
                ),
                var_sf
            )

        variations[var] = awkward.unflatten(var_sf, n_jets) # make jagged again

    for var in variations.keys():
        # Set SFs = 1 for jets which are not applicable (pt <= 20 or |eta| >= 2.5)
        variations[var] = awkward.where(
                (jets.pt <= 20.0) | (abs(jets.eta) >= 2.5),
                awkward.ones_like(variations[var]),
                variations[var]
        )

    return variations


def dummy_jes_syst(events, is_data):
    """
    Dummy function illustrating a jet energy scale uncertainty that results in new jet collections with Jet.pt varied.
    Should be deleted once real examples are implemented.
    """
    jets = events.Jet 

    variations = {}
    variations["central"] = jets.pt + (2 * awkward.ones_like(jets.pt))
    if not is_data:
        variations["up"] = jets.pt + (12 * awkward.ones_like(jets.pt))
        variations["down"] = jets.pt - (8 * awkward.ones_like(jets.pt))
    return variations

def WvsQCD_medium_jes_syst(event,year):
    """
    See:
        - https://github.com/cms-jet/ParticleNetSF/blob/ParticleNet_TopW_SFs_NanoV9/Csv/W_MD_Run2_SF.csv

    Note: Particle Net WvsQCD Nanov9 SFs are applied to the W jets in the event.
    """
    weight_PNet_WvsQCDW1_central = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW2_central = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW1_up = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW2_up = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW1_down = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW2_down = awkward.ones_like(event.category)
    W1ptbin1=(event["fatjet_W_1_pt"]>200 )&( event["fatjet_W_1_pt"]<=300 )
    W1ptbin2=(event["fatjet_W_1_pt"]>300 )&( event["fatjet_W_1_pt"]<=400 )
    W1ptbin3=(event["fatjet_W_1_pt"]>400 )&( event["fatjet_W_1_pt"]<=800 )
    W2ptbin1=(event["fatjet_W_2_pt"]>200 )&( event["fatjet_W_2_pt"]<=300 )
    W2ptbin2=(event["fatjet_W_2_pt"]>300 )&( event["fatjet_W_2_pt"]<=400 )
    W2ptbin3=(event["fatjet_W_2_pt"]>400 )&( event["fatjet_W_2_pt"]<=800 )
    SF_W1_ptbin1=((event.category==1)|(event.category==7)) & W1ptbin1
    SF_W1_ptbin2=((event.category==1)|(event.category==7)) & W1ptbin2
    SF_W1_ptbin3=((event.category==1)|(event.category==7)) & W1ptbin3
    SF_W2_ptbin1=(event.category==6) & W2ptbin1
    SF_W2_ptbin2=(event.category==6) & W2ptbin2
    SF_W2_ptbin3=(event.category==6) & W2ptbin3
    if year=="2016UL_preVFP":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.904261, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.936451, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.872640, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.875718, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.910996, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.841575, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.902272, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.969975, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.836782, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.904261, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.936451, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.872640, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.875718, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.910996, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.841575, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.902272, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.969975, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.836782, weight_PNet_WvsQCDW2_down)
    elif year=="2016UL_postVFP":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.887065, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.934490, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.843186, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.866233, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.904293, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.829220, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.786696, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.854186, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.720922, weight_PNet_WvsQCDW1_down)   
             
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.887065, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.934490, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.843186, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.866233, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.904293, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.829220, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.786696, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.854186, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.720922, weight_PNet_WvsQCDW2_down)        
    elif year=="2017":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.918743, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.961281, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.878062, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.898412, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.923404, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.873408, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.903940, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.949428, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.858823, weight_PNet_WvsQCDW1_down)  
              
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.918743, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.961281, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.878062, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.898412, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.923404, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.873408, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.903940, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.949428, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.858823, weight_PNet_WvsQCDW2_down)        
    elif year=="2018":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.863467, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.883999, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.843701, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.861461, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.883182, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.840135, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.815265, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.949428, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.853422, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.863467, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.883999, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.843701, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.861461, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.883182, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.840135, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.815265, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.949428, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.853422, weight_PNet_WvsQCDW2_down)
    event["weight_PNet_WvsQCDW1_up"]=weight_PNet_WvsQCDW1_up
    event["weight_PNet_WvsQCDW1_central"]=weight_PNet_WvsQCDW1_central
    event["weight_PNet_WvsQCDW1_down"]=weight_PNet_WvsQCDW1_down
    
    event["weight_PNet_WvsQCDW2_up"]=weight_PNet_WvsQCDW2_up
    event["weight_PNet_WvsQCDW2_central"]=weight_PNet_WvsQCDW2_central
    event["weight_PNet_WvsQCDW2_down"]=weight_PNet_WvsQCDW2_down
    return event
def WvsQCD_loose_jes_syst(event,year):
    """
    See:
        - https://github.com/cms-jet/ParticleNetSF/blob/ParticleNet_TopW_SFs_NanoV9/Csv/W_MD_Run2_SF.csv

    Note: Particle Net WvsQCD Nanov9 SFs are applied to the W jets in the event.
    """
    weight_PNet_WvsQCDW1_central = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW2_central = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW1_up = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW2_up = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW1_down = awkward.ones_like(event.category)
    weight_PNet_WvsQCDW2_down = awkward.ones_like(event.category)
    W1ptbin1=(event["fatjet_W_1_pt"]>200 )&( event["fatjet_W_1_pt"]<=300 )
    W1ptbin2=(event["fatjet_W_1_pt"]>300 )&( event["fatjet_W_1_pt"]<=400 )
    W1ptbin3=(event["fatjet_W_1_pt"]>400 )&( event["fatjet_W_1_pt"]<=800 )
    W2ptbin1=(event["fatjet_W_2_pt"]>200 )&( event["fatjet_W_2_pt"]<=300 )
    W2ptbin2=(event["fatjet_W_2_pt"]>300 )&( event["fatjet_W_2_pt"]<=400 )
    W2ptbin3=(event["fatjet_W_2_pt"]>400 )&( event["fatjet_W_2_pt"]<=800 )
    SF_W1_ptbin1=((event.category==1)|(event.category==7)) & W1ptbin1
    SF_W1_ptbin2=((event.category==1)|(event.category==7)) & W1ptbin2
    SF_W1_ptbin3=((event.category==1)|(event.category==7)) & W1ptbin3
    SF_W2_ptbin1=(event.category==6) & W2ptbin1
    SF_W2_ptbin2=(event.category==6) & W2ptbin2
    SF_W2_ptbin3=(event.category==6) & W2ptbin3
    if year=="2016UL_preVFP":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9195379553371, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.955807058879442, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.853997286203603, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9463957277711702, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.981673394066279, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.911701826690222, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9398161345034564, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*1.003032920932323, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.879995031378807, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9195379553371, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.955807058879442, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.853997286203603, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9463957277711702, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.981673394066279, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.911701826690222, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9398161345034564, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*1.003032920932323, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.879995031378807, weight_PNet_WvsQCDW2_down)
    elif year=="2016UL_postVFP":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9450578338842212, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.995323626432011, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.88139954116051, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9092346655321338, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.94696027685031, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.872197862861603, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.8925475131701243, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.965027816841837, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.823833325878663, weight_PNet_WvsQCDW1_down)   
             
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9450578338842212, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.995323626432011, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.88139954116051, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9092346655321338, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.94696027685031, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.872197862861603, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.8925475131701243, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.965027816841837, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.823833325878663, weight_PNet_WvsQCDW2_down)        
    elif year=="2017":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.982730281923659, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*1.018665093373757, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.946056549531482, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9594886915963721, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.98423246182761, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.935112279809698, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9892794381973279, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.944301786656573, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*1.036482653739614, weight_PNet_WvsQCDW1_down)  
              
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.982730281923659, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*1.018665093373757, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.946056549531482, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9594886915963721, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.98423246182761, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.935112279809698, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9892794381973279, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.944301786656573, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*1.036482653739614, weight_PNet_WvsQCDW2_down)        
    elif year=="2018":
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.9190298705570141, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.94124308044617, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.899896452243954, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.91497012004329, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.935173054970369, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.893927977587993, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW1_central = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_central)*0.8647360917703466, weight_PNet_WvsQCDW1_central)
        weight_PNet_WvsQCDW1_up = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_up)*0.905865333361337, weight_PNet_WvsQCDW1_up)
        weight_PNet_WvsQCDW1_down = awkward.where(SF_W1_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW1_down)*0.826121510983113, weight_PNet_WvsQCDW1_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.9190298705570141, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.94124308044617, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin1, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.899896452243954, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.91497012004329, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.935173054970369, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin2, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.893927977587993, weight_PNet_WvsQCDW2_down)
        
        weight_PNet_WvsQCDW2_central = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_central)*0.8647360917703466, weight_PNet_WvsQCDW2_central)
        weight_PNet_WvsQCDW2_up = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_up)*0.905865333361337, weight_PNet_WvsQCDW2_up)
        weight_PNet_WvsQCDW2_down = awkward.where(SF_W2_ptbin3, awkward.ones_like(weight_PNet_WvsQCDW2_down)*0.826121510983113, weight_PNet_WvsQCDW2_down)
    event["weight_PNet_WvsQCDW1_up"]=weight_PNet_WvsQCDW1_up
    event["weight_PNet_WvsQCDW1_central"]=weight_PNet_WvsQCDW1_central
    event["weight_PNet_WvsQCDW1_down"]=weight_PNet_WvsQCDW1_down
    
    event["weight_PNet_WvsQCDW2_up"]=weight_PNet_WvsQCDW2_up
    event["weight_PNet_WvsQCDW2_central"]=weight_PNet_WvsQCDW2_central
    event["weight_PNet_WvsQCDW2_down"]=weight_PNet_WvsQCDW2_down
    return event
