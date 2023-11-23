import awkward
import numpy

from correctionlib import _core

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.utils import systematic_from_bins

ELECTRON_ID_SF_FILE = {
    "2016" : "jsonpog-integration/POG/EGM/2016postVFP_UL/electron.json",
    "2016UL_preVFP" : "jsonpog-integration/POG/EGM/2016preVFP_UL/electron.json",
    "2016UL_postVFP" : "jsonpog-integration/POG/EGM/2016postVFP_UL/electron.json",
    "2017" : "jsonpog-integration/POG/EGM/2017_UL/electron.json",
    "2018" : "jsonpog-integration/POG/EGM/2018_UL/electron.json"
}
MUON_ID_SF_FILE = {
    "2016" : "jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json",
    "2016UL_preVFP" : "jsonpog-integration/POG/MUO/2016preVFP_UL/muon_Z.json",
    "2016UL_postVFP" : "jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json",
    "2017" : "jsonpog-integration/POG/MUO/2017_UL/muon_Z.json",
    "2018" : "jsonpog-integration/POG/MUO/2018_UL/muon_Z.json"
}

ELECTRON_ID_SF = {
    "2016" : "2016postVFP",
    "2016UL_preVFP" : "2016preVFP",
    "2016UL_postVFP" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018"
}
MUON_ID_SF = {
    "2016" : "2016preVFP_UL",
    "2016UL_preVFP" : "2016preVFP_UL",
    "2016UL_postVFP" : "2016postVFP_UL",
    "2017" : "2017_UL",
    "2018" : "2018"
}

def electron_id_sf(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(ELECTRON_ID_SF_FILE[year]))

    electrons = events[input_collection]

    # Flatten electrons then convert to numpy for compatibility with correctionlib
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    ele_eta = numpy.clip(
        awkward.to_numpy(electrons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )

    ele_pt = numpy.clip(
        awkward.to_numpy(electrons_flattened.pt),
        10.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    sf = evaluator["UL-Electron-ID-SF"].evalv(
            ELECTRON_ID_SF[year],
            "sf",
            working_point,
            ele_eta,
            ele_pt
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst_vars = ["sfup", "sfdown"] 
        for syst_var in syst_vars:
            syst = evaluator["UL-Electron-ID-SF"].evalv(
                    ELECTRON_ID_SF[year],
                    syst_var,
                    working_point,
                    ele_eta,
                    ele_pt
            )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_electrons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (electrons.pt < 10.0) | (electrons.pt >= 500.0) | (abs(electrons.eta) >= 2.5),
                awkward.ones_like(variations[var]),
                variations[var]
        )

    return variations




def muon_tightid_tightiso_sf(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/MUO_2017_UL_muon_Z.html
    """

    required_fields = [
         (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)
    logger.debug(year)
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(MUON_ID_SF_FILE[year]))
    # evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(MUON_ID_SF_FILE["2016UL_preVFP"]))

    muon = events[input_collection]

    # Flatten muon then convert to numpy for compatibility with correctionlib
    n_muon = awkward.num(muon)
    muon_flattened = awkward.flatten(muon)
    muon_pt = numpy.clip(
        awkward.to_numpy(muon_flattened.pt),
        15.0, # SFs only valid for pT >= 15.0
        2000 # and pT < Inf.
    )
    muon_eta = numpy.abs(numpy.clip(
        awkward.to_numpy(muon_flattened.eta),
        -2.39999,
        2.39999 # SFs only valid up to eta 2.4
    ))
    # Calculate SF and syst
    variations = {}


    sf = evaluator["NUM_TightRelIso_DEN_TightIDandIPCut"].evalv(
            MUON_ID_SF[year],
            muon_eta,
            muon_pt,
            "sf",
    )
    variations["central"] = awkward.unflatten(sf, n_muon)

    if not central_only:
        syst_vars = ["systup", "systdown"] 
        for syst_var in syst_vars:
            syst = evaluator["NUM_TightRelIso_DEN_TightIDandIPCut"].evalv(
                    MUON_ID_SF[year],
                    muon_eta,
                    muon_pt,
                    syst_var
            )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_muon)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (muon.pt < 15.0) |(abs(muon.eta) >= 2.4),
                awkward.ones_like(variations[var]),
                variations[var]
        )
    return variations






DUMMY_LEPTON_SFs = {
    "variables" : ["lepton_pt", "lepton_eta"],
    "bins" : [
        {
            "lepton_pt" : [25.0, 50.0],
            "lepton_eta" : [0.0, 1.5],
            "value" : 0.5,
            "uncertainty" : 0.25
        },
        {
            "lepton_pt" : [25.0, 50.0],
            "lepton_eta" : [1.5, 999],
            "value" : 1.5,
            "uncertainty" : 0.25
        },
        {
            "lepton_pt" : [50.0, 999999.],
            "lepton_eta" : [0.0, 999999.],
            "value" : 0.5,
            "uncertainty" : 0.25
        },
 
    ]
}

def dummy_lepton_sf(events, central_only, input_collection):
    """
    Dummy function illustrating a lepton scale factor.
    Should be deleted once real examples are implemented.
    """
    required_fields = [
            (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    leptons = events[input_collection]
    variations = systematic_from_bins(
            bins = DUMMY_LEPTON_SFs,
            variables = {
                "lepton_pt" : leptons.pt,
                "lepton_eta" : leptons.eta
            },
            central_only = central_only
    )

    return variations 

def muon_id_sfSTAT(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ID_SF

    if not central_only:
        variations["up"] = leptons.ID_SF+leptons.ID_SFstat
        variations["down"] = leptons.ID_SF-leptons.ID_SFstat

    return variations

def muon_id_sfSYS(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ID_SF

    if not central_only:
        variations["up"] = leptons.ID_SF+leptons.ID_SFsyst
        variations["down"] = leptons.ID_SF-leptons.ID_SFsyst

    return variations

def muon_iso_sfSTAT(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ISO_SF

    if not central_only:
        variations["up"] = leptons.ISO_SF+leptons.ISO_SFstat
        variations["down"] = leptons.ISO_SF-leptons.ISO_SFstat

    return variations

def muon_iso_sfSYS(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ISO_SF

    if not central_only:
        variations["up"] = leptons.ISO_SF+leptons.ISO_SFsyst
        variations["down"] = leptons.ISO_SF-leptons.ISO_SFsyst

    return variations

def tauJ_sf(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.sfDeepTau2017v2p1VSjet_Loose_ext

    if not central_only:
        variations["up"] = leptons.sfDeepTau2017v2p1VSjet_LooseUp_ext
        variations["down"] = leptons.sfDeepTau2017v2p1VSjet_LooseDown_ext

    return variations

def tauM_sf(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.sfDeepTau2017v2p1VSmu_VLoose_ext

    if not central_only:
        variations["up"] = leptons.sfDeepTau2017v2p1VSmu_VLooseUp_ext
        variations["down"] = leptons.sfDeepTau2017v2p1VSmu_VLooseDown_ext

    return variations

def tauE_sf(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.sfDeepTau2017v2p1VSe_VVLoose_ext

    if not central_only:
        variations["up"] = leptons.sfDeepTau2017v2p1VSe_VVLooseUp_ext
        variations["down"] = leptons.sfDeepTau2017v2p1VSe_VVLooseDown_ext #NB according to Tau POG we should be using Tight for the vsJ scale factor to be valid!! check twiki 

    return variations
def highptmuonsf(event):
    weight_muon_highptid_central = awkward.ones_like(event.category)
    weight_muon_highptid_up = awkward.ones_like(event.category)
    weight_muon_highptid_down = awkward.ones_like(event.category)
    loweta=abs(event["muon_noniso_eta"])<1.6
    higheta=abs(event["muon_noniso_eta"])>1.6
    ptbin1=(event["muon_noniso_Tunept"]>50 )&( event["muon_noniso_Tunept"]<=100 )
    ptbin2=(event["muon_noniso_Tunept"]>100 )&( event["muon_noniso_Tunept"]<=150 )
    ptbin3=(event["muon_noniso_Tunept"]>150 )&( event["muon_noniso_Tunept"]<=200 )
    ptbin4=(event["muon_noniso_Tunept"]>200 )&( event["muon_noniso_Tunept"]<=300 )
    ptbin5=(event["muon_noniso_Tunept"]>300 )&( event["muon_noniso_Tunept"]<=400 )
    ptbin6=(event["muon_noniso_Tunept"]>400 )&( event["muon_noniso_Tunept"]<=600 )
    ptbin7=(event["muon_noniso_Tunept"]>600 )&( event["muon_noniso_Tunept"]<=1500 )
    ptbin8=(event["muon_noniso_Tunept"]>1500 )&( event["muon_noniso_Tunept"]<=3500)
    
    sf_loweta_ptbin1=((event.category==3)|(event.category==4)) & loweta & ptbin1
    sf_loweta_ptbin2=((event.category==3)|(event.category==4)) & loweta & ptbin2
    sf_loweta_ptbin3=((event.category==3)|(event.category==4)) & loweta & ptbin3
    sf_loweta_ptbin4=((event.category==3)|(event.category==4)) & loweta & ptbin4
    sf_loweta_ptbin5=((event.category==3)|(event.category==4)) & loweta & ptbin5
    sf_loweta_ptbin6=((event.category==3)|(event.category==4)) & loweta & ptbin6
    sf_loweta_ptbin7=((event.category==3)|(event.category==4)) & loweta & ptbin7
    sf_loweta_ptbin8=((event.category==3)|(event.category==4)) & loweta & ptbin8
    sf_higheta_ptbin1=((event.category==3)|(event.category==4)) & higheta & ptbin1
    sf_higheta_ptbin2=((event.category==3)|(event.category==4)) & higheta & ptbin2
    sf_higheta_ptbin3=((event.category==3)|(event.category==4)) & higheta & ptbin3
    sf_higheta_ptbin4=((event.category==3)|(event.category==4)) & higheta & ptbin4
    sf_higheta_ptbin5=((event.category==3)|(event.category==4)) & higheta & ptbin5
    sf_higheta_ptbin6=((event.category==3)|(event.category==4)) & higheta & ptbin6
    sf_higheta_ptbin7=((event.category==3)|(event.category==4)) & higheta & ptbin7
    sf_higheta_ptbin8=((event.category==3)|(event.category==4)) & higheta & ptbin8
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin1, awkward.ones_like(weight_muon_highptid_central)*0.9938, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin1, awkward.ones_like(weight_muon_highptid_up)*0.9946, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin1, awkward.ones_like(weight_muon_highptid_down)*0.9932, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin2, awkward.ones_like(weight_muon_highptid_central)*0.9950, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin2, awkward.ones_like(weight_muon_highptid_up)*0.9957, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin2, awkward.ones_like(weight_muon_highptid_down)*0.9943, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin3, awkward.ones_like(weight_muon_highptid_central)*0.996, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin3, awkward.ones_like(weight_muon_highptid_up)*0.997, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin3, awkward.ones_like(weight_muon_highptid_down)*0.995, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin4, awkward.ones_like(weight_muon_highptid_central)*0.996, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin4, awkward.ones_like(weight_muon_highptid_up)*0.997, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin4, awkward.ones_like(weight_muon_highptid_down)*0.995, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin5, awkward.ones_like(weight_muon_highptid_central)*0.994, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin5, awkward.ones_like(weight_muon_highptid_up)*0.995, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin5, awkward.ones_like(weight_muon_highptid_down)*0.993, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin6, awkward.ones_like(weight_muon_highptid_central)*1.003, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin6, awkward.ones_like(weight_muon_highptid_up)*1.009, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin6, awkward.ones_like(weight_muon_highptid_down)*0.997, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin7, awkward.ones_like(weight_muon_highptid_central)*0.987, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin7, awkward.ones_like(weight_muon_highptid_up)*0.990, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin7, awkward.ones_like(weight_muon_highptid_down)*0.984, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_loweta_ptbin8, awkward.ones_like(weight_muon_highptid_central)*0.9, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_loweta_ptbin8, awkward.ones_like(weight_muon_highptid_up)*1.0, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_loweta_ptbin8, awkward.ones_like(weight_muon_highptid_down)*0.8, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin1, awkward.ones_like(weight_muon_highptid_central)*1.0, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin1, awkward.ones_like(weight_muon_highptid_up)*1.0, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin1, awkward.ones_like(weight_muon_highptid_down)*1.0, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin2, awkward.ones_like(weight_muon_highptid_central)*0.993, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin2, awkward.ones_like(weight_muon_highptid_up)*0.994, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin2, awkward.ones_like(weight_muon_highptid_down)*0.992, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin3, awkward.ones_like(weight_muon_highptid_central)*0.989, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin3, awkward.ones_like(weight_muon_highptid_up)*0.990, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin3, awkward.ones_like(weight_muon_highptid_down)*0.988, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin4, awkward.ones_like(weight_muon_highptid_central)*0.986, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin4, awkward.ones_like(weight_muon_highptid_up)*0.987, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin4, awkward.ones_like(weight_muon_highptid_down)*0.985, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin5, awkward.ones_like(weight_muon_highptid_central)*0.989, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin5, awkward.ones_like(weight_muon_highptid_up)*0.990, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin5, awkward.ones_like(weight_muon_highptid_down)*0.988, weight_muon_highptid_down) 
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin6, awkward.ones_like(weight_muon_highptid_central)*0.983, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin6, awkward.ones_like(weight_muon_highptid_up)*0.986, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin6, awkward.ones_like(weight_muon_highptid_down)*0.980, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin7, awkward.ones_like(weight_muon_highptid_central)*0.986, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin7, awkward.ones_like(weight_muon_highptid_up)*0.992, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin7, awkward.ones_like(weight_muon_highptid_down)*0.980, weight_muon_highptid_down)
    
    weight_muon_highptid_central = awkward.where(sf_higheta_ptbin8, awkward.ones_like(weight_muon_highptid_central)*1.01, weight_muon_highptid_central)
    weight_muon_highptid_up = awkward.where(sf_higheta_ptbin8, awkward.ones_like(weight_muon_highptid_up)*1.02, weight_muon_highptid_up)
    weight_muon_highptid_down = awkward.where(sf_higheta_ptbin8, awkward.ones_like(weight_muon_highptid_down)*1.00, weight_muon_highptid_down)
    
    event["weight_muon_highptid_up"]=weight_muon_highptid_up
    event["weight_muon_highptid_central"]=weight_muon_highptid_central
    event["weight_muon_highptid_down"]=weight_muon_highptid_down
    return event


    

