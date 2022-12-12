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

ELECTRON_ID_SF = {
    "2016" : "2016postVFP",
    "2016UL_preVFP" : "2016preVFP",
    "2016UL_postVFP" : "2016postVFP",
    "2017" : "2017",
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
