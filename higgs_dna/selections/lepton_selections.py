# TODO items
# 1. Make the selection more configurable
# 2. Should this (and selections more generally) be implemented as a class, where a Tagger will own multiple Selection instances?
#   - this would allow for TagSequence to more efficiently manage computation and avoid repeated computation (e.g. multiple taggers using the same jet selection)

import awkward

from higgs_dna.selections import object_selections

def select_electrons(electrons, diphotons, options, tagger = None):
    """
    TODO
    """
    pt_cut = electrons.pt > options["pt"]
    eta_cut = abs(electrons.eta) < 2.4
    id_cut = electrons.mvaFall17V2Iso_WP90 == True

    dr_pho_cut = object_selections.delta_R(electrons, diphotons.Photon, 0.2)

    electron_cut = pt_cut & eta_cut & id_cut & dr_pho_cut

    # For diagnostic info
    if tagger is not None:
        tagger.register_cuts(
                names = ["pt", "eta", "id", "dr_photons", "all"],
                results = [pt_cut, eta_cut, id_cut, dr_pho_cut, electron_cut],
                cut_type = "electron"
        )

    return electron_cut


def select_muons(muons, diphotons, options, tagger = None):
    """
    TODO
    """
    pt_cut = muons.pt > 25
    eta_cut = abs(muons.eta) < 2.4
    id_cut = muons.mediumId == True

    dr_pho_cut = object_selections.delta_R(muons, diphotons.Photon, 0.2)

    muon_cut = pt_cut & eta_cut & id_cut & dr_pho_cut

    if tagger is not None:
        tagger.register_cuts(
                names = ["pt", "eta", "id", "dr_photons", "all"],
                results = [pt_cut, eta_cut, id_cut, dr_pho_cut, muon_cut],
                cut_type = "muon"
        )

    return muon_cut


