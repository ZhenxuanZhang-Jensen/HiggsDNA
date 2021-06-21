# TODO items
# 1. Make the selection more configurable
# 2. Should this (and selections more generally) be implemented as a class, where a Tagger will own multiple Selection instances?
#   - this would allow for TagSequence to more efficiently manage computation and avoid repeated computation (e.g. multiple taggers using the same jet selection)

import awkward

from higgs_dna.selections import object_selections

def select_jets(jets, diphotons, muons, electrons, options, tagger = None):
    pt_cut = jets.pt > 25
    eta_cut = abs(jets.eta) < 2.4

    dr_lead_pho_cut = object_selections.delta_R(jets, diphotons.LeadPhoton, 0.4)
    dr_sublead_pho_cut = object_selections.delta_R(jets, diphotons.SubleadPhoton, 0.4)

    dr_electrons_cut = object_selections.delta_R(jets, electrons, 0.4)
    dr_muons_cut = object_selections.delta_R(jets, muons, 0.4)

    jet_cut = pt_cut & eta_cut & dr_lead_pho_cut & dr_sublead_pho_cut & dr_electrons_cut & dr_muons_cut

    if tagger is not None:
        tagger.register_cuts(
                names = ["pt", "eta", "dr_lead_photon", "dr_sublead_photon", "dr_electrons", "dr_muons", "all"],
                results = [pt_cut, eta_cut, dr_lead_pho_cut, dr_sublead_pho_cut, dr_electrons_cut, dr_muons_cut, jet_cut],
                cut_type = "jet"
        ) 

    return jet_cut
