# TODO items
#   1. Sync diphoton preselection with flashgg:
#       1.1 Calculate level of agreement using custom nanoAOD with necessary branches added.
#       1.2 Calculate level of agreement using central nanoAOD (which is missing some needed branches).
#       1.3 Decide how much of diphoton preselection should be configurable and how much should be hard-coded (both to improve readability and ensure consistency, i.e. not allowing users to change things they should not change) 
#   2. Decide how we want to handle events with more than one diphoton candidate.
#       - currently we keep all possible diphoton candidates in "awkward" format, meaning we store the diphotons as an events.Diphoton record, with potentially more than one per event
#       - we could either pick the "best" diphoton candidate per-event here (this makes life easier later down the line)
#       - OR we could allow taggers to have the option to consider multiple diphoton candidates per event but only pick one.
#           - if we do this, we need to have TagSequence (or each Tagger) come up with some criteria for picking the "best" diphoton candidate out of the selected ones
#           - do we want to allow multiple diphoton candidates in an event to be tagged by different taggers? E.g. TTHTagger selects the first diphoton candidate in the event and THQTagger selects the second diphoton candidate in the event?


import awkward
import time
import numpy
import numba
import vector

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG 
from higgs_dna.utils import awkward_utils

# NOTE: pre-compiled numba functions should be defined outside the class,
# as numba does not like when a class instance is passed to a function
@numba.njit
def produce_diphotons(photons, n_photons, lead_pt_cut):
    n_events = len(photons)

    diphotons_offsets = numpy.zeros(n_events + 1, numpy.int64)
    diphotons_contents = []
    lead_photons_idx = []
    sublead_photons_idx = []

    # Loop through events and select diphoton pairs
    for i in range(n_events):
        #diphotons_evt = []
        n_diphotons_event = 0
        # Enumerate photon pairs
        if n_photons[i] >= 2: # only try if there are at least 2 photons
            for j in range(n_photons[i]):
                for k in range(j, n_photons[i]):
                    # Choose who is the leading photon
                    lead_idx = j if photons[i][j].pt > photons[i][k].pt else k
                    sublead_idx = k if photons[i][j].pt > photons[i][k].pt else j
                    lead_photon = photons[i][lead_idx]
                    sublead_photon = photons[i][sublead_idx]

                    # Lead photon must satisfy lead pt cut
                    if lead_photon.pt < lead_pt_cut:
                        continue

                    # Construct four vectors
                    lead_photon_vec = vector.obj(
                            pt = lead_photon.pt,
                            eta = lead_photon.eta,
                            phi = lead_photon.phi,
                            mass = lead_photon.mass
                    )
                    sublead_photon_vec = vector.obj(
                            pt = sublead_photon.pt,
                            eta = sublead_photon.eta,
                            phi = sublead_photon.phi,
                            mass = sublead_photon.mass
                    )
                    diphoton = lead_photon_vec + sublead_photon_vec
                    
                    # Enforce diphoton mass requirements
                    if diphoton.mass < 100 or diphoton.mass > 180:
                        continue

                    # This diphoton candidate passes
                    n_diphotons_event += 1

                    # This doesn't work with numba
                    # so, we need to do each by hand
                    diphotons_contents.append([
                        diphoton.pt,
                        diphoton.eta,
                        diphoton.phi,
                        diphoton.mass
                    ])
                    
                    lead_photons_idx.append(lead_idx)
                    sublead_photons_idx.append(sublead_idx)

        diphotons_offsets[i+1] = diphotons_offsets[i] + n_diphotons_event 

    return diphotons_offsets, numpy.array(diphotons_contents), numpy.array(lead_photons_idx), numpy.array(sublead_photons_idx)

# NOTE: implementation still in progress.
# TODO: sync diphoton preselection with flashgg
class DiphotonTagger(Tagger):
    def calculate_selection(self, syst_tag, syst_events):
        """
        Select photons and create diphoton pairs.
        Add a record "Diphoton" to events array with relevant information about each diphoton pair.
        In principle, there can be more than one Diphoton pair per event.
        TODO: Should we select the "best" diphoton pair here (currently just storing this as an awkward array)? If so, how?
        """
        if "photon" not in syst_tag and not syst_tag == NOMINAL_TAG: # for non-applicable systs (e.g. leptons, jets), return the nominal selection and copy the nominal diphotons
            diphoton_selection, nominal_events = self.get_selection(NOMINAL_TAG, syst_events)
            syst_events = awkward.with_field( # make a copy so that subsequent taggers do not modify the same events twice
                    syst_events,
                    nominal_events.Diphoton,
                    "Diphoton"
            )

        else:
            photon_selection = self.select_photons(
                    photons = syst_events.Photon,
                    options = self.options["photons"]
            )

            diphoton_selection = self.produce_and_select_diphotons(
                    events = syst_events,
                    photons = syst_events.Photon[photon_selection],
                    options = self.options["diphotons"]
            )

        return diphoton_selection, syst_events


    def produce_and_select_diphotons(self, events, photons, options):
        """
        Perform diphoton preselection.
        For events with more than 2 photons, more than 1 diphoton candidate
        per event is possible.
        TODO: sync with flashgg

        :param events: events array to calculate diphoton candidates from
        :type events: awkward.highlevel.Array
        :param photons: array of selected photons from events
        :type photons: awkward.highlevel.Array
        :param options: dictionary containing configurable options for the diphoton preselection
        :type options: dict
        :return: boolean array indicating which events pass the diphoton preselection
        :rtype: awkward.highlevel.Array
        """

        start = time.time()

        # Produce diphotons
        n_photons = awkward.num(photons)
        lead_pt_cut = options["lead_pt"]
        diphoton_offsets, diphoton_contents, lead_photon_idx, sublead_photon_idx = produce_diphotons(photons, n_photons, lead_pt_cut) 

        # Zip diphotons into events array as records
        diphotons = {}
        diphotons["p4"] = awkward_utils.create_four_vectors(events, diphoton_offsets, diphoton_contents)
        #awkward_utils.create_record(events, diphotons, "Diphoton") 
        awkward_utils.add_field(events, "Diphoton", diphotons)

        lead_photon_idx = awkward_utils.construct_jagged_array(diphoton_offsets, lead_photon_idx)
        sublead_photon_idx = awkward_utils.construct_jagged_array(diphoton_offsets, sublead_photon_idx)

        # Zip lead/sublead photons into events.Diphoton as sub-records
        lead_photons = {}
        sublead_photons = {}
        for field in photons.fields:
            lead_photons[field] = photons[lead_photon_idx][field]
            sublead_photons[field] = photons[sublead_photon_idx][field]
        awkward_utils.add_field(events, ("Diphoton", "LeadPhoton"), lead_photons)
        awkward_utils.add_field(events, ("Diphoton", "SubleadPhoton"), sublead_photons)

        # Select events with at least 1 diphoton candidate
        dipho_presel_cut = awkward.num(events.Diphoton) >= 1
        self.register_cuts(
                names = "diphoton preselection",
                results = dipho_presel_cut
        )

        elapsed_time = time.time() - start
        logger.debug("[DiphotonTagger] %s, syst variation : %s, total time to execute select_diphotons: %.6f s" % (self.name, self.current_syst, elapsed_time))

        return dipho_presel_cut 

    def select_photons(self, photons, options):
        """
        Enforces all photon cuts that are commmon to both
        leading and subleading photons in the diphoton preselection.
        Cuts specific to a diphoton pair are not enforced here.

        :param photons: input collection of photons to use for calculating selected photons
        :type photons: awkward.highlevel.Array
        :param options: dictionary containing configurable options for the photon selection
        :type options: dict
        :return: boolean array indicating which photons pass the photon selection
        :rtype: awkward.highlevel.Array
        """
        # pt
        pt_cut = photons.pt > options["pt"]

        # eta
        eta_cut = Tagger.get_range_cut(abs(photons.eta), options["eta"])

        # electron veto
        e_veto_cut = photons.electronVeto > options["e_veto"]

        use_central_nano = options["use_central_nano"] # indicates whether we are using central nanoAOD (with some branches that are necessary for full diphoton preselection missing) or custom nanoAOD (with these branches added)

        # r9/isolation cut
        r9_cut = photons.r9 > options["r9"]

        if not use_central_nano:
            charged_iso_cut = photons.chargedHadronIso < options["charged_iso"]
            charged_rel_iso_cut = (photons.chargedHadronIso / photons.pt) < options["charged_rel_iso"]

        else:
            charged_iso_cut = (photons.pfRelIso03_chg * photons.pt) < options["charged_iso"]
            charged_rel_iso_cut = photons.pfRelIso03_chg < options["charged_rel_iso"]

        r9_iso_cut = r9_cut | charged_iso_cut | charged_rel_iso_cut

        # hadronic over em energy cut
        hoe_cut = photons.hoe < options["hoe"]

        # hlt-mimicking cuts
        photons_eb_high_r9 = photons.isScEtaEB & (photons.r9 > options["hlt"]["eb_high_r9"]["r9"]) 
        photons_eb_low_r9 = photons.isScEtaEB & (photons.r9 > options["hlt"]["eb_low_r9"]["r9"]) & (photons.r9 < options["hlt"]["eb_high_r9"]["r9"])

        photons_ee_high_r9 = photons.isScEtaEE & (photons.r9 > options["hlt"]["ee_high_r9"]["r9"])
        photons_ee_low_r9 = photons.isScEtaEE & (photons.r9 > options["hlt"]["ee_low_r9"]["r9"]) & (photons.r9 < options["hlt"]["ee_high_r9"]["r9"])

        if not use_central_nano:
            eb_low_r9_track_pt_cut = photons.trkSumPtHollowConeDR03 < options["hlt"]["eb_low_r9"]["track_sum_pt"]
            ee_low_r9_track_pt_cut = photons.trkSumPtHollowConeDR03 < options["hlt"]["ee_low_r9"]["track_sum_pt"]

        else:
            eb_low_r9_track_pt_cut = photons.pt > 0
            ee_low_r9_track_pt_cut = photons.pt > 0

        eb_low_r9_sigma_ieie_cut = photons.sieie < options["hlt"]["eb_low_r9"]["sigma_ieie"]
        ee_low_r9_sigma_ieie_cut = photons.sieie < options["hlt"]["eb_low_r9"]["sigma_ieie"]

        low_eta = abs(photons.eta) < options["hlt"]["eta_rho_corr"]

        if not use_central_nano:
            eb_low_r9_pho_iso_low_eta_cut = low_eta & (photons.photonIso * options["hlt"]["low_eta_rho_corr"] < options["hlt"]["eb_low_r9"]["pho_iso"]) 
            eb_low_r9_pho_iso_high_eta_cut = ~low_eta & (photons.photonIso * options["hlt"]["high_eta_rho_corr"] < options["hlt"]["eb_low_r9"]["pho_iso"]) 
        else:
            eb_low_r9_pho_iso_low_eta_cut = low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["low_eta_rho_corr"]) < options["hlt"]["eb_low_r9"]["pho_iso"])
            eb_low_r9_pho_iso_high_eta_cut = ~low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["high_eta_rho_corr"]) < options["hlt"]["eb_low_r9"]["pho_iso"])

        eb_low_r9_pho_iso_cut = eb_low_r9_pho_iso_low_eta_cut | eb_low_r9_pho_iso_high_eta_cut

        if not use_central_nano:
            ee_low_r9_pho_iso_low_eta_cut = low_eta & (photons.photonIso * options["hlt"]["low_eta_rho_corr"] < options["hlt"]["ee_low_r9"]["pho_iso"])
            ee_low_r9_pho_iso_high_eta_cut = ~low_eta & (photons.photonIso * options["hlt"]["high_eta_rho_corr"] < options["hlt"]["ee_low_r9"]["pho_iso"]) 
        else:
            ee_low_r9_pho_iso_low_eta_cut = low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["low_eta_rho_corr"]) < options["hlt"]["ee_low_r9"]["pho_iso"])
            ee_low_r9_pho_iso_high_eta_cut = ~low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["high_eta_rho_corr"]) < options["hlt"]["ee_low_r9"]["pho_iso"])
    
        ee_low_r9_pho_iso_cut = ee_low_r9_pho_iso_low_eta_cut | ee_low_r9_pho_iso_high_eta_cut

        hlt_cut = photons.pt < 0 # initialize to all False
        hlt_cut = hlt_cut | photons_eb_high_r9
        hlt_cut = hlt_cut | (photons_eb_low_r9 & eb_low_r9_track_pt_cut & eb_low_r9_sigma_ieie_cut & eb_low_r9_pho_iso_cut)
        hlt_cut = hlt_cut | photons_ee_high_r9
        hlt_cut = hlt_cut | (photons_ee_low_r9 & ee_low_r9_track_pt_cut & ee_low_r9_sigma_ieie_cut & ee_low_r9_pho_iso_cut)

        all_cuts = pt_cut & eta_cut & e_veto_cut & r9_iso_cut & hoe_cut & hlt_cut

        self.register_cuts(
                names = ["pt", "eta", "e_veto", "r9", "hoe", "hlt", "all"],
                results = [pt_cut, eta_cut, e_veto_cut, r9_iso_cut, hoe_cut, hlt_cut, all_cuts],
                cut_type = "photon"
        )

        return all_cuts 
