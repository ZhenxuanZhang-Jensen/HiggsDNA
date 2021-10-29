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

DEFAULT_OPTIONS = {
    "photons" : {
        "use_central_nano" : True,
        "pt" : 25.0,
        "eta" : [
            [0.0, 1.4442],
            [1.566, 2.5]
        ],
        "e_veto" : 0.5,
        "hoe" : 0.08,
        "r9" : 0.8,
        "charged_iso" : 20.0,
        "charged_rel_iso" : 0.3,
        "hlt" : {
            "eta_rho_corr" : 1.5,
            "low_eta_rho_corr" : 0.16544,
            "high_eta_rho_corr" : 0.13212,
            "eb_high_r9" : {
                "r9" : 0.85
            },
            "eb_low_r9" : {
                "r9" : 0.5,
                "pho_iso" : 4.0,
                "track_sum_pt" : 6.0,
                "sigma_ieie" : 0.015
            },
            "ee_high_r9" : {
                "r9" : 0.9
            },
            "ee_low_r9" : {
                "r9" : 0.8,
                "pho_iso" : 4.0,
                "track_sum_pt" : 6.0,
                "sigma_ieie" : 0.035
            }
        }
    },
    "diphotons" : {
        "lead_pt" : 35.0,
        "sublead_pt" : 25.0,
        "lead_pt_mgg" : 0.33,
        "sublead_pt_mgg" : 0.25,
        "select_highest_pt_sum" : True
    }

}

# NOTE: pre-compiled numba functions should be defined outside the class,
# as numba does not like when a class instance is passed to a function
@numba.njit
def produce_diphotons(photons, n_photons, lead_pt_cut, lead_pt_mgg_cut, sublead_pt_mgg_cut):
    n_events = len(photons)

    diphotons_offsets = numpy.zeros(n_events + 1, numpy.int64)
    diphotons_contents = []
    lead_photons_idx = []
    sublead_photons_idx = []

    # Loop through events and select diphoton pairs
    for i in range(n_events):
        n_diphotons_event = 0
        # Enumerate photon pairs
        if n_photons[i] >= 2: # only try if there are at least 2 photons
            sum_pt = 0
            for j in range(n_photons[i]):
                for k in range(j+1, n_photons[i]):
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
                    diphoton = vector.obj(px = 0., py = 0., pz = 0., E = 0.) # IMPORTANT NOTE: you need to initialize this to an empty vector first. Otherwise, you will get ZeroDivisionError exceptions for like 1 out of a million events (seemingly only with numba). 
                    diphoton = diphoton + lead_photon_vec + sublead_photon_vec

                    if (diphoton.mass < 100) | (diphoton.mass > 180):
                        continue

                    if lead_photon.pt / diphoton.mass < lead_pt_mgg_cut:
                        continue 

                    if sublead_photon.pt / diphoton.mass < sublead_pt_mgg_cut:
                        continue

                    # This diphoton candidate passes
                    n_diphotons_event += 1

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
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(DiphotonTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
    
    def calculate_selection(self, syst_tag, syst_events):
        """
        Select photons and create diphoton pairs.
        Add a record "Diphoton" to events array with relevant information about each diphoton pair.
        In principle, there can be more than one Diphoton pair per event.
        TODO: Should we select the "best" diphoton pair here (currently just storing this as an awkward array)? If so, how?
        """
        photon_selection = self.select_photons(
                photons = syst_events.Photon,
                rho = syst_events.fixedGridRhoAll if not self.options["photons"]["use_central_nano"] else awkward.ones_like(syst_events.Photon), # FIXME: to be deleted once fixedGridRhoAll is added to central nanoAOD
                options = self.options["photons"]
        )

        diphoton_selection, diphotons = self.produce_and_select_diphotons(
                events = syst_events,
                photons = syst_events.Photon[photon_selection],
                options = self.options["diphotons"]
        )

        return diphoton_selection, diphotons 


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
        diphoton_offsets, diphoton_contents, lead_photon_idx, sublead_photon_idx = produce_diphotons(photons, n_photons, options["lead_pt"], options["lead_pt_mgg"], options["sublead_pt_mgg"]) 

        # Zip diphotons into events array as records
        diphotons = {}
        diphotons["p4"] = awkward_utils.create_four_vectors(events, diphoton_offsets, diphoton_contents)

        lead_photon_idx = awkward_utils.construct_jagged_array(diphoton_offsets, lead_photon_idx)
        sublead_photon_idx = awkward_utils.construct_jagged_array(diphoton_offsets, sublead_photon_idx)

        # Zip lead/sublead photons into events.Diphoton as sub-records
        lead_photons = photons[lead_photon_idx]
        sublead_photons = photons[sublead_photon_idx]

        # Reshape events to have shape [n_events, n_diphotons_per_event]
        events = awkward.broadcast_arrays(events, lead_photon_idx)[0] # events now has shape [n_events, n_diphotons_per_event]

        # Add photon/diphoton fields
        awkward_utils.add_field(events, "Diphoton", diphotons)
        awkward_utils.add_field(events, ("Diphoton", "sumPt"), lead_photons.pt + sublead_photons.pt)
        awkward_utils.add_field(events, "LeadPhoton", lead_photons)
        awkward_utils.add_field(events, "SubleadPhoton", sublead_photons)

        # lead/sublead photons have shape [n_events, n_diphotons_per_event], but to merge them we need to give them shape [n_events, n_diphotons_per_event, 1]
        lead_photons = awkward.unflatten(
                lead_photons,
                counts = 1, # new dimension has length 1
                axis = -1 # add new dimension to innermost axis
        )
        sublead_photons = awkward.unflatten(
                sublead_photons,
                counts = 1,
                axis = -1
        )

        awkward_utils.add_field(
                events = events,
                name = ("Diphoton", "Photon"),
                data = awkward.concatenate(
                    [lead_photons, sublead_photons],
                    axis = 2 # keep constant along first two dimensions ([n_events, n_diphotons]), and merge along remaining dimensions
                )
        )

        n_events = len(events)

        # NOTE: the default option (following flashgg) is to pick the diphoton pair from each event with the highest sum_pt of the two photons.
        # For analyses with electrons in the final state, there can be a non-neglible fraction (O(5%)) of events with more than one diphoton pair (usually when the electron is reconstructed as a photon) and you may want to revisit this assumption.

        # Select highest pt sum diphoton
        if options["select_highest_pt_sum"]:
            logger.info("[DiphotonTagger : produce_and_select_diphotons] %s, selecting the highest pt_sum diphoton pair from each event." % (self.name))
            # Sort diphotons by sum_pt
            events = events[awkward.argsort(events.Diphoton.sumPt, ascending=False, axis=1)]
            # Take the highest sum_pt diphoton for each event
            diphotons = awkward.firsts(events) 
            # Remove events with no diphotons
            diphotons = diphotons[~awkward.is_none(diphotons)]

            n_diphotons_total = awkward.sum(awkward.num(events.Diphoton))
            n_diphotons_selected = len(diphotons)

            logger.debug("[DiphotonTagger : produce_and_select_diphotons] %s, syst variation : %s. Number of total diphoton candidates: %d, number of diphoton candidates after selecting candidate with highest pt sum in each event: %d (%.2f percent of diphoton events removed)." % (self.name, self.current_syst, n_diphotons_total, n_diphotons_selected, 100. * (float(n_diphotons_total - n_diphotons_selected) / float(n_diphotons_total))))
            

        # Otherwise flatten, so events has shape [n_diphotons]
        else:
            logger.info("[DiphotonTagger : produce_and_select_diphotons] %s, keeping all diphoton candidates in each event." % (self.name))
            n_events = len(events)
            diphotons = awkward.flatten(events)
            n_diphotons = len(diphotons)
            avg_diphotons = awkward.mean(awkward.num(events.Diphoton[awkward.num(events.Diphoton) >= 1]))
            logger.debug("[DiphotonTagger : produce_and_select_diphotons] %s, syst variation : %s. Number of events (before flattening): %d. Number of diphotons (after flattening): %d (%.4f diphotons per event for events with at least 1 diphoton)." % (self.name, self.current_syst, n_events, n_diphotons, avg_diphotons)) 

        dipho_presel_cut = (diphotons.Diphoton.mass > 100.) & (diphotons.Diphoton.mass < 180.) 
        self.register_cuts(
                names = "diphoton preselection",
                results = dipho_presel_cut
        )

        elapsed_time = time.time() - start
        logger.debug("[DiphotonTagger] %s, syst variation : %s, total time to execute select_diphotons: %.6f s" % (self.name, self.current_syst, elapsed_time))

        return dipho_presel_cut, diphotons 


    def select_photons(self, photons, rho, options):
        """
        Enforces all photon cuts that are commmon to both
        leading and subleading photons in the diphoton preselection.
        Cuts specific to a diphoton pair are not enforced here.

        :param photons: input collection of photons to use for calculating selected photons
        :type photons: awkward.highlevel.Array
        :param rho: energy density in each event, used for corrections to photon isolation cuts
        :type rho: awkward.highlevel.Array
        :param options: dictionary containing configurable options for the photon selection
        :type options: dict
        :return: boolean array indicating which photons pass the photon selection
        :rtype: awkward.highlevel.Array
        """
        # pt
        pt_cut = photons.pt > options["pt"]

        # eta
        #eta_cut = Tagger.get_range_cut(abs(photons.eta), options["eta"]) | (photons.isScEtaEB | photons.isScEtaEE)
        eta_cut = (photons.isScEtaEB | photons.isScEtaEE)

        # electron veto
        e_veto_cut = photons.electronVeto > options["e_veto"]

        use_central_nano = options["use_central_nano"] # indicates whether we are using central nanoAOD (with some branches that are necessary for full diphoton preselection missing) or custom nanoAOD (with these branches added)

        # r9/isolation cut
        r9_cut = photons.r9 > options["r9"]

        if not use_central_nano:
            charged_iso_cut = photons.chargedHadronIso < options["charged_iso"]
            charged_rel_iso_cut = (photons.chargedHadronIso / photons.pt) < options["charged_rel_iso"]

        else: # FIXME: to be deleted once Photon.chargedHadronIso is added to central nanoAOD
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

        else: # FIXME: to be deleted once Photon.trkSumPtHollowConeDR03 is added to central nanoAOD
            eb_low_r9_track_pt_cut = photons.pt > 0
            ee_low_r9_track_pt_cut = photons.pt > 0

        eb_low_r9_sigma_ieie_cut = photons.sieie < options["hlt"]["eb_low_r9"]["sigma_ieie"]
        ee_low_r9_sigma_ieie_cut = photons.sieie < options["hlt"]["ee_low_r9"]["sigma_ieie"]

        low_eta = abs(photons.eta) < options["hlt"]["eta_rho_corr"]

        if not use_central_nano:
            eb_low_r9_pho_iso_low_eta_cut = low_eta & (photons.photonIso - rho * options["hlt"]["low_eta_rho_corr"] < options["hlt"]["eb_low_r9"]["pho_iso"]) 
            eb_low_r9_pho_iso_high_eta_cut = ~low_eta & (photons.photonIso - rho * options["hlt"]["high_eta_rho_corr"] < options["hlt"]["eb_low_r9"]["pho_iso"]) 
        else: # FIXME: to be deleted once Photon.photonIso is added to central nanoAOD
            eb_low_r9_pho_iso_low_eta_cut = low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["low_eta_rho_corr"]) < options["hlt"]["eb_low_r9"]["pho_iso"])
            eb_low_r9_pho_iso_high_eta_cut = ~low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["high_eta_rho_corr"]) < options["hlt"]["eb_low_r9"]["pho_iso"])

        eb_low_r9_pho_iso_cut = eb_low_r9_pho_iso_low_eta_cut | eb_low_r9_pho_iso_high_eta_cut

        if not use_central_nano:
            ee_low_r9_pho_iso_low_eta_cut = low_eta & (photons.photonIso - rho * options["hlt"]["low_eta_rho_corr"] < options["hlt"]["ee_low_r9"]["pho_iso"]) 
            ee_low_r9_pho_iso_high_eta_cut = ~low_eta & (photons.photonIso - rho * options["hlt"]["high_eta_rho_corr"] < options["hlt"]["ee_low_r9"]["pho_iso"]) 
        else: # FIXME: to be deleted once Photon.photonIso is added to central nanoAOD
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
