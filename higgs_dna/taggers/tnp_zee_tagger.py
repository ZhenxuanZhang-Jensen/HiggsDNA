import awkward
import time
import numpy
import numba
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG 
from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.selections import gen_selections, object_selections

DEFAULT_OPTIONS = {
    "photons" : {
        "use_central_nano" : True,
        "pt" : 25.0,
        "eta" : [
            [0.0, 1.4442],
            [1.566, 2.5]
        ],
        "e_veto" : 0.5,
        "pixel_seed" : 0.5,
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
    "zee_cands" : {
        "mass" : [50., 200.]
    },
    "tag" : {
        "pt" : 40.0,
        "id" : "WP80",
        "trig_obj_dr" : 0.1
    },
    "trigger" : {
        "2016" : ["HLT_Ele27_WPTight_Gsf"],
        "2016UL_preVFP" : ["HLT_Ele27_WPTight_Gsf"],
        "2016UL_postVFP" : ["HLT_Ele27_WPTight_Gsf"],
        "2017" : ["HLT_Ele35_WPTight_Gsf"],
        "2018" : ["HLT_Ele32_WPTight_Gsf"],
    },
    "gen_info" : {
        "calculate" : False,
        "max_dr" : 0.2,
        "max_pt_diff" : 15.
    }  
}


# Diphoton preselection below synced with flashgg, see details in:
#   - https://indico.cern.ch/event/1071721/contributions/4551056/attachments/2320292/3950844/HiggsDNA_DiphotonPreselectionAndSystematics_30Sep2021.pdf

class TnPZeeTagger(Tagger):
    def __init__(self, name = "default_diphoton_tagger", options = {}, is_data = None, year = None):
        super(TnPZeeTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, events):
        """
        Select photons and create diphoton pairs.
        Add a record "Diphoton" to events array with relevant information about each diphoton pair.
        In principle, there can be more than one Diphoton pair per event.
        """
        
        # Determine what type of rho variable is available in nanoAOD
        # To be deleted once a standard rho variable is added to central nanoAOD
        if not self.options["photons"]["use_central_nano"]:
            if "fixedGridRhoAll" in events.fields:
                rho = events.fixedGridRhoAll
            elif "fixedGridRhoFastjetAll" in events.fields:
                rho = events.fixedGridRhoFastjetAll
            elif "Rho_fixedGridRhoAll" in events.fields:
                rho = events.Rho_fixedGridRhoAll
            else:
                logger.exception("[DiphotonTagger : calculate_selection] Did not find valid 'rho' field.")
                raise RuntimeError()
        else:
            rho = awkward.ones_like(events.Photon)

        photon_selection = self.select_photons(
                photons = events.Photon,
                rho = rho,
                options = self.options["photons"]
        )

        photons = events.Photon[photon_selection]
        photons = photons[awkward.argsort(photons.pt, ascending=False, axis=1)]
        photons = awkward.Array(photons, with_name = "Momentum4D")
    
        # Tag definition
        ele_cut = photons.electronIdx >= 0
        pt_cut = photons.pt > self.options["tag"]["pt"]

        # Grab electron ID and assign it to photons
        for field in events.Electron.fields:
            if field in photons.fields: # don't overwrite existing fields like pt/eta/phi
                continue

            electron_idx_safe = awkward.where(photons.electronIdx >= 0, photons.electronIdx, awkward.zeros_like(photons.electronIdx))
            electrons_padded = awkward.pad_none(events.Electron, 1)
            photons[field] = awkward.where(
                    photons.electronIdx >= 0,
                    electrons_padded[electron_idx_safe][field],
                    awkward.ones_like(photons.pt) * -999.
            )
                    
        if self.options["tag"]["id"] == "WP80":
            id_cut = photons.mvaFall17V2Iso_WP80 == True
        elif self.options["tag"]["id"] == "WP90":
            id_cut = photons.mvaFall17V2Iso_WP90 == True
        else:
            id_cut = awkward.ones_like(photons, dtype=bool) 

        # Require tags to be dR matched to a TrigObj which fired the trigger
        electron_trig_obj = events.TrigObj[events.TrigObj.id == 11] # require the trig obj to be an electron
        electron_trig_obj["bitlist"] = awkward_utils.to_bitlist(electron_trig_obj.filterBits, 18)
        electron_trig_obj = electron_trig_obj[electron_trig_obj.bitlist[:,:,1] == True] # require the trig obj to be matched to the hltEle*WPTight*TrackIsoFilter* bit (bit 2)
        electron_trig_obj["mass"] = awkward.zeros_like(electron_trig_obj.pt)
        
        trig_obj_dr_cut = object_selections.delta_R(photons, electron_trig_obj, self.options["tag"]["trig_obj_dr"], mode = "max") 

        photons["pass_tag"] = ele_cut & pt_cut & id_cut & trig_obj_dr_cut
       
        self.register_cuts(
                names = ["in electron collection", "pt", "ID", "dR matched to trigger object", "all"],
                results = [ele_cut, pt_cut, id_cut, trig_obj_dr_cut, photons.pass_tag],
                cut_type = "tag"
        )

        photons = photons[awkward.argsort(photons.pass_tag * photons.pt, ascending=False, axis=1)] # photons that pass tag selection come first, then sorted by pT

        zee_pairs = awkward.combinations(photons, 2, fields=["TagPhoton", "ProbePhoton"]) # we will require at least one photon passing the tag and mark the highest pT photon as the tag in events with more than one
        zee_pairs["Diphoton"] = zee_pairs.TagPhoton + zee_pairs.ProbePhoton

        # Z candidate cuts
        n_tags = awkward.num(photons[photons.pass_tag])
        tag_cut = n_tags >= 1 # at least 1 tag

        n_photons = awkward.num(photons)
        n_pho_cut = n_photons == 2 # exactly 2 photons

        # mass on Z peak
        mass_cut = (zee_pairs.Diphoton.mass >= self.options["zee_cands"]["mass"][0]) & (zee_pairs.Diphoton.mass <= self.options["zee_cands"]["mass"][1]) 

        # trigger
        if self.year is not None:
            if len(self.options["trigger"][self.year]) == 0: # no triggers:
                trigger_cut = awkward.ones_like(events.run, dtype=bool) # dummy cut, all True
            else:
                trigger_cut = awkward.zeros_like(events.run, dtype=bool) # dummy cut, all False
                for hlt in self.options["trigger"][self.year]: # logical OR of all triggers
                    trigger_cut = (trigger_cut) | (events[hlt] == True)
        else:
            trigger_cut = awkward.ones_like(events.run, dtype=bool) # dummy cut, all True

        zee_cut = tag_cut & n_pho_cut & mass_cut & trigger_cut

        zee_pairs[("Diphoton", "mass")] = zee_pairs.Diphoton.mass
        zee_pairs[("Diphoton", "eta")] = zee_pairs.Diphoton.eta
        zee_pairs[("Diphoton", "phi")] = zee_pairs.Diphoton.phi
        zee_pairs[("Diphoton", "pt")] = zee_pairs.Diphoton.pt

        # Event level cut
        presel_cut = awkward.num(zee_pairs[zee_cut]) == 1

        zee_pairs = awkward.firsts(zee_pairs[zee_cut])

        for field in ["Diphoton", "TagPhoton", "ProbePhoton"]:
            events[field] = zee_pairs[field]

        self.register_cuts(
                names = ["passing tag", "2 photons", "m_ee", "trigger", "all"],
                results = [tag_cut, n_pho_cut, mass_cut, trigger_cut, presel_cut],
                cut_type = "event"
        ) 

        return presel_cut, events 


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
        e_veto_cut = photons.electronVeto < options["e_veto"]
        psv_cut = photons.pixelSeed < options["pixel_seed"] 

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

        # Check which type of photonIso name is present in these nanoAODs
        if not use_central_nano: # FIXME: to be deleted once Photon.photonIso is added to central nanoAOD
            if "pfPhoIso03" in photons.fields:
                photon_iso = photons.pfPhoIso03
            elif "photonIso" in photons.fields:
                photon_iso = photons.photonIso
            else:
                logger.warning("[TnPZeeTagger : select_photons] Did not find an appropriate photon isolation variable in events array.")

        if not use_central_nano:
            
            eb_low_r9_pho_iso_low_eta_cut = low_eta & (photon_iso - rho * options["hlt"]["low_eta_rho_corr"] < options["hlt"]["eb_low_r9"]["pho_iso"]) 
            eb_low_r9_pho_iso_high_eta_cut = ~low_eta & (photon_iso - rho * options["hlt"]["high_eta_rho_corr"] < options["hlt"]["eb_low_r9"]["pho_iso"]) 
        else: # FIXME: to be deleted once Photon.photonIso is added to central nanoAOD
            eb_low_r9_pho_iso_low_eta_cut = low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["low_eta_rho_corr"]) < options["hlt"]["eb_low_r9"]["pho_iso"])
            eb_low_r9_pho_iso_high_eta_cut = ~low_eta & ((photons.pfRelIso03_all * photons.pt * options["hlt"]["high_eta_rho_corr"]) < options["hlt"]["eb_low_r9"]["pho_iso"])

        eb_low_r9_pho_iso_cut = eb_low_r9_pho_iso_low_eta_cut | eb_low_r9_pho_iso_high_eta_cut

        if not use_central_nano:
            ee_low_r9_pho_iso_low_eta_cut = low_eta & (photon_iso - rho * options["hlt"]["low_eta_rho_corr"] < options["hlt"]["ee_low_r9"]["pho_iso"]) 
            ee_low_r9_pho_iso_high_eta_cut = ~low_eta & (photon_iso - rho * options["hlt"]["high_eta_rho_corr"] < options["hlt"]["ee_low_r9"]["pho_iso"]) 
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

