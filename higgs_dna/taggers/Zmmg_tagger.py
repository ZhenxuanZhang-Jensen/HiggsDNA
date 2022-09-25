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
#from higgs_dna.selections import gen_selections

DEFAULT_OPTIONS = {
    "muons": {
        "pt" : 10.0,
        "eta" : 2.5,
        "id" : "tight",
        "pfRelIso03_chg" : 0.2,
        "global" : True
    },
    "photons" : {
        "use_central_nano" : True,
        "pt" : 20.0,
        "eta" : [
            [0.0, 1.4442],
            [1.566, 2.5]
        ],
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
    "dimuons" : {
        "lead_pt" : 20.0,
        "oppositeCharge" : True,
        "mass" : [35., 180.],
        "select_highest_pt_sum" : False
    },
    "Zmmg" : {
        "farmu_pt" : 20.0,
        "mass" : [60., 120.],
        "min_dRmg" : 0.8,
        "mm_mmg_summass" : [60., 180.],
        "select_closest_PDGmass" : True
    },
    "trigger" : {
        "2016" : ["HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"],
        "2017" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2018" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"]
    } 
}


# Selections of Zmmg Tagger:

class ZmmgTagger(Tagger):
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(ZmmgTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, syst_tag, syst_events):
        """
        Select muons and photons, and create dimuon pairs and mmg objects.
        Add a record "Zmmg" to events array with relevant information about each Zmmg candidate.
        In principle, there can be more than one mmg pair per event.
        """
        muon_selection   = self.select_muons(
                muons = syst_events.Muon,
                options = self.options["muons"]
        ) 
        photon_selection = self.select_photons(
                photons = syst_events.Photon,
                rho = syst_events.fixedGridRhoAll if not self.options["photons"]["use_central_nano"] else awkward.ones_like(syst_events.Photon), # FIXME: to be deleted once fixedGridRhoAll is added to central nanoAOD
                options = self.options["photons"]
        )
        #dimu
        dimuon_selection, dimuons = self.produce_and_select_dimuons(
                events = syst_events,
                muons = syst_events.Muon[muon_selection],
                options = self.options["dimuons"]
        )
        #mmg
        Zmmg_selection, Zmmgs = self.produce_and_select_Zmmg(
                events = syst_events,
                dimuons,
                photons = syst_events.Photon[photon_selection],
                options = self.options["Zmmg"]
        )

        return Zmmg_selection, Zmmgs 


    def produce_and_select_dimuons(self, events, muons, options):
        start = time.time()
        # Sort muons by pt
        muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        muons = awkward.Array(muons, with_name = "Momentum4D")

        # Get all combinations of two muons in each event
        muons = awkward.combinations(muons, 2, fields=["LeadMuon", "SubleadMuon"])
        dimuons["DiMuon"] = dimuons.LeadMuon + dimuons.SubleadMuon

        # Add sumPt and dR for convenience
        dimuons[("DiMuon", "sumPt")] = dimuons.LeadMuon.pt + dimuons.SubleadMuon.pt
        dimuons[("DiMuon", "dR")] = dimuons.LeadMuon.deltaR(dimuons.SubleadMuon)        

        # Add lead/sublead muons to additionally be accessible together as dimuons.DiMuon.Muon
        lead_muons = awkward.unflatten(
                dimuons.LeadMuon,
                counts = 1, # new dimension has length 1
                axis = -1 # add new dimension to innermost axis
        )

        sublead_muons = awkward.unflatten(
                dimuons.SubleadMuon,
                counts = 1,
                axis = -1
        )

        awkward_utils.add_field(
                events = dimuons,
                name = ("DiMuon", "Muon"),
                data = awkward.concatenate(
                    [lead_muons, sublead_muons],
                    axis = 2 # keep constant along first two dimensions ([n_events, n_dimuons]), and merge along remaining dimensions
                )
        )

        lead_pt_cut = dimuons.LeadMuon.pt >= options["lead_pt"]
        if options["oppositeCharge"]:
           charge_cut = (dimuons.LeadMuon.charge * dimuons.SubleadMuon.charge) == -1
        else:
           charge_cut = True # dummy cut, all True
        mass_cut = (dimuons.DiMuon.mass >= options["mass"][0]) & (dimuons.DiMuon.mass <= options["mass"][1])
        all_cuts = lead_pt_cut & charge_cut & mass_cut

        self.register_cuts(
            names = ["lead pt cut",  "charge cut", "mass cut", "all cuts"],
            results = [lead_pt_cut, charge_cut, mass_cut, all_cuts],
            cut_type = "dimuon"
        )

        dimuons = dimuons[all_cuts]

        # Sort by sumPt
        dimuons = dimuons[awkward.argsort(dimuons.DiMuon.sumPt, ascending=False, axis=1)]
        if options["select_highest_pt_sum"]:
            logger.debug("[ZmmgTagger : produce_and_select_dimuons] %s, selecting the highest pt_sum dimuon pair from each event." % (self.name))
            n_dimuons_total = awkward.sum(awkward.num(dimuons))
            dimuons = awkward.singletons(awkward.firsts(dimuons)) # ak.firsts takes the first dimuon in each event (highest sumpt) and ak.singletons makes it variable length again
            n_dimuons_selected = awkward.sum(awkward.num(dimuons)) 
            logger.debug("[ZmmgTagger : produce_and_select_dimuons] %s, syst variation : %s. Number of total dimuon candidates: %d, number of dimuon candidates after selecting candidate with highest pt sum in each event: %d (%.2f percent of diphoton events removed)." % (self.name, self.current_syst, n_dimuons_total, n_dimuons_selected, 100. * (float(n_dimuons_total - n_dimuons_selected) / float(n_dimuons_total))))
        # Otherwise, keep all diphotons
        else:
            logger.debug("[ZmmgTagger : produce_and_select_dimuons] %s, keeping all dimuon candidates in each event." % (self.name))
            avg_dimuons = awkward.mean(awkward.num(dimuons))
            logger.debug("[ZmmgTagger : produce_and_select_dimuons] %s, syst variation : %s, %.4f dimuons per event for events with at least 1 dimuon)." % (self.name, self.current_syst, avg_dimuons)) 
        # Reshape events to have shape [n_events, n_dimuons_per_event] and then flatten.
        # This is necessary so event-level variables are properly copied for each dimuon per event
        dimu_events = awkward.broadcast_arrays(events, dimuons.DiMuon.mass)[0]

        # Add to events
        for field in ["DiMuon", "LeadMuon", "SubleadMuon"]:
            dimu_events[field] = dimuons[field] 
            dimu_events[(field, "pt")] = dimuons[field].pt
            dimu_events[(field, "eta")] = dimuons[field].eta
            dimu_events[(field, "phi")] = dimuons[field].phi
            dimu_events[(field, "mass")] = dimuons[field].mass

        dimu_presel_cut = awkward.num(dimu_events.DiMuon) >= 1
        if self.is_data and self.year is not None:
            trigger_cut = awkward.num(dimu_events.DiMuon) < 0 # dummy cut, all False
            for hlt in self.options["trigger"][self.year]: # logical OR of all triggers
                trigger_cut = (trigger_cut) | (dimu_events[hlt] == True)
        else:
            trigger_cut = awkward.num(dimu_events.DiMuon) >= 0 # dummy cut, all True

        presel_cut = dimu_presel_cut & trigger_cut

        self.register_cuts(
            names = ["At least 1 dimuon pair", "HLT Trigger", "all"],
            results = [dimu_presel_cut, trigger_cut, presel_cut]
        )

#        dimu_events = dimu_events[dimu_presel_cut]
        dimu_events = dimu_events[presel_cut]

        dimu_events = awkward.flatten(dimu_events)

        elapsed_time = time.time() - start
        logger.debug("[ZmmgTagger] %s, syst variation : %s, total time to execute select_dimuons: %.6f s" % (self.name, self.current_syst, elapsed_time))

        dummy_cut = dimu_events.DiMuon.pt > 0
        return dummy_cut, dimu_events 


    def produce_and_select_Zmmg(self, events, dimuons, photons, options):
        """
        Perform Zmmg selections.
        
        """

        start = time.time()

        Zmmgs["ZmmgDiMu"] =  dimuons
        Zmmgs["ZmmgPho"] =  photons  
        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        photons = awkward.Array(photons, with_name = "Momentum4D")
        dimuons = awkward.Array(dimuons, with_name = "Momentum4D")

        # Get all combinations of dimuon +  photon in each event        
        Zmmgs["Zmmg"] = photons + dimuons
        ## 
        lead_dR = dimuons.LeadMuon.deltaR(photons) 
        sublead_dR = dimuons.SubleadMuon.deltaR(photons)  
        min_DR = sublead_dR 
        FarMuon = dimuons.LeadMuon
        NearMuon =  dimuons.SubleadMuon
        if lead_dR < sublead_dR:
           min_DR = lead_dR
           FarMuon = dimuons.SubleadMuon
           NearMuon = dimuons.LeadMuon

        Zmmgs[("Zmmg", "min_dR_mg")] = min_DR
        Zmmgs[("Zmmg", "DMassPDG")] = DeltaMass2PDG
        Zmmgs["ZmmgFarMu"] = FarMuon 
        Zmmgs["ZmmgNearMu"] = NearMuon  

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
    
        # Add 
        Zmmg_photons = awkward.unflatten(
                photons,
                counts = 1, # new dimension has length 1
                axis = -1 # add new dimension to innermost axis
        )
        Zmmg_dimuons = awkward.unflatten(
                dimuons,
                counts = 1,
                axis = -1
        )
        Zmmg_farmuons = awkward.unflatten(
                FarMuon,
                counts = 1,
                axis = -1
        )
        Zmmg_nearmuons = awkward.unflatten(
                NearMuon,
                counts = 1,
                axis = -1
        )

        awkward_utils.add_field(
                events = Zmmgs,
                name = ("ZmmgDiMu", "ZmmgPho", "ZmmgFarMu", "ZmmgNearMu"),
                data = awkward.concatenate(
                    [Zmmg_dimuons, Zmmg_photons, Zmmg_farmuons, Zmmg_nearmuons],
                    axis = 4 # keep constant along first two dimensions ([n_events, n_diphotons]), and merge along remaining dimensions
                )
        )

        minDR_cut  =  min_DR  < options["min_dRmg"]                            
        farmuPT_cut = FarMuon.pt  >= options["farmu_pt"]
        mass_cut = (Zmmgs.Zmmg.mass >= options["mass"][0]) & (diphotons.Diphoton.mass <= options["mass"][1])
        mm_plus_mmg_mass = Zmmgs.Zmmg.mass + dimuons.DiMuon.mass
        sum_mass_cut = (mm_plus_mmg_mass >= options["mm_mmg_summass"][0]) & (mm_plus_mmg_mass <= options["mm_mmg_summass"][1])

        all_cuts = minDR_cut & farmuPT_cut & mass_cut & sum_mass_cut

        self.register_cuts(
            names = ["min DRmg cut", "farmu pt cut", "mmg mass cut", "mmgAmm mass cut", "all cuts"],
            results = [minDR_cut, farmuPT_cut, mass_cut, sum_mass_cut, all_cuts],
            cut_type = "Zmmg"
        )

        Zmmgs = Zmmgs[all_cuts]

        # Sort by sumPt
        DeltaMass2PDG = abs(Zmmgs.Zmmg.mass - 91.188)
        Zmmgs = Zmmgs[awkward.argsort(Zmmgs.Zmmg.DMassPDG, ascending=True, axis=1)]

        # Select closest mass to PDG value
        if options["select_closest_PDGmass"]:
            logger.debug("[ZmmgTagger : produce_and_select_Zmmgs] %s, selecting the Zmmg with closest PDG mass value from each event." % (self.name))

            n_Zmmgs_total = awkward.sum(awkward.num(Zmmgs))
            Zmmgs = awkward.singletons(awkward.firsts(Zmmgs)) # ak.firsts takes the first object in each event and ak.singletons makes it variable length again
            n_Zmmgs_selected = awkward.sum(awkward.num(Zmmgs)) 

            logger.debug("[ZmmgTagger : produce_and_select_Zmmgs] %s, syst variation : %s. Number of total Zmmg candidates: %d, number of Zmmg candidates after selecting candidate with closest PDG mass value in each event: %d (%.2f percent of Zmmg events removed)." % (self.name, self.current_syst, n_Zmmgs_total, n_Zmmgs_selected, 100. * (float(n_Zmmgs_total - n_Zmmgs_selected) / float(n_Zmmgs_total))))

        # Otherwise, keep all tags
        else:
            logger.debug("[ZmmgTagger : produce_and_select_Zmmgs] %s, keeping all Zmmgs candidates in each event." % (self.name))
            avg_Zmmgs = awkward.mean(awkward.num(Zmmgs))
            logger.debug("[ZmmgTagger : produce_and_select_Zmmgs] %s, syst variation : %s, %.4f Zmmgs per event for events with at least 1 Zmmg)." % (self.name, self.current_syst, avg_Zmmgs)) 

        # Reshape events to have shape [n_events, n_Zmmgs_per_event] and then flatten.
        # This is necessary so event-level variables are properly copied for each Zmmg per event
        Zmmg_events = awkward.broadcast_arrays(events, Zmmgs.Zmmg.mass)[0]

        # Add to events
        for field in ["Zmmg", "ZmmgDiMu", "ZmmgPho", "ZmmgFarMu", "ZmmgNearMu"]:
            Zmmg_events[field] = Zmmgs[field] 
            Zmmg_events[(field, "pt")] = Zmmgs[field].pt
            Zmmg_events[(field, "eta")] = Zmmgs[field].eta
            Zmmg_events[(field, "phi")] = Zmmgs[field].phi
            Zmmg_events[(field, "mass")] = Zmmgs[field].mass

        Zmmg_presel_cut = awkward.num(Zmmg_events.Zmmg) >= 1
 
        Zmmg_events = Zmmg_events[Zmmg_presel_cut]

        Zmmg_events = awkward.flatten(Zmmg_events)

        elapsed_time = time.time() - start
        logger.debug("[ZmmgTagger] %s, syst variation : %s, total time to execute select_Zmmg: %.6f s" % (self.name, self.current_syst, elapsed_time))

        dummy_cut = Zmmg_events.Zmmg.pt > 0
        return dummy_cut, Zmmg_events 

    
    
        
    def select_muons(self, muons, options):
         """
        Enforces all muon cuts          
        """
        # pt
        pt_cut = muons.pt > options["pt"]
        eta_cut = abs(muons.eta) < options["eta"]
 
        if options["id"] == "medium":
            id_cut = muons.mediumId == True
        elif options["id"] == "tight":
            id_cut = muons.tightId == True
        elif not options["id"] or options["id"].lower() == "none":
            id_cut = muons.pt > 0.
        else:
            logger.warning("[select_muons] : Tagger '%s', id cut '%s' not recognized, not applying an ID cut." % (str(tagger), options["id"]))
            id_cut = muons.pt > 0.
        id_cut = id_cut  & (muons.pfRelIso03_chg < options["pfRelIso03_chg"])

        if options["global"]:
            global_cut = muons.isGlobal == True
        else:
            global_cut = muons.pt > 0

        all_cuts = pt_cut & eta_cut & id_cut & global_cut
        return all_cuts

    def select_photons(self, photons, rho, options):
        """
        Enforces all photon cuts 
        """
        # pt
        pt_cut = photons.pt > options["pt"]
    
        # eta
        #eta_cut = Tagger.get_range_cut(abs(photons.eta), options["eta"]) | (photons.isScEtaEB | photons.isScEtaEE)
        eta_cut = (photons.isScEtaEB | photons.isScEtaEE)

        # electron veto
        #e_veto_cut = photons.electronVeto > options["e_veto"]

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

