from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.selections import lepton_selections, jet_selections, fatjet_selections
from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
import logging
import awkward
from matplotlib.pyplot import jet
import vector

vector.register_awkward()

logger = logging.getLogger(__name__)


DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {
    "electrons": {
        "pt": 20.0,
        "dr_photons": 0.2
    },
    "muons": {
        "pt": 21.0,
        "dr_photons": 0.2
    },
    "jets": {
        "pt": 25.0,
        "eta": 2.4,
        "dr_photons": 0.4,
        "dr_electrons": 0.4,
        "dr_muons": 0.4,
    },
    "fatjets": {
        "pt": 250.0,
        "eta": 2.4,
        "dr_photons": 0.4,
        "dr_electrons": 0.4,
        "dr_muons": 0.4,
    },
    "photon_id": 0.0,
    "btag_wp": {
        "2016": 0.3093,
        "2017": 0.3040,
        "2018": 0.2783
    }
}


class HHWW_Preselection(Tagger):
    """
    HHWW Preselection tagger for tutorial
    """

    def __init__(self, name, options={}, is_data=None, year=None):
        super(HHWW_Preselection, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                original=DEFAULT_OPTIONS,
                new=options
            )

    def calculate_selection(self, events):
        # Electrons
        electron_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["electrons"]["dr_photons"]
                }
            },
            name="SelectedElectron",
            tagger=self
        )

        electrons = awkward_utils.add_field(
            events=events,
            name="SelectedElectron",
            data=events.Electron[electron_cut]
        )

        # Muons
        muon_cut = lepton_selections.select_muons(
            muons=events.Muon,
            options=self.options["muons"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["muons"]["dr_photons"]
                }
            },
            name="SelectedMuon",
            tagger=self
        )

        muons = awkward_utils.add_field(
            events=events,
            name="SelectedMuon",
            data=events.Muon[muon_cut]
        )
        # --------------------- if fatjet branches are not empty --------------------- #
        # if len(events.FatJet.pt > 0 ):
        # Fat jets
        fatjet_cut = fatjet_selections.select_fatjets(
            fatjets = events.FatJet,
            options = self.options["fatjets"],
            clean = {
                "photons" : {
                    "objects" : events.Diphoton.Photon,
                    "min_dr" : self.options["jets"]["dr_photons"]
                },
                "electrons" : {
                    "objects" : events.SelectedElectron,
                    "min_dr" : self.options["jets"]["dr_electrons"]
                },
                "muons" : {
                    "objects" : events.SelectedMuon,
                    "min_dr" : self.options["jets"]["dr_muons"]
                    }
                },
            name = "SelectedFatJet",
            tagger = self
        )
        fatjets = awkward_utils.add_field(
            events = events,
            name = "SelectedFatJet",
            data = awkward.Array(events.FatJet[fatjet_cut], with_name = "Momentum4D")
        )   
        # Jets
        jet_cut = jet_selections.select_jets(
            jets=events.Jet,
            options=self.options["jets"],
            clean={
                "photons": {
                    "objects": events.Diphoton.Photon,
                    "min_dr": self.options["jets"]["dr_photons"]
                },
                "electrons": {
                    "objects": events.SelectedElectron,
                    "min_dr": self.options["jets"]["dr_electrons"]
                },
                "muons": {
                    "objects": events.SelectedMuon,
                    "min_dr": self.options["jets"]["dr_muons"]
                }
            },
            name="SelectedJet",
            tagger=self
        )

        jets = awkward_utils.add_field(
            events=events,
            name="SelectedJet",
            data=events.Jet[jet_cut]
        )
        # jets_7WithEachEvent = awkward.pad_none(jets, 7, clip=True)
        # FIXME: I have no other method, but just loop event by event to sort 2jets combined with 4 jets mass and get the W1 and W2 mass.
        # for i in range(len(jets_7WithEachEvent)):
        #     for j in range(len(jets_7WithEachEvent)):
        #         jets_7WithEachEvent[i] + jets_7WithEachEvent[j]

        # ---------------------------------------------------------------------------- #
        #              TODO try to add order with WinvM jets in output.parquet              #
        # ---------------------------------------------------------------------------- #
        awkward_utils.add_object_fields(
            events=events,
            name="jet",
            objects=jets,
            n_objects=7,
            dummy_value=-999
        )
        bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis=1, ascending=False)]
        bjets = bjets[bjets.btagDeepFlavB > self.options["btag_wp"][self.year]]

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        electrons = awkward.Array(electrons, with_name="Momentum4D")
        muons = awkward.Array(muons, with_name="Momentum4D")

        # Preselection
        n_electrons = awkward.num(electrons)
        n_muons = awkward.num(muons)
        n_leptons = n_electrons + n_muons

        n_jets = awkward.num(jets)
        # n_bjets = awkward.num(bjets)

        photon_id_cut = (events.LeadPhoton.mvaID > self.options["photon_id"]) & (
            events.SubleadPhoton.mvaID > self.options["photon_id"])

        # Hadronic presel
        hadronic = (n_leptons == 0) & (n_jets >= 4)

        # Semi-Leptonic presel
        # Semileptonic = (n_leptons == 1) & (n_jets >= 2)

        # Fully-Leptonic presel
        # FulllyLeptonic = (n_leptons >= 2)

        # presel_cut = (hadronic | Semileptonic | FulllyLeptonic)  & photon_id_cut
        presel_cut = (hadronic) & photon_id_cut

        self.register_cuts(
            names=["hadronic presel", "photon ID cut", "all"],
            results=[hadronic, photon_id_cut, presel_cut]
        )
        return presel_cut, events
