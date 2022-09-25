from higgs_dna.utils.logger_utils import setup_logger
logger = setup_logger("DEBUG", None)

# Get sample
from higgs_dna.samples.sample_manager import SampleManager
sample_manager = SampleManager(
        catalog = "metadata/samples/tth_tutorial.json",
        sample_list = ["VH_M125"],
        years = ["2017"]
)

vh_sample = sample_manager.get_samples()[0]

# Load a tth file
import uproot, awkward

vh_file = vh_sample.files[0].name
with uproot.open(vh_file) as f:
    tree = f["Events"]
    events = tree.arrays(library = "ak", how = "zip")

events = vh_sample.prep(events)

# Run diphoton preselection
from higgs_dna.taggers.diphoton_tagger import DiphotonTagger

diphoton_tagger = DiphotonTagger()
diphoton_events = diphoton_tagger.select(events)

# Apply L1-prefiring SFs
from higgs_dna.systematics.systematic import EventWeightSystematic

l1_prefiring = EventWeightSystematic(
        name = "L1Prefire",
        method = "from_branch",
        branches = {
            "central" : "L1PreFiringWeight_Nom",
            "up" : "L1PreFiringWeight_Up",
            "down" : "L1PreFiringWeight_Dn"
        },
        modify_central_weight = True,
        sample = vh_sample
)
diphoton_events = l1_prefiring.produce(diphoton_events) # adds branches
diphoton_events = l1_prefiring.apply(diphoton_events) # modifies central weight value

# can now access:
#   diphoton_events.weight_L1Prefire_{central/up/down}

# Calculate electron ID SFs
from higgs_dna.systematics.systematic import ObjectWeightSystematic
from higgs_dna.systematics import lepton_systematics

electron_id_sf = ObjectWeightSystematic(
        name = "ElectronID",
        method = "from_function",
        function = {
            "module_name" : "higgs_dna.systematics.lepton_systematics",
            "name" : "electron_id_sf"
        },
        input_collection = "Electron",
        target_collection = "SelectedElectron",
        modify_central_weight = True,
        sample = vh_sample
)
electron_id_sf.working_point = "wp90iso"

diphoton_events = electron_id_sf.produce(diphoton_events)

# Select electrons
from higgs_dna.selections.lepton_selections import select_electrons

diphoton_events["SelectedElectron"] = diphoton_events.Electron[
        select_electrons(
            electrons = diphoton_events.Electron,
            options = {
                "pt" : 20.0,
                "eta" : 2.4,
                "id" : "WP90" # 90% efficiency working point
            },
            clean = {
                "photons" : {
                    "objects" : diphoton_events.Diphoton.Photon,
                    "min_dr" : 0.2
                }
            }
        )
]

# Apply event-level SF based on selected electrons
diphoton_events = electron_id_sf.apply(diphoton_events)

# can now access:
#   diphoton_events.weight_ElectronID_{central/up/down}

# Calculate Z->ee candidates
import vector
vector.register_awkward()

diphoton_events["SelectedElectron"] = awkward.Array(diphoton_events.SelectedElectron, with_name = "Momentum4D")
diphoton_events["dielectron_pairs"] = awkward.combinations(
        diphoton_events.SelectedElectron,
        2,
        fields = ["LeadElectron", "SubleadElectron"]
)
diphoton_events[("dielectron_pairs", "ZCandidate")] = diphoton_events.dielectron_pairs.LeadElectron + diphoton_events.dielectron_pairs.SubleadElectron

# Apply OS and m_ll cuts
os_cut = diphoton_events.dielectron_pairs.LeadElectron.charge * diphoton_events.dielectron_pairs.SubleadElectron.charge == -1
mass_cut = (diphoton_events.dielectron_pairs.ZCandidate.mass > 86.) & (diphoton_events.dielectron_pairs.ZCandidate.mass < 96.)
diphoton_events["dielectron_pairs"] = diphoton_events.dielectron_pairs[os_cut & mass_cut]

# Get gen-level Z->ee info
from higgs_dna.selections.gen_selections import select_x_to_yz
diphoton_events["GenZee"] = select_x_to_yz(
        gen_part = diphoton_events.GenPart,
        x_pdgId = 23,
        y_pdgId = 11,
        z_pdgId = 11
)

awkward.to_parquet(diphoton_events, "my_output.parquet")
