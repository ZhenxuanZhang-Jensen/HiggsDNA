import awkward
import numpy
import json

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger
from higgs_dna.constants import GOLDEN_JSON
from higgs_dna.utils import misc_utils

class GoldenJsonTagger(Tagger):
    """
    Tagger to select events in the golden json.
    Should only be applied to data.
    NOTE: TagSequence will automatically add this for data if a corresponding json file for the relevant year is listed in higgs_dna.constants.GOLDEN_JSON. You should not have to add it yourself.
    """
    def __init__(self, name = "golden_json_tagger", options = {}, is_data = True, year = None):
        super(GoldenJsonTagger, self).__init__(name, options, is_data, year)

        if not self.is_data:
            logger.exception("[GoldenJsonTagger : __init__] Golden json should only be applied on data.")
            raise RuntimeError()

        if self.year not in GOLDEN_JSON.keys():
            logger.warning("[A corresponding json file was not found in higgs_dna.constants.GOLDEN_JSON for year '%s'. This tagger will apply a dummy selection." % (self.year))
            self.json_file = None
            self.golden_json = None

        else:
            self.json_file = misc_utils.expand_path(GOLDEN_JSON[self.year])
            with open(self.json_file, "r") as f_in:
                self.golden_json = json.load(f_in)

        logger.info("Created golden json tagger")
    
    def calculate_selection(self, events):
        logger.info("Running golden json tagger")
        
        # If no valid golden json file, return all True
        if self.json_file is None:
            cut = awkward.ones_like(events.run, dtype=bool)#abs(events.run) >= 0
        else:
            # Start with dummy all False cut and mark passing run x lumis as True
            cut = awkward.zeros_like(events.run, dtype=bool) 
 
            # Find all runs present in events
            runs = numpy.unique(awkward.to_numpy(events.run))
            logger.info("Runs: %s" % str(runs))

            # Loop through each run x lumi pair
            for run in runs:
                # If run is not in json file, all lumis are bad
                if str(run) not in self.golden_json.keys():
                    continue # cut is initialized to all False, so no need to update
                
                for lumi_range in self.golden_json[str(run)]:
                    cut = cut | ((events.run == run) & (events.luminosityBlock >= lumi_range[0]) & (events.luminosityBlock <= lumi_range[1]))

        self.register_cuts(
                names = ["golden json"],
                results = [cut]
        )

        return cut, events


