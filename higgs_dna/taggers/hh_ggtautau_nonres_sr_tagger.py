import awkward
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG
from higgs_dna.utils import awkward_utils, misc_utils

DEFAULT_OPTIONS = {}

class HHggTauTauNonResSRTagger(Tagger):
    """
    Signal region tagger for the non-resonant HH->ggTauTau analysis.
    """
    def __init__(self, name, options = {}, is_data = None, year = None):
        super(HHggTauTauNonResSRTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, syst_tag, syst_events):
        #####################################
        ### HH->ggTauTau Non-resonant SRs ###
        #####################################

        # BDT evaluation

        return presel_cut, syst_events
