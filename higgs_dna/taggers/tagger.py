import awkward
import vector
import json
import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import misc_utils

from higgs_dna.constants import NOMINAL_TAG

class Tagger():
    """
    Abstract base class for taggers
    :param name: Name to identify this tagger
    :type name: str
    :param options: options with selection info
    :type options: str, dict
    :param is_data: whether this tagger is being run on data
    :type is_data: bool
    :param year: which year this tagger is being run on
    :type year: str
    """
    def __init__(self, name = "tagger", options = {}, is_data = None, year = None):
        self.name = name
        self.options = options
        self.is_data = is_data
        self.year = year

        self.selection = {}
        self.events = {}
        self.cut_summary = {}

        self.options = misc_utils.load_config(options)


    def select(self, events):
        """
        Convenience function for running tagger in standalone contexts.

        :param events: awkward array of events
        :type events: awkward.Array
        :returns: awkward array of selected events
        :rtype: awkward.Array
        """
        self.current_syst = NOMINAL_TAG
        selection, events_updated = self.calculate_selection(events)
        return events_updated[selection]


    def run(self, events, syst_tag = NOMINAL_TAG): 
        """
        Return dictionary of boolean arrays of events to be selected
        by this tagger for each systematic variation,
        along with updated dictionary of events containing any additional fields
        added by this tagger.
        :param events: sets of events to perform selection for (nominal + any independent collections for syst variations
        :type events: dict
        :return: boolean arrays of events to be selected by this tagger for each systematic variation, events dict with added fields computed by this tagger
        :rtype: dict, dict 
        """

        self.current_syst = syst_tag
   
        if not len(events) >= 1:
            logger.debug("[Tagger] %s : event set : %s : 0 events passed to tagger, skipping running this tagger." % (self.name, syst_tag))
            selection = awkward.ones_like(events, dtype=bool)
            self.selection[syst_tag] = selection

        else:
            selection, events = self.calculate_selection(events)
            self.selection[syst_tag] = selection
            logger.debug("[Tagger] %s : event set : %s : %d (%d) events before (after) selection" % (self.name, syst_tag, len(selection), awkward.sum(selection)))

        return selection, events


    def calculate_selection(self, events): 
        """
        Abstract function that should be reimplemented for each tagger.
        """

        raise NotImplementedError()


    @staticmethod
    def get_range_cut(array, ranges):
        """
        Return mask corresponding to selecting elements
        in <array> that are within the ranges (inclusive)
        specified in <ranges>.
        :param array: array of quantity to be cut on
        :type array: awkward.Array
        :param ranges: list of 2-element lists corresponding to allowed ranges of variable in <array>
        :type ranges: list
        """
        cuts = []

        for range in ranges:
            if not len(range) == 2:
                message = "Range of allowed values should be 2 numbers ([lower bound, upper bound]), but you gave %s" % str(range)
                logger.exception(message)
                raise AssertionError(message)
            cut_low = array >= range[0]
            cut_high = array <= range[1]
            cuts.append(cut_low & cut_high)

        for idx, cut in enumerate(cuts):
            if idx == 0:
                final_cut = cut
            else:
                final_cut = final_cut | cut

        return final_cut


    def register_cuts(self, names, results, cut_type = "event"):
        """
        Record a given cut in the tagger instance.
        cut_type could be event-level (default) or object-level ("photon", "muon", etc)
        :param names: names to identify cuts
        :type names: list of str or str
        :param results: boolean arrays with results of applying given cut
        :type results: list of awkward.Array
        :param cut_type: flag to indicate the type of cut
        :type cut_type: str, optional
        """
        if cut_type not in self.cut_summary.keys():
            self.cut_summary[cut_type] = {}

        if not isinstance(names, list):
            names = [names]
        if not isinstance(results, list):
            results = [results]

        for name, result in zip(names, results):
            if awkward.count(result) > 0:
                individual_eff = float(awkward.sum(result)) / float(awkward.count(result))
            else:
                individual_eff = 0.
            self.cut_summary[cut_type][name] = {
                    "individual_eff" : float(individual_eff) 
                    #TODO: add eff as N-1 cut
            }

            logger.debug("[Tagger] : %s, syst variation : %s, cut type : %s, cut : %s, efficiency : %.4f"
                    % (self.name, self.current_syst, cut_type, name, individual_eff))


    def get_summary(self):
        """
        Get dictionary of tag diagnostic info.
        :return: A dictionary of tag diagnostic info
        :rtype: dict
        """
        summary = {
                "options" : self.options,
                "cut_summary" : self.cut_summary
        }
        return summary
