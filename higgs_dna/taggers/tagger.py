import awkward
import numpy
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
        selection, events_updated = self.get_selection(NOMINAL_TAG, events)
        return events_updated[selection]


    def run(self, events): 
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

        for syst_tag, syst_events in events.items():
            self.current_syst = syst_tag
            selection, syst_events_updated = self.get_selection(syst_tag, syst_events)
            self.selection[syst_tag] = selection
            self.events[syst_tag] = syst_events_updated

            logger.debug("[Tagger] %s : event set : %s : %d (%d) events before (after) selection" % (self.name, syst_tag, len(syst_events), len(syst_events_updated[self.selection[syst_tag]])))


        return self.selection, self.events


    def set_selection(self, selection):
        """
        Update selection (TagSequence will do this to remove overlap
        with other tags)
        :param selection: boolean array for each set of events
        :type selection: dict
        """
        for name in selection.keys():
            self.selection[name] = selection[name]


    def get_selection(self, syst_tag, syst_events): 
        """
        Return boolean array of events to be selected for a given event set.
        Checks if the selection has already been calculated and
        calculates it if not.
        :param syst_tag: name of current systematic variation with independent collection
        :type syst_tag: str
        :param syst_events: events for current systematic variation with independent collection
        :type syst_events: awkward.Array
        """

        if syst_tag in self.selection.keys():
            return self.selection[syst_tag], self.events[syst_tag] 

        else:
            return self.calculate_selection(syst_events) 

        
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

        iter = 0
        for name, result in zip(names, results):
            iter += 1
            if awkward.count(result) > 0:
                individual_eff = float(awkward.sum(result)) / float(awkward.count(result))
            else:
                individual_eff = 0.
            self.cut_summary[cut_type][name] = {
                    "individual_eff" : float(individual_eff)
                    #TODO: add eff as N-1 cut
            }
            logger.debug("[Tagger] : %s, syst variation : %s, cut type : %s, cut : %s, indiviual efficiency : %.4f"% (self.name, self.current_syst, cut_type, name, individual_eff))
            # get the combined cuts and names
            if(iter==1):
                _tmp_cut = result
                _tmp_name = name
            else:
                _tmp_cut = numpy.logical_and(_tmp_cut, result)
                _tmp_name += " & " + name
            if awkward.count(_tmp_cut) > 0:
                ncandi_per_event = awkward.num(_tmp_cut[_tmp_cut==True],axis=-1) 
                candi_event=_tmp_cut[ncandi_per_event!=0]
                if type(candi_event) == bool:
                    print(_tmp_cut)
                # for event selection level, like "at least one diphoton pair", _tmp_cut is 1D array, ncandi_per_event is the number of the events which contians at least one diphoton pair
                    combined_eff = float(ncandi_per_event) / float(len(_tmp_cut))
                else:
                    n_candi_event = len(candi_event)
                    combined_eff = float(n_candi_event) / float(len(_tmp_cut))
                    # combined_candieff = float(awkward.sum(_tmp_cut)) / float(awkward.count(_tmp_cut))
            else:
                combined_eff = 0.
            self.cut_summary[cut_type][_tmp_name]={
                "combined eff": float(combined_eff)
            }
            # logger.debug("[Tagger] : %s, syst variation : %s, cut type : %s, cut : %s, combined candi/sum_candi efficiency : %.4f"% (self.name, self.current_syst, cut_type, _tmp_name, combined_candieff))
            logger.debug("[Tagger] :  \nc_cut : %s\neff : %.4f"% (_tmp_name, combined_eff))



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
