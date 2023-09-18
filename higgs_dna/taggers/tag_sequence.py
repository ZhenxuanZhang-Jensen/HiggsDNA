import sys
import numpy
import json
import awkward
import os
import importlib

import logging
logger = logging.getLogger(__name__)

from higgs_dna.constants import NOMINAL_TAG
from higgs_dna.taggers.tagger import Tagger
from higgs_dna.utils import awkward_utils

class TagSequence():
    """
    Accepts a list of tags (Tagger instances) and applies each tag's
    selection. Priority is given in the order tags are passed.
    :param events: sets of events to perform selection on for each tagger
    :type events: dict or awkward.Array
    :param weight_variations: list of weight_variations to save for nominal set of events
    :type weight_variations: list
    :param tag_list: list of lists of higgs_dna.taggers.tagger.Tagger objects or list of lists of dicts 
    :type tag_list: list
    :param ext: string to identify the output files from this :class:`higgs-dna.taggers.tag_sequence.TagSequence` object
    :type ext: str, optional
    """
    def __init__(self, tag_list, name = "default", sample = None,output_dir=None):
        self.name = name
        self.sample = sample

        self.tag_list = []
        self.selections = {}
        self.summary = {}
        self.output_dir = output_dir

            
        logger.info("[TagSequence : __init__] Creating tag sequence with the following tags and priority:")
        for i, tag_set in enumerate(tag_list):
            m_tag_set = []
            if not isinstance(tag_set, list): # if not already given as a list, cast as a list of one tagger
                tag_set = [tag_set]
            # Make sure we only have >1 taggers if it is the last step in tag sequence
            if len(tag_set) > 1:
                if i != len(tag_list) - 1:
                    logger.exception("[TagSequence : __init__] A list of tags should only be given as the very last step of the tag sequence.")
                    raise ValueError()

            for j, tagger in enumerate(tag_set):
                if isinstance(tagger, dict):
                    m_tagger = self.create_tagger(tagger)
                elif isinstance(tagger, Tagger):
                    m_tagger = tagger
                else:
                    logger.exception("[TagSequence : __init__] Each entry in a tag set must either be a dict (from which we construct a Tagger) or a Tagger instance, not '%s' as you have provided for the %d-th tagger in the %d-th tag set." % (str(type(tagger)), j, i))
                    raise TypeError()
                m_tag_set.append(m_tagger)

            self.tag_list.append(m_tag_set)
            logger.info("[TagSequence : __init__] The %d-th tag set has taggers %s (listed in order priority will be given)." % (i, str([t.name for t in m_tag_set])))


    def create_tagger(self, config):
        module = importlib.import_module(config["module_name"])
        
        if "kwargs" not in config.keys():
            config["kwargs"] = { "name" : "my_" + config["tagger"] }

        if self.sample is not None:
            config["kwargs"]["is_data"] = self.sample.is_data
            config["kwargs"]["year"] = self.sample.year
            config['kwargs']['output_dir']=self.output_dir

        tagger = getattr(
            module,
            config["tagger"]
        )(**config["kwargs"])

        return tagger


    def run(self, events, syst_tag = NOMINAL_TAG):
        """
        Get initial tag selections (potentially with overlap),remove overlap between tags,
        add a field 'tag_idx' indicating which tag each event belongs to, and return selected events. 

        :param events: either a single events array or a dictionary of arrays corresponding to systematic variations with independent collections
        :type events: awkward.highlevel.Array or dict
        :return: events selected by the final set of taggers and a dictionary of tagger name : idx corresponding to that tag
        :rtype: dict, dict
        """

        for i, tag_set in enumerate(self.tag_list):
            n_taggers = len(tag_set)

            events = self.run_taggers(events, syst_tag, tag_set)
            if n_taggers > 1:
                events = self.orthogonalize_tags(events, syst_tag, tag_set)
            events = self.select_events(events, syst_tag, tag_set)


        self.summarize()

        return events, self.tag_idx_map


    def run_taggers(self, events, syst_tag, tag_list):
        """
        Get initial selection for each tag and allow taggers to add fields to events array.
        There is still overlap between tags at this step (if more than one tagger present).

        :param tag_list: list of taggers to run selections for
        :type tag_list: list of higgs_dna.taggers.tagger.Tagger
        """
        for idx, tagger in enumerate(tag_list):
            if tagger.name not in self.selections.keys():
                self.selections[tagger.name] = {}
            self.selections[tagger.name][syst_tag], events = tagger.run(events, syst_tag)
            self.summary[tagger.name] = { "priority" : idx }

            self.summary[tagger.name][syst_tag] = int(awkward.sum(self.selections[tagger.name][syst_tag]))
        return events


    def orthogonalize_tags(self, tag_list):
        """
        Remove overlap between tags, assigning priority as the
        order in which tags were passsed to tag_list.

        :param tag_list: list of taggers to run selections for
        :type tag_list: list of higgs_dna.taggers.tagger.Tagger
        """

        for idx, tagger in enumerate(tag_list):
            tagger_name = tagger.name
            selection = self.selections[tagger_name][syst_tag]
            if idx == 0: # first tag in sequence, no overlap removal needed 
                prev_selected_events = selection
            else:
                selection = selection & (~self.prev_selected_events)
                self.selections[tagger_name][syst_tag] = selection
                self.summary[tagger_name][syst_tag]["n_events_post_overlap_removal"] = int(awkward.sum(selection))
                logger.debug("[TagSequence : othogonalize_tags] tag : %s, syst tag : %s, %d (%d) events before (after) overlap removal" % (tagger_name, syst_tag, self.summary[tagger_name][syst_tag]["n_events_pre_overlap_removal"], self.summary[tagger_name][syst_tag]["n_events_post_overlap_removal"]))

                prev_selected_events = selection & prev_selected_events


    def select_events(self, events, syst_tag, tag_list):
        """
        Add a column to indicate which tag each event belongs to.
        Create trimmed awkward.Array for each systematic variation,
        containing only the selected columns to write to disk.

        :param tag_list: list of taggers to run selections for
        :type tag_list: list of higgs_dna.taggers.tagger.Tagger 
        """
        self.tag_idx_map = {}
        for idx, tagger in enumerate(tag_list):
            self.tag_idx_map[tagger.name] = idx


        tag_idx = numpy.ones_like(awkward.to_numpy(events.run)) * -1
        if len(events) >= 1:
            for tagger_name, idx in self.tag_idx_map.items():
                tag_idx[self.selections[tagger_name][syst_tag]] = idx
        awkward_utils.add_field(events, "tag_idx", tag_idx, overwrite = True)
        selection = events.tag_idx >= 0

        return events[selection] 


    def summarize(self):
        """
        Grab summary info from each tag and create a json file
        with TagSequence diagnostic info.
        """
        for tag_set in self.tag_list: 
            for tagger in tag_set:
                self.summary[tagger.name]["summary"] = tagger.get_summary()

        return self.summary
        file_name = "summary%s.json" % self.ext
        with open("output/" + file_name, "w") as f_out:
            json.dump(self.summary, f_out, indent=4)