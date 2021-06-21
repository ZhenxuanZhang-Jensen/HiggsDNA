import awkward
import vector
import json
import time
import fnmatch
import copy

import logging
logger = logging.getLogger(__name__)

from higgs_dna.systematics.systematic import EventWeightSystematic, ObjectWeightSystematic, SystematicWithIndependentCollection
from higgs_dna.utils import awkward_utils, misc_utils

class SystematicsProducer():
    """
    Class to organize systematics variations.
    :param options: information for producing various systematic variations
    :type options: str, dict
    """
    def __init__(self, options):
        self.options = misc_utils.load_config(options)

        logger.info("Creating systematics producer with the following systematic variations:")
        for syst, syst_info in self.options["weights"].items():
            logger.info("\tWeight variation: %s" % syst)
            logger.debug("\tDetails: %s" % syst_info)

        for syst, syst_info in self.options["independent_collections"].items():
            logger.info("\tIndependent collection: %s" % syst)
            logger.debug("\tDetails: %s" % syst_info)

        self.initialize_systematics()

    def initialize_systematics(self):
        """
        Construct Systematic objects from options dict
        """
        weights = {}
        for syst, syst_info in self.options["weights"].items():
            """
            A single weight systematic in the input config can in general correspond to multiple weight systematics in the output file.
            This happens when multiple output collections are specified.
            We create a list of weight systematics each sharing the same input collection, but in principle with different output collections,
            and run the 'produce' method just once, but the 'apply' method for all of them.
            """
            weights[syst] = [] # a single weight syst in the config can result in multiple systematics if it has multiple target collections

            if syst_info["type"] == "object": 
                if "input_collection" in syst_info.keys():
                    syst_info["branches"] = None
                    if not "target_collections" in syst_info.keys():
                        logger.debug("[SystematicsProducer : initialize_systematics] No target collections specified for syst '%s', using the input collection '%s' as the target" % (syst, syst_info["input_collection"]))
                        syst_info["target_collections"] = [syst_info["input_collection"]]

                else:
                    syst_info["input_collection"] = None
                    syst_info["target_collections"] = [None]

                for target_collection in syst_info["target_collections"]:
                    if isinstance(target_collection, tuple):
                        name = syst + "_" + "_".join(list(target_collection))
                    elif isinstance(target_collection, str):
                        name = syst + "_" + target_collection

                    weight_syst = ObjectWeightSystematic(
                            name = name,
                            method = syst_info["method"],
                            modify_central_weight = syst_info["modify_central_weight"],
                            branches = syst_info["branches"],
                            function = syst_info["function"],
                            input_collection = syst_info["input_collection"],
                            target_collection = target_collection 
                    )

                    if "modifies_taggers" in syst_info.keys():
                        if isinstance(syst_info["modifies_taggers"], dict): # you gave a map of target_collection : taggers
                            weight_syst.modifies_taggers = syst_info["modifies_taggers"][target_collection]
                        elif isinstance(syst_info["modifies_taggers"], list): # you gave a list of taggers
                            weight_syst.modifies_taggers = syst_info["modifies_taggers"]

                    weights[syst].append(weight_syst)           

            elif syst_info["type"] == "event":
                if syst_info["method"] == "from_function":
                    syst_info["branches"] = None
                elif syst_info["method"] == "from_branch":
                    syst_info["function"] = None

                weight_syst = EventWeightSystematic(
                        name = syst,
                        method = syst_info["method"],
                        modify_central_weight = syst_info["modify_central_weight"],
                        branches = syst_info["branches"],
                        function = syst_info["function"]
                )

                if "modifies_taggers" in syst_info.keys():
                    weight_syst.modifies_taggers = syst_info["modifies_taggers"]

                weights[syst].append(weight_syst)

        independent_collections = {}
        for syst, syst_info in self.options["independent_collections"].items():
            if syst_info["method"] == "from_branch":
                syst_info["function"] = None
            if syst_info["method"] == "from_function":
                syst_info["branches"] = None

            independent_collections[syst] = SystematicWithIndependentCollection(
                    name = syst,
                    method = syst_info["method"],
                    branch_modified = syst_info["branch_modified"],
                    branches = syst_info["branches"],
                    function = syst_info["function"]
            )

        self.weights = weights
        self.independent_collections = independent_collections


    def produce(self, events):
        """
        Produces systematic variations *before* the tag sequence has been run:
            1. Produce and apply weight systematics. For each weight systematic:
                1.1 Produce: compute the variations (currently, the necessary branches for computing the initial weights are assumed to always be present in the input file. This may not always be true in the future...) 
                1.2 Apply: if the necessary branches for applying variations are present, translate this to a per-event weight (for ObjectWeightSystematics) and modify the central weight (if specified)
            2. Produce systematics with independent collections.
                - note: independent collections are assumed to always be computable with the branches present in the input nanoAOD file (I can't think of a case where this wouldn't be true)

        :param events: events array to use for calculating systematic variations
        :type events: awkward.highlevel.Array
        :return: dictionary of the format "syst_with_ic_variation" : events_with_syst_variation, where each events array has weight systematics produced and (possibly) applied
        :rtype: dict of str : awkward.highlevel.Array
        """
        
        events_with_syst = {}
 
        for name, weight_systs in self.weights.items():
            logger.debug("[SystematicsProducer : produce] Producing weights for weight variation: %s" % name)
            events = weight_systs[0].produce(events)
            for idx, weight_syst in enumerate(weight_systs):
                weight_syst.is_produced = True
                if idx >= 1: # need to copy over the calculated branches
                    weight_syst.branches = weight_systs[0].branches

                # Try to apply this weight syst, otherwise we wait til after the taggers
                if hasattr(weight_syst, "modifies_taggers"):  
                    continue
                if not hasattr(weight_syst, "target_collection"): 
                    events = weight_syst.apply(events)
                elif weight_syst.target_collection is None:
                    events = weight_syst.apply(events)
                else:
                    missing_fields = awkward_utils.missing_fields(events, [weight_syst.target_collection])
                    if not missing_fields:
                        events = weight_syst.apply(events)
            
        events_with_syst["nominal"] = events

        for name, ic_syst in self.independent_collections.items():
            ics = ic_syst.produce(events)
            for ic_name, ic in ics.items():
                events_with_syst[name + "_" + ic_name] = ic

        return events_with_syst

    
    def apply_remaining_weight_systs(self, events_with_syst, tag_idx_map):
        """
        Applies weight systematics which were not able to be applied before running the tag sequence.
        This should only be run *after* running a tag sequence.

        :param events_with_syst: dictionary of the format "syst_with_ic_variation" : events_with_syst_variation
        :type events_with_syst: dict of str : awkward.highlevel.Array
        :param tag_idx_map: dictionary of the format 'tagger_name' : tag_idx. This allows us to apply weight systematics only to events which were tagged by specific taggers.
        :type tag_idx_map: dict of str : int
        """
        
        for name, syst_events in events_with_syst.items():
            for weight_name, weight_systs in self.weights.items():
                for weight_syst in weight_systs:
                    if name in weight_syst.is_applied.keys():
                        if weight_syst.is_applied[name]:
                            continue

                    if hasattr(weight_syst, "modifies_taggers"):
                        mask = syst_events.tag_idx < 0 # initialize to all False
                        for tagger_name in weight_syst.modifies_taggers:
                            if tagger_name not in tag_idx_map.keys():
                                message = "[SystematicsProducer : apply_remaining_weight_syst] Weight systematic: %s was specified to modify the tagger '%s', but it is not in the list of tagger names we got from TagSequence. List is %s." % (weight_name, tagger_name, str(tag_idx_map.keys()))
                                logger.exception(message)
                                raise ValueError(message)
                            mask = (mask) | (syst_events.tag_idx == tag_idx_map[tagger_name])

                    else:
                        mask = None

                    syst_events = weight_syst.apply(
                            events = syst_events,
                            syst_tag = name,
                            central_only = name != "nominal",
                            mask = mask
                    )
                    events_with_syst[name] = syst_events

        self.summary = { "weights" : {}, "independent_collections" : {} }
        for weight_name, weight_systs in self.weights.items():
            for weight_syst in weight_systs:
                for variation, info in weight_syst.summary["nominal"].items():
                    self.summary["weights"][info["branch"]] = weight_syst.summary

        return events_with_syst
