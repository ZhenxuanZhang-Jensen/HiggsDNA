import awkward
import vector
import json
import time
import fnmatch
import copy
import inspect

import logging
logger = logging.getLogger(__name__)

from higgs_dna.systematics.systematic import EventWeightSystematic, ObjectWeightSystematic, SystematicWithIndependentCollection
from higgs_dna.systematics import photon_systematics, lepton_systematics, jet_systematics, pileup_systematics
from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.constants import NOMINAL_TAG

class SystematicsProducer():
    """
    Class to organize systematics variations.
    :param options: information for producing various systematic variations
    :type options: str, dict
    :param name: name to identify this SystematicsProducer
    :type name: str
    :param sample: Sample object for this SystematicsProducer
    :type sample: higgs_dna.samples.sample.Sample
    """
    def __init__(self, options, name = "default", sample = None):
        self.name = name
        self.options = copy.deepcopy(misc_utils.load_config(options))
        self.sample = sample

        self.weights = {}
        self.independent_collections = {}
        self.add_systematics(self.options)


    def syst_is_relevant(self, syst, syst_info):
        """
        Check if this systematic is applicable for the current sample/year

        :param syst: name of this systematic
        :type syst: str
        :param syst_info: dictionary of options for this systematic
        :type syst_info: dict
        :return: True if this systematic should be applied on this sample/year, False otherwise 
        :rtype: bool
        """
        # If neither specific years nor specific samples are specified, it is always relevant
        if "years" not in syst_info.keys() and "samples" not in syst_info.keys():
            return True

        elif "years" in syst_info.keys() and "samples" not in syst_info.keys():
            if self.sample.year in syst_info["years"]:
                return True

        elif "years" not in syst_info.keys() and "samples" in syst_info.keys():
            if self.sample.process in syst_info["samples"]:
                return True

        elif "years" in syst_info.keys() and "samples" in syst_info.keys():
            if self.sample.year in syst_info["years"] and self.sample.process in syst_info["samples"]:
                return True

        logger.debug("[SystematicsProducer : syst_is_relevant] SystematicsProducer '%s' found that syst '%s' is not applicable for year '%s' and sample '%s'." % (self.name, syst, self.sample.year, self.sample)) 

        return False

    @staticmethod
    def convert_lists_to_tuples(syst):
        """

        """
        for field in ["branches", "branch_modified", "input_collection", "target_collections"]:
            if field in syst.keys():
                if syst[field] is not None:
                    if field in ["branch_modified", "input_collection"]:
                        if isinstance(syst[field], list):
                            syst[field] = tuple(syst[field])
                    elif field == "branches":
                        for var, branch in syst[field].items():
                            if isinstance(branch, list):
                                syst[field][var] = tuple(branch)
                    elif field == "target_collections":
                        for idx, coll in enumerate(syst[field]):
                            if isinstance(coll, list):
                                syst[field][idx] = tuple(coll)

        return syst


    def add_systematics(self, systematics):
        """
        Construct Systematic objects from options dict
        """
        if "weights" in systematics.keys():
            for syst, syst_info in systematics["weights"].items():
                """
                A single weight systematic in the input config can in general correspond to multiple weight systematics in the output file.
                This happens when multiple output collections are specified.
                We create a list of weight systematics each sharing the same input collection, but in principle with different output collections,
                and run the 'produce' method just once, but the 'apply' method for all of them.
                """
                if syst in self.weights.keys(): # check if it has already been added
                    continue

                if self.sample is not None:
                    if self.sample.is_data: # weight systs should never be applicable to data
                        continue

                    if not self.syst_is_relevant(syst, syst_info): # otherwise some weight systs might only be applicable to certain MC samples
                        continue

                syst_info = self.convert_lists_to_tuples(syst_info)

                logger.info("[SystematicsProducer : add_systematics] SystematicsProducer '%s', adding weight systematic '%s'." % (self.name, syst))
                logger.debug("[SystematicsProducer : add_systematics] Details: %s" % syst_info)

                self.weights[syst] = [] # a single weight syst in the config can result in multiple systematics if it has multiple target collections

                if syst_info["type"] == "object": 
                    if syst_info["method"] == "from_function":
                        syst_info["branches"] = None
                    elif syst_info["method"] == "from_branch":
                        syst_info["function"] = None
                    if "input_collection" in syst_info.keys():
                        if not "target_collections" in syst_info.keys():
                            logger.debug("[SystematicsProducer : add_systematics] No target collections specified for syst '%s', using the input collection '%s' as the target" % (syst, syst_info["input_collection"]))
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
                                target_collection = target_collection,
                                sample = self.sample
                        )

                        if "modifies_taggers" in syst_info.keys():
                            if isinstance(syst_info["modifies_taggers"], dict): # you gave a map of target_collection : taggers
                                weight_syst.modifies_taggers = syst_info["modifies_taggers"][target_collection]
                            elif isinstance(syst_info["modifies_taggers"], list): # you gave a list of taggers
                                weight_syst.modifies_taggers = syst_info["modifies_taggers"]
                        if "kwargs" in syst_info.keys():
                            for kwarg, val in syst_info["kwargs"].items():
                                setattr(weight_syst, kwarg, val)

                        self.weights[syst].append(weight_syst)           

                elif syst_info["type"] == "event":
                    if syst_info["method"] == "from_function":
                        syst_info["branches"] = None
                    elif syst_info["method"] == "from_branch":
                        syst_info["function"] = None

                    if not "requires_branches" in syst_info.keys():
                        syst_info["requires_branches"] = None

                    weight_syst = EventWeightSystematic(
                            name = syst,
                            method = syst_info["method"],
                            modify_central_weight = syst_info["modify_central_weight"],
                            branches = syst_info["branches"],
                            function = syst_info["function"],
                            requires_branches = syst_info["requires_branches"],
                            sample = self.sample
                    )

                    if "modifies_taggers" in syst_info.keys():
                        weight_syst.modifies_taggers = syst_info["modifies_taggers"]

                    if "kwargs" in syst_info.keys():
                        for kwarg, val in syst_info["kwargs"].items():
                            setattr(weight_syst, kwarg, val) 

                    self.weights[syst].append(weight_syst)

        if "independent_collections" in systematics.keys():
            for syst, syst_info in systematics["independent_collections"].items():
                if syst in self.independent_collections.keys(): # check if it has already been added
                    continue

                if not self.syst_is_relevant(syst, syst_info):
                    continue 

                if self.sample is not None:
                    if self.sample.is_data: # FIXME: need to allow for ic systs to modify central value in data but not result in ICs with variations 
                        continue 

                syst_info = self.convert_lists_to_tuples(syst_info)

                logger.info("[SystematicsProducer : add_systematics] SystematicsProducer '%s', adding systematic with independent collection '%s'." % (self.name, syst))
                logger.debug("[SystematicsProducer : add_systematics] Details: %s" % syst_info) 

                if syst_info["method"] == "from_branch":
                    syst_info["function"] = None
                if syst_info["method"] == "from_function":
                    syst_info["branches"] = None

                if "additive" in syst_info.keys():
                    additive = syst_info["additive"]
                else:
                    additive = False

                self.independent_collections[syst] = SystematicWithIndependentCollection(
                        name = syst,
                        method = syst_info["method"],
                        branch_modified = syst_info["branch_modified"],
                        branches = syst_info["branches"],
                        function = syst_info["function"],
                        sample = self.sample,
                        additive = additive
                )

                # Add any relevant arguments for calculating this systematic
                if "kwargs" in syst_info.keys():
                    for kwarg, val in syst_info["kwargs"].items():
                        setattr(self.independent_collections[syst], kwarg, val)


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
           
            # Check if the input collection is already produced (e.g. trigger SF relies on Lead/Sublead photons which are not determined until after running the diphoton tagger
            if hasattr(weight_systs[0], "input_collection"):
                missing_fields = awkward_utils.missing_fields(events, [weight_systs[0].input_collection])
                if missing_fields:
                    continue

            if hasattr(weight_systs[0], "requires_branches"):
                if weight_systs[0].requires_branches is not None:
                    missing_fields = awkward_utils.missing_fields(events, weight_systs[0].requires_branches)
                    if missing_fields:
                        continue

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

                if NOMINAL_TAG in weight_syst.is_applied.keys(): 
                    if weight_syst.is_applied[NOMINAL_TAG]: # if we already apply the weight syst before creating ICs, mark it as applied for all ICs (to avoid double-applying later on)
                        weight_syst.is_applied_all = True


        for name, ic_syst in self.independent_collections.items():
            ics = ic_syst.produce(events)
            if NOMINAL_TAG in ics.keys():
                events = ics[NOMINAL_TAG]
            for ic_name, ic in ics.items():
                if ic_name == NOMINAL_TAG:
                    continue
                events_with_syst[name + "_" + ic_name] = ic

        events_with_syst[NOMINAL_TAG] = events

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
                    if weight_syst.is_applied_all:
                        continue
                    if name in weight_syst.is_applied.keys():
                        if weight_syst.is_applied[name]:
                            continue

                    reset = False
                    if not weight_syst.is_produced:
                        reset = True
                        syst_events = weight_syst.produce(syst_events, central_only = name != NOMINAL_TAG)

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
                            central_only = name != NOMINAL_TAG,
                            mask = mask
                    )

                    events_with_syst[name] = syst_events
                
                    if reset:
                        weight_syst.is_produced = False

        self.summary = { "weights" : {}, "independent_collections" : {} }
        for weight_name, weight_systs in self.weights.items():
            for weight_syst in weight_systs:
                for variation, info in weight_syst.summary[NOMINAL_TAG].items():
                    self.summary["weights"][info["branch"]] = weight_syst.summary

        return events_with_syst
