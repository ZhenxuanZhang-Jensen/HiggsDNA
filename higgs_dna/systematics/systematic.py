import sys
import awkward

import logging
logger = logging.getLogger(__name__)

from inspect import signature

from higgs_dna.utils import awkward_utils
from higgs_dna.constants import CENTRAL_WEIGHT, NOMINAL_TAG

class Systematic():
    """
    Abstract base class for any systematic: scale factors/corrections + their variations, and per-event/object variations,
    There are two main types of systematic:
        1. WeightSystematic: a correction and/or uncertainty on a weight
            1.1 EventWeightSystematic: a correction and/or uncertainty on the *central event weight*
                Examples:
                    - event weight correction: pileup reweighting scale factor
                    - event weight uncertainty: QCD scale/PDF/alpha_s weights. Up/down variations for the weight in each event, but does not modify the central event weight
                    - event weight correction and uncertainty: pileup reweighting scale factor and up/down variations
            1.2 ObjectWeightSystematic: a correction and/or uncertainty on a *per-object weight*. Almost always translated into a per-event weight.
                Examples:
                    - lepton ID SFs: the scale factor ("correction") is calculated on a per-lepton basis. This gets translated into a per-event weight by multiplying the scale factors of each lepton selected by a tagger in each event 
                    - jet b-tag SFs
        2. SystematicWithIndependentcollection: varies a per-event (e.g. MET resolution up/down) or per-object (e.g. photon scales/smearings up/down) quantity, and results in an entirely new set of events. 

    There are two ways of calculating a systematic, specified through the 'method' arg:
        1. 'from_branch' : read the variations ('central','up','down', etc) for this systematic directly from a branch in the input ntuple
        2. 'from_function' : calculate this systematic on-the-fly. If 'from_function' is specified, you will also need to implement the calculation in the 'produce' method.

    :param name: name to identify this systematic
    :type name: str
    :param method: method of calculating this systematic, either "from_branch" or "from_function"
    :type method: str
    :param branches: which branches to read the central/up/down variations from, defaults to None
    :type branches: dict (, optional)
    :param function: function to use to calculate the central/up/down variations for this syst. Provided in the form of { "module" : <module under higgs_dna.systematics, "name" : <name of function within module> }, defaults to None
    :type function: dict (, optional)
    :param sample: Sample object for this Systematic
    :type sample: higgs_dna.samples.sample.Sample
    """
    def __init__(self, name, method, branches = None, function = None, sample = None):
        self.name = name
        self.method = method
        self.branches = branches
        self.function = function
        self.sample = sample
        if sample is not None:
            self.process = sample.process
            self.is_data = sample.is_data
            self.year = sample.year


        if not (self.method == "from_branch" or self.method == "from_function"):
            message = "[Systematic : __init__] Systematic: %s, supported methods for constructing a Systematic are 'from_branch' and 'from_function', not '%s' as you passed." % (self.name, self.method)
            logger.exception(message)
            raise TypeError(message)

        if self.method == "from_branch":
            if self.branches is None:
                message = "[Systematic : __init__] Systematic: %s, method was selected as 'from_branch', but 'branches' argument was not provided." % self.name
                logger.exception(message)
                raise ValueError(message)

            if not isinstance(self.branches, dict):
                message = "[Systematic : __init__] Systematic: %s, for systematics with method 'from_branch', the 'branches' argument must be a dictionary of the form { <variation> : <branch>, ... }, not %s as you have passed." % (self.name, str(type(self.branches)))
                logger.exception(message)
                raise TypeError(message)

        if self.method == "from_function":
            if self.function is None:
                message = "[Systematic : __init__] Systematic: %s, method was selected as 'from_function' but 'function' argument was not provide." % self.name
                logger.exception(message)
                raise ValueError(message)

            if not isinstance(self.function["module_name"], str):
                message = "[Systematic : __init__] Systematic: %s, when specifying 'from_function', the module containing the function should be given as a string, e.g. 'higgs_dna.systematics.lepton_systematics', not type '%s' as you have given" % (self.name, str(type(self.function["module_name"])))
                logger.exception(message)
                raise TypeError(message)

            if self.function["module_name"] not in sys.modules:
                message = "[Systematic : __init__] Systematic: %s, method was selected as 'from_function' with module '%s', but the module has not been imported! Add an import statement to your analysis script." % (self.name, self.function["module_name"])
                logger.exception(message)
                raise ImportError(message)
            else:
                self.function["module"] = sys.modules[self.function["module_name"]]

            if not hasattr(self.function["module"], self.function["name"]):
                message = "[Systematic : __init__] Systematic: %s, method was selected as 'from_function' but function '%s' was not found in module %s." % (self.name, self.function["name"], self.function["module"])
                logger.exception(message)
                raise NotImplementedError(message)
 

        self.is_produced = False

        self.summary = {}


    def check_fields(self, events):
        """
        Check if necessary fields for calculating this Systematic are present in an events array.

        :param events: events array to check for fields in
        :type events: awkward.highlevel.Array

        :raises ValueError: if the fields are not present
        """
        
        missing_fields = awkward_utils.missing_fields(events, [field for var, field in self.branches.items()])
        if missing_fields:
            message = "[Systematic : check_fields] Systematic: %s, the events array is missing the following fields: %s, which were specified as the input branches to calculate this systematic." % (self.name, str(missing_fields))
            logger.exception(message)
            raise ValueError(message)


    def get_function_kwargs(self, events, central_only = None):
        """
        Infer what arguments should be passed to the function which calculates the central/up/down variations for a Systematic.
        Only applicable to Systematic instances with method = "from_function"

        :param events: events to calculate variations from
        :type events: awkward.highlevel.Array
        :central_only: whether to calculate only the central variation, defaults to None
        :type central_only: bool (, optional)
        :return: args to pass to function which calculates the variations for this Systematic
        :rtype: dict
        """
        
        kwargs = {}
        function_args = signature(
                getattr(
                    self.function["module"],
                    self.function["name"]
                )
        ).parameters

        for arg in function_args:
            if arg == "events":
                kwargs[arg] = events
            elif arg == "central_only":
                kwargs[arg] = central_only
            else:
                if hasattr(self, arg):
                    kwargs[arg] = getattr(self, arg)
                else:
                    message = "[Systematic : get_function_kwargs] Systematic: %s, for function '%s' from module <%s>, we found an argument '%s' that is not present as an attribute of this class and we did not know how to otherwise set this argument. This may lead to unintended behavior or crashes!" % (self.name, self.function["name"], self.function["module"], arg)
                    logger.warning(message)

        logger.debug("[Systematic : get_function_kwargs] Systematic: %s, for function '%s' from module <%s>, we are passing arguments as:" % (self.name, self.function["name"], self.function["module"]))
        for arg, val in kwargs.items():
            if arg == "events":
                continue
            logger.debug("\t\t %s : %s" % (arg, str(val)))

        return kwargs


    def produce(self, events):
        """
        Abstract function that should be implemented for each systematic.
        Reads (in the case of 'method' = 'from_branch') or calculates (in the case of 'method' = 'from_function')
        the relevant quantities for this systematic and either adds these as fields in the events array (for WeightSystematic)
        or returns one or more independent collections (for SystematicWithIndependentCollection).

        :param events: events array for which this systematic should be calculated on and/or applied to
        :type events: awkward.highlevel.Array
        """
        raise NotImplementedError()


class WeightSystematic(Systematic):
    """
    Abstract class for a systematic which affects weights but does not change any other event content. 
    Derives from higgs_dna.systematics.systematic.Systematic

    :param modify_central_weight: whether to multiply the central event weight by the central weight for this systematic
    :type modify_central_weight: bool
    :param requires_branches: branch(es) that are required to compute this systematic
    :type requires_branches: list of str, defaults to None
    """
    def __init__(self, name, method, modify_central_weight, branches = None, requires_branches = None, function = None, sample = None):
        super(WeightSystematic, self).__init__(name, method, branches, function, sample)
        self.modify_central_weight = modify_central_weight
        self.requires_branches = requires_branches

        self.is_applied = {} # a weight will often be applied to multiple different sets of events (corresponding to different independent collections)
        self.is_applied_all = False # will be set to true if the weight is applied directly on the nominal events (i.e. before producing the independent collections) to avoid duplicate application on systematics with independent collections

    def produce(self, events, central_only = False):
        """
        Calculate the central/up/down variations for this WeightSystematic and add these as fields to the events array.

        :param events: events array which will have fields added to it
        :type events: awkward.highlevel.Array
        :param central_only: whether to only calculate the central variation for this WeightSystematic, defaults to False
        :type central_only: bool (, optional)
        :return: events array with variations added as fields
        :rtype: awkward.highlevel.Array
        """
        if self.method == "from_branch":
            self.check_fields(events)

        elif self.method == "from_function":
            self.branches = {}
            variations = self.calculate_variations(events, central_only)
            for variation, values in variations.items():
                if isinstance(variation, str): # plain field in array, e.g. events.weight_syst1_up
                    name = "weight" + "_" + self.name + "_" + variation
                    self.branches[variation] = name
                elif isinstance(variation, tuple): # nested field in array, e.g. events.Photon.weight_syst1_up
                    name = tuple((variation[0], "weight" + "_" + self.name + "_" + variation[1])) 
                    self.branches[variation[1]] = name
                logger.debug("[WeightSystematic : produce] WeightSystematic: %s, adding field %s to events array" % (self.name, name))
                awkward_utils.add_field(events, name, values)

        for variation, branch in self.branches.items():
            if central_only and not variation == "central":
                continue
            self.summary[branch] = {
                    "mean" : awkward.mean(events[branch]),
                    "std" : awkward.std(events[branch])
            }
            logger.debug("[WeightSystematic : produce] WeightSystematic: %s, variation: %s has mean: %.4f and std. dev. %.4f" % (self.name, variation, self.summary[branch]["mean"], self.summary[branch]["std"]))

        self.is_produced = True
        return events


    def calculate_variations(self, events, central_only):
        """
        Call the function which was specified to calculate the central/up/down variations for this WeightSystematic,
        and check to make sure that the output is of the expected format.

        :param events: events array to use for calculating variations for this WeightSystematic
        :type events: awkward.highlevel.Array
        :param central_only: whether to only calculate the central variation
        :type central_only: bool
        :return: dictionary of variations in the format "weight_var" : weight_var where "var" would be "up"/"down"/"central". In the case of an ObjectWeightSystematic, the keys will be tuples with the input collection as the first entries and the actual weight name as the last entry 
        :rtype: dict
        """
        kwargs = self.get_function_kwargs(events, central_only)
        variations = getattr(
                self.function["module"],
                self.function["name"]
        )(**kwargs)

        for variation, values in variations.items():
            if not (isinstance(variation, str) or isinstance(variation, tuple)):
                message = "[WeightSystematic : produce] WeightSystematic: %s, each key in variations dict should be a <str> or <tuple> of <str>, not %s as we got." % (self.name, str(type(variation)))
                logger.exception(message)
                raise TypeError(message)

            if not isinstance(values, awkward.highlevel.Array):
                message = "[WeightSystematic : produce] WeightSystematic: %s, each value in variations dict should be an <awkward.highlevel.Array>, not %s as we got." % (self.name, str(type(values)))
                logger.exception(message)
                raise TypeError(message)

        return variations

    def apply(self, events, syst_tag = NOMINAL_TAG, central_only = False, mask = None, add_fields = False):
        """
        Apply this weight systematic to the events array.
        Applying consists of two steps:
            1. Adding the per-event variation to the events array. For EventWeightSystematics, this is already done in the produce() function. For ObjectWeightSystematics, the per-object weights computed in the produce() function will be translated to per-event weights.
            2. Modify central weight: if specified, multiply the central weight in the events array by the central variation of this weight systematic 

        :param events: events for current systematic variation with independent collection
        :type events: awkward.highlevel.Array
        :param syst_tag: name of current systematic variation with independent collection, defaults to 'nominal'
        :type syst_tag: str (, optional)
        :param central_only: whether to only calculate the central variation
        :type central_only: bool
        :param mask: boolean array indicating which events to modify the central weight for. If no mask is provided, will modify the central weight for all events, defaults to None
        :type mask: awkward.highlevel.Array (, optional)
        :param add_fields: whether to add this weight variation as an additional field in the events array, defaults to False
        :type add_fields: bool (, optional)
        :return: events for the current systematic variation with the per-event variations added as fields in the array (if applicable) and the central weight modified (if applicable)
        :rtype: awkward.highlevel.Array
        """
        if not self.is_produced:
            message =  "[WeightSystematic : apply] WeightSystematic: %s must be produced before applying it!" % self.name
            logger.exception(message)
            raise ValueError(message)

        if syst_tag in self.is_applied.keys():
            if self.is_applied[syst_tag]:
                message =  "[WeightSystematic : apply] WeightSystematic: %s has already been applied to events set: %s" % (self.name, syst_tag)
                logger.exception(message)
                raise ValueError(message)

        weights = self.calculate_event_weights(events, central_only)

        self.summary[syst_tag] = {} 
        for variation, values in weights.items():
            mean = awkward.mean(values)
            std = awkward.std(values)
            self.summary[syst_tag][variation] = {
                    "branch" : self.name + "_" + variation,
                    "mean" : mean, 
                    "std" : std 
            }
            logger.debug("[WeightSystematic : apply] WeightSystematic: %s, variation: %s, independent collection: %s per-event weight has mean %.4f and std. dev. %.4f" % (self.name, variation, syst_tag, mean, std))

            if add_fields:
                awkward_utils.add_field(events, self.name + "_" + variation, values)

        if self.modify_central_weight:
            if not ("central" in weights.keys()):
                message = "[WeightSystematic : apply] WeightSystematic: %s, you selected to modify the central event weight, but a central value for this weight systematic was not found (variations are: %s)." % (self.name, str(weights.keys()))
                logger.exception(message)
                raise ValueError(message)

            weight = weights["central"]
            logger.debug("[WeightSystematic : apply] WeightSystematic: %s, modifying central weight '%s' by multiplying it by the per-event weight." % (self.name, CENTRAL_WEIGHT))

            if mask is None:
                mask = awkward.ones_like(weight)
            else:
                logger.debug("[WeightSystematic : apply] WeightSystematic: %s, a mask was provided, so we will modify central weight only where mask evaluates True." % (self.name)) 

            weight_masked = awkward.where(
                mask,
                weight,
                awkward.ones_like(weight)
            )

            if len(weight) > 0:
                self.summary[syst_tag]["central"]["frac"] = float(awkward.sum(mask)) / float(len(weight))
            else:
                self.summary[syst_tag]["central"]["frac"] = 0.

            if self.summary[syst_tag]["central"]["frac"] < 1.:
                logger.debug("[WeightSystematic : apply] WeightSystematic: %s, independent collection: %s, per-event weight is applied to %.2f percent of events." % (self.name, syst_tag, self.summary[syst_tag]["central"]["frac"] * 100.))

            self.summary[syst_tag]["central"]["before"] = awkward.mean(events[CENTRAL_WEIGHT])

            awkward_utils.add_field(
                    events = events,
                    name = CENTRAL_WEIGHT,
                    data = events[CENTRAL_WEIGHT] * weight_masked,
                    overwrite = True
            )
            logger.debug("event.central_weight%s"%events[CENTRAL_WEIGHT])

            self.summary[syst_tag]["central"]["after"] = awkward.mean(events[CENTRAL_WEIGHT])

            logger.debug("[WeightSystematic : apply] WeightSystematic: %s, independent collection: %s, mean value of central weight before/after applying: %.4f/%.4f" % (self.name, syst_tag, self.summary[syst_tag]["central"]["before"], self.summary[syst_tag]["central"]["after"]))
     
        self.is_applied[syst_tag] = True
        return events

    def calculate_event_weights(self, events, central_only):
        """
        Abstract function that is reimplemented for EventWeightSystematic and ObjectWeightSystematic.
        Translates the weight variations into a per-event weight.
            - In the case of an EventWeightSystematic, this is trivial.
            - In the case of an ObjectWeightSystematic, the per-object weights are multiplied together for the specified target collection to create the per-event weight. 

        :param events: events array to calculate per-event weights
        :type events: awkward.highlevel.Array
        :param central_only: whether to only calculate the central variation
        :type central_only: bool
        """
        raise NotImplementedError()
    

class ObjectWeightSystematic(WeightSystematic):
    """
    Class for calculating a weight systematic which is initially calculated for an object quantity, where there may be more than one object per event.
    The per-object weight is assumed to be translated to a per-event weight by multiplying the weights of all of the relevant objects for a given event.
    The objects that the per-object weight is calculated on need not be the same objects that the per-event weight is calculated from.

    :param input_collection: name of record that contains the objects for calculating the per-object weights
    :type input_collection: str
    :param target_collection: name of record that contains the objects for translating the per-object weights into per-event weights. May be passed as a string, e.g. "SelectedElectron" or as a tuple, in the case of nested records, e.g. ("Diphoton", "LeadPhoton")
    :type target_collection: str or tuple
    """
    def __init__(self, name, method, modify_central_weight, input_collection, branches = None, function = None, sample = None, target_collection = None):
        super(ObjectWeightSystematic, self).__init__(name, method, modify_central_weight, branches = branches, function = function, sample = sample)

        self.input_collection = input_collection
        if target_collection is None:
            self.target_collection = self.input_collection
        else:
            self.target_collection = target_collection

        # If user specified an input collection, update the branches to be in tuple format so it can be accessed directly from events as events[full_branch]
        if self.method == "from_branch" and self.input_collection is not None:
            for variation, branch in self.branches.items():
                if isinstance(branch, str):
                    full_branch = (self.input_collection, branch)
                elif isinstance(branch, tuple):
                    full_branch = (self.input_collection,) + branch
                self.branches[variation] = full_branch


    def calculate_variations(self, events, central_only):
        """
        Calls the base class WeightSystematic.calculate_variations method and also renames these to include the input collection name (so these can be added as fields for the relevant record). 
        """
        variations = super(ObjectWeightSystematic, self).calculate_variations(events, central_only)
        if self.method == "from_branch":
            return variations
        if self.method == "from_function":
            variations_renamed = {}
            for variation, values in variations.items():
                name = (self.input_collection, variation)
                variations_renamed[name] = values
            return variations_renamed


    def calculate_event_weights(self, events, central_only):
        weights = {}

        for variation, branch in self.branches.items():
            if central_only and not variation == "central":
                continue

            if self.input_collection == self.target_collection:
                target_branch = branch
            elif isinstance(self.target_collection, tuple):
                if isinstance(branch, tuple):
                    target_branch = self.target_collection + branch[1:]
                else:
                    target_branch = self.target_collection + (branch)
            else:
                target_branch = (self.target_collection,) + branch[1:]

            weights[variation] = awkward.prod(events[target_branch], axis = 1)
   
            name = "weight_" + self.name + "_" + variation
            awkward_utils.add_field(events, name, weights[variation], overwrite = True)
        logger.debug("add variation weight: %s"%weights)
        return weights

    def apply(self, events, syst_tag = NOMINAL_TAG, central_only = False, mask = None, add_fields = True):
        return super(ObjectWeightSystematic, self).apply(events, syst_tag, central_only, mask, add_fields = add_fields)


class EventWeightSystematic(WeightSystematic):
    """
    Class for calculating a weight systematic which is calculated from per-event quantities (as opposed to per-object quantities) 
    """

    def calculate_event_weights(self, events, central_only):
        weights = {}

        for variation, branch in self.branches.items():
            if central_only and not variation == "central":
                continue
            weights[variation] = events[branch]

        return weights

    def apply(self, events, syst_tag = NOMINAL_TAG, central_only = False, mask = None, add_fields = False):
        return super(EventWeightSystematic, self).apply(events, syst_tag, central_only, mask, add_fields = add_fields) 


class SystematicWithIndependentCollection(Systematic):
    """
    Abstract class for a systematic which modifies event content, i.e. a field in the original events array has its values modified. 
    Derives from higgs_dna.systematics.systematic.Systematic

    :param branch_modified: the field in the events array which will be modified by the up/down variations of this systematic
    :type branch_modified: str or tuple
    :param modify_nominal: whether to modify the nominal value of <branch_modified> 
    :type modify_nominal: bool 
    :param nominal_only: whether to calculate only the nominal correction (and not up/down variations)
    :type nominal_only: bool
    :param additive: for ICs from branches, whether the branches with variations should be added to the nominal branch
    :type additive: bool, defaults to False 
    """
    def __init__(self, name, method, branch_modified, modify_nominal = True, nominal_only = False, branches = None, function = None, sample = None, additive = False):
        super(SystematicWithIndependentCollection, self).__init__(name, method, branches, function, sample)

        self.branch_modified = branch_modified
        self.modify_nominal = modify_nominal
        self.nominal_only = nominal_only
        self.additive = additive

    def produce(self, events):
        """
        Produces copies of the original events array where the branch_modified field is updated with the corresponding variations.

        :param events: events array to use for calculating variations for this SystematicWithIndependentCollection
        :type events: awkward.highlevel.Array
        :return: dictionary of independent collections of the format "variation" : events_with_variation. For example: varying photon pt would result in { "up" : events_with_photon_pt_varied_up, "down" : ... }
        :rtype: dict 
        """
        if self.method == "from_branch":
            self.check_fields(events)
            independent_collections = self.produce_from_branch(events)

        elif self.method == "from_function":
            independent_collections = self.produce_from_function(events)


        if not isinstance(independent_collections, dict):
            message = "[SystematicWithIndependentCollection : produce] SystematicWithIndependentCollection: %s the 'produce' method should return a dictionary of the form { 'variation_name' : modified_events, ... }, not %s as we got." % (self.name, str(type(independent_collections)))
            logger.exception(message)
            raise TypeError(message)

        # Run various checks to make sure the independent collections dict is structured properly
        for variation, events_var in independent_collections.items():
            if not isinstance(variation, str):
                message = "[SystematicWithIndependentCollection : produce] SystematicWithIndependentCollection: %s, each key in the independent_collections dict should be a <str>, not %s as we got." % (self.name, str(type(variation)))
                logger.exception(message)
                raise TypeError(message)

            if not isinstance(events_var, awkward.highlevel.Array):
                message = "[SystematicWithIndependentCollection : produce] SystematicWithIndependentCollection: %s, each value in the independent_collections dict should be an <awkward.highlevel.Array>, not %s as we got." % (self.name, str(type(events_var)))
                logger.exception(message)
                raise TypeError(message)

            # Check that events_var has same structure as input array
            if not (len(events) == len(events_var)):
                message = "[SystematicWithIndependentCollection : produce] SystematicWithIndependentCollection: %s, each independent collection should have the same number of events as the input events array. We got len(events) = %d, len(independent_collection) = %d" % (self.name, len(events), len(events_var))
                logger.exception(message)
                raise ValueError(message)

            # Check that events_var has same fields as input array
            if not (set(events.fields) == set(events_var.fields)):
                message = "[SystematicWithIndependentCollection : produce] SystematicWithIndependentCollection: %s, each independent collection should have the same fields as the input events array. We got events.fields = %s, independent_collection.fields = %s" % (self.name, events.fields, events_var.fields)
                logger.exception(message)
                raise ValueError(message)

        self.is_produced = True
        return independent_collections


    def produce_from_branch(self, events):
        """
        Produce independent collections from a branch in the input ntuple.

        :param events: events array to use for calculating variations for this SystematicWithIndependentCollection
        :type events: awkward.highlevel.Array
        :return: dictionary of independent collections of the format "variation" : events_with_variation. For example: varying photon pt would result in { "up" : events_with_photon_pt_varied_up, "down" : ... }
        :rtype: dict 
        """
        independent_collections = {}
        for variation, branch in self.branches.items():
            if variation == NOMINAL_TAG and not self.modify_nominal:
                continue
            if not variation == NOMINAL_TAG and self.nominal_only:
                continue
            name = self.name + "_" + variation
            if self.additive:
                new_branch_content = events[self.branch_modified] + events[branch]
            else:
                new_branch_content = events[branch]
            independent_collections[variation] = awkward.with_field(
                    events,
                    new_branch_content,
                    self.branch_modified
            )

        return independent_collections

    
    def produce_from_function(self, events):
        """
        Produce independent collections from a function.

        :param events: events array to use for calculating variations for this SystematicWithIndependentCollection
        :type events: awkward.highlevel.Array
        :return: dictionary of independent collections of the format "variation" : events_with_variation. For example: varying photon pt would result in { "up" : events_with_photon_pt_varied_up, "down" : ... }
        :rtype: dict 
        """
        independent_collections = {}
        kwargs = self.get_function_kwargs(events)
        variations = getattr(
                self.function["module"],
                self.function["name"]
        )(**kwargs)

        for variation, branch in variations.items():
            if variation == NOMINAL_TAG and not self.modify_nominal:
                continue
            if not variation == NOMINAL_TAG and self.nominal_only:
                continue 
            name = self.name + "_" + variation
            independent_collections[variation] = awkward.with_field(
                    events,
                    branch,
                    self.branch_modified
            )

        return independent_collections

