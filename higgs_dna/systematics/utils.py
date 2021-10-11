import awkward
import numpy
import logging
logger = logging.getLogger(__name__)

def systematic_from_bins(bins, variables, central_only = False):
    """
    Function for calculating weight variations from bins which give the value and uncertainty.
    :param bins: dictionary containing variables used for binning and a list of the bins and their corresponding values
    :type bins: dict
    :param variables: dictionary of variables used for binning and the arrays that should be used to evaluate them
    :type variables: dict
    :param central_only: only compute the central variation (not up/down)
    :type central_only: bool
    :return: dictionary of variations in the format "var" : weight_var where "var" would be "up"/"down"/"central"
    :rtype: dict
    """

    # TODO: check bins and variables and make sure everything is passed properly

    # Initialize arrays of ones
    for variable, array in variables.items():
        weight = awkward.ones_like(array)
        if not central_only:
            weight_up = awkward.ones_like(array)
            weight_down = awkward.ones_like(array)
        break

    # Loop through bins and update central weight + up/down variations
    for bin in bins["bins"]:
        bin_cuts = []
        # Loop through variables which define bins
        for variable, contents in variables.items():
            cut_lower = contents >= bin[variable][0]
            cut_upper = contents < bin[variable][1]
            bin_cuts.append(cut_lower & cut_upper)

        # Create mask for selecting events in this bin
        for idx, cut in enumerate(bin_cuts):
            if idx == 0:
                bin_cut = cut
            else:
                bin_cut = bin_cut & cut

        # Update values of events in this bin
        weight = awkward.where(
                bin_cut,
                weight * bin["value"],
                weight
        )

        if not central_only:
            weight_up = awkward.where(
                    bin_cut,
                    weight_up * (bin["value"] + bin["uncertainty"]),
                    weight_up
            )

            weight_down = awkward.where(
                    bin_cut,
                    weight_down * (bin["value"] - bin["uncertainty"]),
                    weight_down
            )

    variations = { "central" : weight }
    if not central_only:
        variations["up"] =  weight_up
        variations["down"] =  weight_down

    return variations


def ic_systematic_from_bins(bins, variables, branch, nominal_only = False, modify_nominal = False, mask = None):
    """
    Function for calculating systematics with independent collections (e.g. correcting and/or varying an energy scale)
    from bins which give the value and uncertainty of the correction.
    NOTE: corrections and uncertainties are assumed to be fractional. I.e. a correction of 1.01 and an uncertainty of 0.05 would give an object with 10 GeV a corrected value of 10.1 GeV and up/down variations of 10 * (1.01 +/- 0.05) GeV

    :param bins: dictionary containing variables used for binning and a list of the bins and their corresponding values
    :type bins: dict
    :param variables: dictionary of variables used for binning and the arrays that should be used to evaluate them
    :type variables: dict
    :param nominal_only: optional, only compute the nominal variation (not up/down)
    :type nominal_only: bool, defaults to False
    :param modify_nominal: optional, modify the nominal value of the branch to be modified (and its up/down variations)
    :type modify_nominal: bool, defaults to False
    :param mask: optional, a boolean mask indicating the subset of events/objects for which the correction and variations should be applied
    :type mask: awkward.highlevel.Array, defaults to None
    :return: dictionary of variations in the format "var" : branch_var where "var" would be "up"/"down"/"nominal"
    :rtype: dict
    """
    variations = {}
    if modify_nominal:
        variations["nominal"] = awkward.copy(branch)

    if not nominal_only:
        variations["up"] = awkward.copy(branch)
        variations["down"] = awkward.copy(branch)

    scale = awkward.ones_like(branch)
    for bin in bins["bins"]:
        bin_cuts = []
        # Loop through variables which define bins
        for variable, contents in variables.items():
            cut_lower = contents >= bin[variable][0]
            cut_upper = contents < bin[variable][1]
            bin_cuts.append(cut_lower & cut_upper)

        # Create mask for selecting events in this bin
        for idx, cut in enumerate(bin_cuts):
            if idx == 0:
                bin_cut = cut
            else:
                bin_cut = bin_cut & cut

        if mask is not None:
            bin_cut = bin_cut & mask

        if modify_nominal:
            scale = awkward.where(
                    bin_cut,
                    bin["value"],
                    scale
            )

            variations["nominal"] = awkward.where(
                    bin_cut,
                    branch * scale,
                    variations["nominal"]
            )

        if not nominal_only:
            variations["up"] = awkward.where(
                    bin_cut,
                    branch * (scale + bin["uncertainty"]),
                    variations["up"]
            )
            variations["down"] = awkward.where(
                    bin_cut,
                    branch * (scale - bin["uncertainty"]),
                    variations["down"]
            )

    return variations
