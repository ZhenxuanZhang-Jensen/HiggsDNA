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
