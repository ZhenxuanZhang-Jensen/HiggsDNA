import awkward
import vector
import numpy
import numba

from higgs_dna.utils import awkward_utils

import logging
logger = logging.getLogger(__name__)

def select_objects(objects, cuts = {}, clean = {}, name = "none", tagger = None):
    """

    """

    tagger_name = "none" if tagger is None else tagger.name

    if cuts:
        logger.debug("[select_objects] : Tagger '%s', selecting objects '%s', with the following requirements:" % (tagger_name, name))
        for cut, value in cuts.items():
            logger.debug("\t '%s' : %s" % (cut, str(value)))

    if clean:
        logger.debug("[select_objects] : Tagger '%s', cross-cleaning objects '%s' with respect to the following sets of objects:" % (tagger_name, name))
        for other_objects, info in clean.items():
            logger.debug("\t '%s', dR > %.2f" % (other_objects, info["min_dr"]))

    cut_names = []
    cut_results = []
    for cut, value in cuts.items():
        cut_ = None
        if cut == "pt":
            cut_ = objects.pt > value
            cut_names.append("pt > %.4f" % value)
        if cut in ["eta", "dxy", "dz"]:
            cut_ = abs(objects[cut]) < value
            cut_names.append("|%s| < %.4f" % (cut, value))
        if cut in ["pfRelIso03_all", "pfRelIso03_chg"]:
            cut_ = objects[cut] < value
            cut_names.append("%s < %.4f" % (cut, value))

        if cut_ is not None:
            cut_results.append(cut_)

    for other_objects, info in clean.items():
        cut_ = delta_R(objects, info["objects"], info["min_dr"])
        cut_names.append("dR with '%s' > %.2f" % (other_objects, info["min_dr"]))
        cut_results.append(cut_)

    all_cuts = objects.pt > 0
    for cut in cut_results:
        all_cuts = (all_cuts) & cut

    if tagger is not None:
        tagger.register_cuts(
                names = cut_names,
                results = cut_results,
                cut_type = name
        )

    return all_cuts


def delta_R(objects1, objects2, min_dr):
    """
    Select objects from objects1 which are at least min_dr away from all objects in objects2.

    :param objects1: objects which are required to be at least min_dr away from all objects in objects2 
    :type objects1: awkward.highlevel.Array
    :param objects2: objects which are all objects in objects1 must be at leats min_dr away from
    :type objects2: awkward.highlevel.Array
    :param min_dr: minimum delta R between objects
    :type min_dr: float
    :return: boolean array of objects in objects1 which pass delta_R requirement
    :rtype: awkward.highlevel.Array
    """
    n_objects1 = awkward.num(objects1)
    n_objects2 = awkward.num(objects2)

    offsets, contents = compute_delta_R(
            objects1, n_objects1,
            objects2, n_objects2,
            min_dr
    )

    selection = awkward_utils.construct_jagged_array(offsets, contents)

    return selection

@numba.njit
def compute_delta_R(objects1, n_objects1, objects2, n_objects2, min_dr):
    """
    :param objects1: objects which are required to be at least min_dr away from all objects in objects2 
    :type objects1: awkward.highlevel.Array
    :param n_objects1: number of objects in objects1 for each event
    :type n_objects1: awkward.highlevel.Array
    :param objects2: objects which are all objects in objects1 must be at leats min_dr away from
    :type objects2: awkward.highlevel.Array
    :param n_objects2: number of objects in objects2 for each event
    :type n_objects2: awkward.highlevel.Array
    :param min_dr: minimum delta R between objects
    :type min_dr: float
    :return: offsets corresponding to the indices which subdivide the contents into lists, contents of the underlying data for all lists
    :rtype: numpy.array, numpy.array
    """
    n_events = len(objects1)

    offsets = numpy.zeros(n_events + 1, numpy.int64)
    contents = []

    for i in range(n_events):
        offsets[i+1] = offsets[i] + n_objects1[i]
        for j in range(n_objects1[i]):
            contents.append(True)
            offset_idx = offsets[i] + j
            for k in range(n_objects2[i]):
                if not contents[offset_idx]:
                    continue
                obj1 = vector.obj(
                        pt = objects1[i][j].pt,
                        eta = objects1[i][j].eta,
                        phi = objects1[i][j].phi,
                        mass = objects1[i][j].mass
                )
                obj2 = vector.obj(
                        pt = objects2[i][k].pt,
                        eta = objects2[i][k].eta,
                        phi = objects2[i][k].phi,
                        mass = objects2[i][k].mass 
                )
                dR = obj1.deltaR(obj2)
                if dR < min_dr:
                    contents[offset_idx] = False

    return offsets, numpy.array(contents)

