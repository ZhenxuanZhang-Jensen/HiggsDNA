import awkward
import vector

vector.register_awkward()

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


def mass_veto(objects1, objects2, mass_range):
    """
    Select objects from objects1 which have invariant mass with all objects in objects2 outside `mass_range`.
    For example, `mass_veto(electrons, photons, [85., 95.])` selects all electrons which have m(e,gamma) outside the Z peak.

    :param objects1: objects which are required to be at least min_dr away from all objects in objects2 
    :type objects1: awkward.highlevel.Array
    :param objects2: objects which are all objects in objects1 must be at leats min_dr away from
    :type objects2: awkward.highlevel.Array
    :param mass_range: mass range to veto
    :type mass_range: list of float 
    :return: boolean array of objects in objects1 which pass mass veto requirement
    :rtype: awkward.highlevel.Array
    """
    if awkward.count(objects1) == 0 or awkward.count(objects2) == 0:
        return objects1.pt < 0. 

    if not isinstance(objects1, vector.Vector4D):
        objects1 = awkward.Array(objects1, with_name = "Momentum4D")
    if not isinstance(objects2, vector.Vector4D):
        objects2 = awkward.Array(objects2, with_name = "Momentum4D")

    obj1 = awkward.unflatten(objects1, counts = 1, axis = -1) # shape [n_events, n_obj, 1]
    obj2 = awkward.unflatten(objects2, counts = 1, axis = 0) # shape [n_events, 1, n_obj]

    mass = (obj1 + obj2).mass
    
    selection = (awkward.all(mass < mass_range[0], axis = -1)) | (awkward.all(mass > mass_range[1], axis = -1)) 
    return selection


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
    #obj1:jet obj2:muon     
    # if awkward.count(objects1) == 0 or awkward.count(objects2) == 0:
    #     return objects1.pt < 0. 
    if awkward.count(objects1) == 0: 
        return objects1.pt < 0. 
    if awkward.count(objects2) == 0:
        return objects1.pt > 0.

    if not isinstance(objects1, vector.Vector4D):
        objects1 = awkward.Array(objects1, with_name = "Momentum4D")
    if not isinstance(objects2, vector.Vector4D):
        objects2 = awkward.Array(objects2, with_name = "Momentum4D")

    obj1 = awkward.unflatten(objects1, counts = 1, axis = -1) # shape [n_events, n_obj, 1]
    obj2 = awkward.unflatten(objects2, counts = 1, axis = 0) # shape [n_events, 1, n_obj]

    dR = obj1.deltaR(obj2) # shape [n_events, n_obj1, n_obj2]

    selection = awkward.all(dR >= min_dr, axis = -1)
    return selection


# This is the same deltaR function as above, but calculated with an explicit loop through both sets of objects (C++ style)
# and compiled with numba. This is just to illustrate how to do things both ways.
def delta_R_numba(objects1, objects2, min_dr):
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


def get_closest(base_objects, target_objects, max_dr = 999999., max_pt_diff = 999999., name = ""):
    """
    For each object in `base_objects` find the closest object from each object in `target_objects`. 
    Cuts on the maximum dR and maximum pt difference are also available. These might be useful for e.g. finding electrons which are al
so in the photons collection, where you don't simply want the nearest photon, but a photon within a certain dR and within a certain pT
 difference.
    If there is more than one possible match, the match with the smallest dR is taken.
    """
    if not isinstance(base_objects, vector.Vector4D):
        base_objects = awkward.Array(base_objects, with_name = "Momentum4D")
    if not isinstance(target_objects, vector.Vector4D):
        target_objects = awkward.Array(target_objects, with_name = "Momentum4D")

    # If either set of objects is flat, make them jagged for consistency
    if base_objects.ndim == 1:
        base_objects = awkward.unflatten(base_objects, counts=1, axis=-1)
    if target_objects.ndim == 1:
        target_objects = awkward.unflatten(target_objects, counts=1, axis=-1)

    base_unflat = awkward.unflatten(base_objects, counts=1, axis=-1) # shape [n_events, n_base_objects, 1]
    target_unflat = awkward.unflatten(target_objects, counts=1, axis=0) # shape [n_events, 1, n_target_objects]

    dR = base_unflat.deltaR(target_unflat) # shape [n_events, n_base_objects, n_target_objects]
    pt_diff = base_unflat.pt - target_unflat.pt

    # Find the best target match for each base particle (if any)
    target_match = awkward.unflatten(target_objects, counts = awkward.num(base_objects), axis = 0)# shape [n_events, n_base_objects, n_target_objects]
    # select only those with dR < max_dr and pt difference < max_pt_diff
    target_match = target_match[(dR < max_dr) & (abs(pt_diff) < max_pt_diff)]

    base_best_target = awkward.firsts(
        target_match[awkward.argsort(target_match.deltaR(base_unflat), ascending=True, axis=-1)],
        axis=-1
    )

    base_objects["%spt_diff" % name] = abs(base_objects.pt - base_best_target.pt)
    base_objects["%sdR" % name] = base_objects.deltaR(base_best_target)

    return base_objects

