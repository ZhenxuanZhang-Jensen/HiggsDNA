import awkward
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

def abs_cos_helicity(p1, p2):
    """
    Compute abs(cos) of the helicity angle of the decay X->p1p2,
    with helicity angle defined as the angle between p1 and X in the rest frame of X.

    Note that abs_cos_helicity(p1,p2) == abs_cos_helicity(p2,p1).

    This follows the helicity angle definition of the CMS Run 2 ttH (H->gg) analysis, HIG-19-013.

    :param p1: four vector of p1
    :type p1: vector.Vector4D or awkward.Array with name "Momentum4D" 
    :param p2: four vector of p2
    :type p2: vector.Vector4D or awkward.Array with name "Momentum4D"
    :returns: absolute value of the cos theta
    :rtype: float or awkward.Array of floats 
    """
    parent = p1 + p2

    # Boost p1 to CM frame of parent
    p1_in_CM_parent = p1.boost_beta3(-parent.to_beta3()) # equivalent to p1.boostCM_of_p4(parent) in vector 0.8.5 and later

    p1_in_CM_parent3d = p1_in_CM_parent.to_Vector3D()
    parent3d = parent.to_Vector3D()

    cos_theta = p1_in_CM_parent3d.dot(parent3d) / (p1_in_CM_parent3d.mag * parent3d.mag)
    return abs(cos_theta)


def abs_cos_theta_parentCM(p1, p2):
    """
    Compute abs(cos(theta)) of p1 in the rest frame of the parent of p1 and p2 (p1 + p2)

    Note that abs_cos_theta_parentCM(p1,p2) == abs_cos_theta_parentCM(p2,p1)

    This follows the helicity angle definition of the CMS Run 2 HH->ggbb analysis, HIG-19-018.

    :param p1: four vector of p1
    :type p1: vector.Vector4D or awkward.Array with name "Momentum4D" 
    :param p2: four vector of p2
    :type p2: vector.Vector4D or awkward.Array with name "Momentum4D"
    :returns: absolute value of the cos theta
    :rtype: float or awkward.Array of floats 
    """
    parent = p1 + p2

    # Boost p1 to CM frame of parent
    p1_in_CM_parent = p1.boost_beta3(-parent.to_beta3()) # equivalent to p1.boostCM_of_p4(parent) in vector 0.8.5 and later
    return abs(p1_in_CM_parent.costheta)
