import logging

import awkward
import numpy
import vector
from higgs_dna.selections import (fatjet_selections, jet_selections,
                                  lepton_selections,gen_selections)
from higgs_dna.taggers.tagger import NOMINAL_TAG, Tagger
from higgs_dna.utils import awkward_utils, misc_utils
from matplotlib.pyplot import jet

vector.register_awkward()

logger = logging.getLogger(__name__)


DUMMY_VALUE = -999.

DEFAULT_OPTIONS = {

}


class GEN_Preselection(Tagger):
    """
    HHWW Preselection tagger for tutorial
    """

    def __init__(self, name, options={}, is_data=None, year=None):
        super(GEN_Preselection, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                original=DEFAULT_OPTIONS,
                new=options
            )

    def calculate_selection(self, events):
        # data will not select gen level infos
        # Gen selection
        # if not self.is_data:    
        gen_part = awkward.Array(events.GenPart,with_name="Momentum4D")
        logger.debug(" debug before gen selection :")        
        W1_candi,W2_candi,H_candi = gen_selections.select_ww_to_qqqq(gen_part)
        logger.debug(" debug after gen selection :")        
        # logger.debug(" the type of the W1candi is %s :"%type(W1_candi))        
        # logger.debug(" W1candi:",W1_candi)        
        W1_candi4D = awkward.zip(
            {
            "pt": numpy.array(W1_candi)[:,0],
            "eta": numpy.array(W1_candi)[:,1],
            "phi": numpy.array(W1_candi)[:,2],
            "mass": numpy.array(W1_candi)[:,3],
            }, 
        with_name="Momentum4D")
        W2_candi4D = awkward.zip(
            {
            "pt": numpy.array(W2_candi)[:,0],
            "eta": numpy.array(W2_candi)[:,1],
            "phi": numpy.array(W2_candi)[:,2],
            "mass": numpy.array(W2_candi)[:,3],
            }, 
        with_name="Momentum4D")
        H_candi4D = awkward.zip(
            {
            "pt": numpy.array(H_candi)[:,0],
            "eta": numpy.array(H_candi)[:,1],
            "phi": numpy.array(H_candi)[:,2],
            "mass": numpy.array(H_candi)[:,3],
            }, 
        with_name="Momentum4D")
        unflatten_W1_candi4D = awkward.unflatten(W1_candi4D,1)
        unflatten_W2_candi4D = awkward.unflatten(W2_candi4D,1)
        unflatten_H_candi4D = awkward.unflatten(H_candi4D,1)


        
        return events
