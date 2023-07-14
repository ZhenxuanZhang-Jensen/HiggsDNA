import os
import uproot
import numpy

import logging
logger = logging.getLogger(__name__)

class File():
    """

    """
    def __init__(self, name, status = None, basepath = None, size_gb = None, n_events = None, sum_weights = None, is_data = False):
        self.name = name
        self.status = status
        self.basepath = basepath
        self.size_gb = size_gb
        self.n_events = n_events
        self.sum_weights = sum_weights
        self.is_data = is_data

        if self.basepath is not None:
            self.name = os.path.join(self.basepath, self.name)


    def calculate(self, treename = "Runs"):
        """
        Calculate number of events and sum of weights
        """
        if self.is_data:
            return

        with uproot.open(self.name) as f:
            runs = f["Runs"]
            if "genEventCount" in runs.keys() and "genEventSumw" in runs.keys():
                self.n_events = int(numpy.sum(runs["genEventCount"].array()))
                if "NMSSM" in name:
                	self.sum_weights = self.n_events
                else:
                	self.sum_weights = numpy.sum(runs["genEventSumw"].array())
            elif "genEventCount_" in runs.keys() and "genEventSumw_" in runs.keys():
                self.n_events = int(numpy.sum(runs["genEventCount_"].array()))
                if "NMSSM" in name:
                	self.sum_weights = self.n_events
                else:
                	self.sum_weights = numpy.sum(runs["genEventSumw_"].array())
            else: 
                self.n_events = 0
                self.sum_weights = 0 
