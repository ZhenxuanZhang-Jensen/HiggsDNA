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
                self.sum_weights = numpy.sum(runs["genEventSumw"].array())
                logger.debug("using genEventCount to check generator level eventNum")
                logger.debug("using genEventSumw to calculate Nanoaod events number")
            elif "genEventCount_" in runs.keys() and "genEventSumw_" in runs.keys():
                self.n_events = int(numpy.sum(runs["genEventCount_"].array()))
                self.sum_weights = numpy.sum(runs["genEventSumw_"].array())
                logger.debug("using genEventCount_ to check generator level eventNum")
                logger.debug("using genEventSumw_ to calculate Nanoaod events number")
            else: 
                self.n_events = 0
                self.sum_weights = 0 
        logger.debug("MC generate n_events =  %s"%self.n_events)
        logger.debug("NanoAod files N events = %s"%self.sum_weights)