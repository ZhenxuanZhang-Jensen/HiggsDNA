import os
import glob
import json
import copy
from tqdm import tqdm 

import logging
logger = logging.getLogger(__name__)

from higgs_dna.samples.sample import Sample
from higgs_dna.samples.file import File
from higgs_dna.utils.misc_utils import load_config 
from higgs_dna.utils import metis_utils, misc_utils

class SampleManager():
    """

    """
    def __init__(self, sample_list, years, catalog = "metadata/samples_catalog.json"):
        self.sample_list = sample_list
        if not isinstance(self.sample_list, list):
            self.sample_list = [self.sample_list]

        self.years = years
        if not isinstance(self.years, list):
            self.years = [self.years]

        self.catalog_name = misc_utils.expand_path(catalog)
        self.catalog = load_config(catalog)

        self.samples = {}
        self.loaded_samples = False
        self.data = None
        self.process_id_map = {}


    def get_samples(self):
        if self.loaded_samples:
            return self.data

        samples = []
        # Loop through each sample
        logger.info("[SampleManager : get_samples] Fetching input files for %d samples." % (len(self.sample_list)))
        for s_idx, sample in enumerate(tqdm(self.sample_list)):
            if sample not in self.catalog.keys():
                logger.exception("[SampleManager : get_samples] Could not find sample '%s' in samples catalog." % (sample))
                raise ValueError()

            self.samples[sample] = {}
            info = self.catalog[sample]
            # Loop through years and create a separate Sample object for each year
            for year in self.years:
                # Get xs and bf info
                if "xs" in info.keys():
                    if isinstance(info["xs"], dict): # different xs for different years
                        xs = info["xs"][year]
                    else:
                        xs = info["xs"]

                    if "bf" in info.keys():
                        if isinstance(info["bf"], dict): # different bf for different years
                            bf = info["bf"][year]
                        else:
                            bf = info["bf"]
                    else:
                        bf = 1.

                else:
                    xs = None
                    bf = None

                if xs is not None:
                    logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', found xs of %.6f pb and bf of %.6f" % (sample, year, xs, bf))
                    is_data = False
                else:
                    logger.debug("[SampleManager : get_samples] For sample '%s', year '%s' no 'xs' info was found, treating this as data." % (sample, year))
                    is_data = True
                
                # Get input files
                files = []
                if year not in info["files"].keys():
                    logger.exception("[SampleManager : get_samples] Could not find any information about 'files' in samples catalog for sample '%s', year '%s'." % (sample, year))
                    raise ValueError()

                self.samples[sample][year] = {}

                # Multiple ways the user can specify files:
                #   0. Directory with wildcards
                #   1. Hard-coded list (local or xrd format)
                #   2a. Single directory (local or xrd format)
                #   2b. List of directories (local or xrd format)
                #   3a/b. DAS dataset name 

                grabbed_files = False

                if isinstance(info["files"][year], str) and "*" in info["files"][year]: # directory with wildcards
                    files = self.get_files_from_wildcard(info["files"][year], is_data)
                    grabbed_files = True

                # Is this a list? Could be a list of hard-coded files (Option 1) or list of dirs (Option 2b) 
                elif isinstance(info["files"][year], list):
                    if info["files"][year][0].endswith(".root"): # Option 1
                        logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', getting files from hard-coded list." % (sample, year))
                        files = [File(name = x, is_data = is_data) for x in info["files"][year]]
                        grabbed_files = True


                if not grabbed_files:
                    if isinstance(info["files"][year], str): # Option 2a/3a 
                        info["files"][year] = [info["files"][year]] # recast this as Option 2b/3b format

                    for path in info["files"][year]:
                        if path.startswith("root://"): # access via xrd
                            logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', we interpreted the specified files '%s' as a directory to be accessed via xrd" % (sample, year, path))

                            # Check for valid proxy
                            proxy = misc_utils.check_proxy()
                            if proxy is None:
                                logger.exception("[CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found to be expired. Since you are accessing files through xrd, a valid proxy is necessary.")
                                raise RuntimeError()
                            files_dir = self.get_files_from_xrd(path, is_data)

                        elif os.path.exists(path): # Option 1: local
                            logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', we interpreted the specified files '%s' as a local directory, whose files will be grabbed with <glob>." % (sample, year, path))
                            files_dir = self.get_files_from_local_dir(path, is_data)

                        elif path.endswith("NANOAOD") or path.endswith("NANOAODSIM"): # Option 3: DAS
                            logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', we interpreted the specified files '%s' as a DAS dataset, whose files will be grabbed with <dasgoclient>." % (sample, year, path))

                            # Check for valid proxy
                            proxy = misc_utils.check_proxy()
                            if proxy is None:
                                logger.exception("[CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found to be expired. Since you are accessing files through dasgoclient, a valid proxy is necessary.")
                                raise RuntimeError()
                            files_dir = self.get_files_from_dasgoclient(path, is_data)

                        files += files_dir
                    grabbed_files = True

                files = sorted(files, key = lambda x : x.name)

                logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', found %d input files:" % (sample, year, len(files)))
                if len(files) < 50: # don't print out if more than 50
                    for file in files:
                        logger.debug("\t %s" % file.name)

                # Check for sample-specific systematics
                if "systematics" in info.keys():
                    if year in info["systematics"].keys(): # year-specific syst
                        systematics = info["systematics"][year]
                    else: # same systs for all years
                        systematics = info["systematics"]
                    logger.debug("[SampleManager : get_samples] For sample '%s', year '%s' adding sample-specific systematics: %s" % (sample, year, str(systematics)))
                else:
                    systematics = None

                self.samples[sample][year] = {
                        "files" : files,
                        "xs" : xs,
                        "bf" : bf,
                }
                if systematics is not None:
                    self.samples[sample][year]["systematics"] = systematics

                # Check if fpo is specified
                if "fpo" in info.keys():
                    fpo = info["fpo"]
                else:
                    fpo = None

                samples.append(
                        Sample(
                            process = sample,
                            year = year,
                            files = files,
                            xs = xs,
                            bf = bf,
                            process_id = s_idx,
                            fpo = fpo,
                            systematics = systematics
                        )
                )
                self.process_id_map[sample] = s_idx

        self.data = samples
        self.loaded_samples = True

        self.update_catalog(samples)

        return samples


    def update_catalog(self, samples):
        """
        Rewrite the input catalog with the full list of files written out explicitly
        This way it is not necessary to rerun the glob/xrootd/dasgoclient commands and we save some time.
        """

        if "_sample_manager_full.json" in self.catalog_name: # we already made the full sample list on a previous run
            return

        catalog_full = copy.deepcopy(self.catalog)

        for sample in samples:
            catalog_full[sample.process]["files"][sample.year] = [x.name for x in sample.files]

        with open(self.catalog_name.replace(".json", "_sample_manager_full.json"), "w") as f_out:
            json.dump(catalog_full, f_out, indent = 4)


    def get_files_from_dasgoclient(self, sample, is_data):
        """
        Inspired by function from Nick Amin in ProjectMetis:
        https://github.com/aminnj/ProjectMetis/blob/f9e71556cb84496731fa71dcab2dfc82b6e3022f/metis/Sample.py#L224-L236

        Get list of files, along with number of events and size in GB for each.

        :param sample: DBS sample in "three-slash" format, e.g. "/EGamma/Run2018A.../NANOAOD"
        :type sample: str
        :return: dictionary of files and metadata
        :rtype: dict
        """
        results = {}

        #cmd = "dasgoclient -query 'file dataset={}' -json".format(sample)
        cmd = "/cvmfs/cms.cern.ch/common/dasgoclient -query 'file dataset={}' -json".format(sample)
        query = json.loads(
                metis_utils.do_cmd(cmd)
        )

        files = []
        for j in query:
            f = j["file"][0]
            files.append(
                    File(
                        name = "root://cmsxrootd.fnal.gov/" + f["name"],
                        is_data = is_data,
                        n_events = f["nevents"],
                        size_gb = round(f["size"]*1e-9,2) 
                    )
            )

        return files


    def get_files_from_wildcard(self, path, is_data):
        """

        """
        files_dir = glob.glob(path)
        logger.debug("[SampleManager : get_files_from_wildcard] Found %d files from specified wildcard '%s'." % (len(files_dir), path))
        files = []
        for f in files_dir:
            if ".root" not in f:
                logger.debug("[SampleManager : get_files_from_wildcard] File '%s' was grabbed under your specified wildcard '%s' but is not a .root file, so we are skipping it." % (f, path))
                continue

            files.append(
                    File(
                        name = f,
                        is_data = is_data,
                        size_gb = round(os.path.getsize(f)*1e-9,2)
                    )
            )

        return files


    def get_files_from_local_dir(self, directory, is_data):
        """

        """
        #files_dir = glob.glob(info["files"][year] + "/*.root")
        files_dir = glob.glob(directory + "/*.root")
        logger.debug("[SampleManager : get_files_from_local_dir] Found %d files in dir '%s'." % (len(files_dir), directory))
        files = []
        for f in files_dir:
            files.append(
                    File(
                        name = f,
                        is_data = is_data,
                        size_gb = round(os.path.getsize(f)*1e-9,2)
                    )
            )

        return files


    def get_files_from_xrd(self, directory, is_data):
        """
        Get list of files from an xrootd directory with xrdfs ls command"

        :param directory: directory (in xrootd format) to glob files from
        :type directory: str
        :return: list of all root files from directory
        :rtype: list of str
        """
        files = []

        idx = directory.find("//store") + 1
        redirector = directory[:idx]
        dir = directory[idx:]

        command = "xrdfs %s ls %s" % (redirector, dir)

        logger.debug("[SampleManager : get_files_from_xrd] We will find files for dir '%s' with the command: \n %s" % (directory, command))

        contents = os.popen(command).read().split("\n")
        for x in contents:
            if x.endswith(".root"):
                files.append(File(name = redirector + x, is_data = is_data))

        logger.debug("[SampleManager : get_files_from_xrd] Found %d files in dir '%s'." % (len(files), directory))
            
        return files
