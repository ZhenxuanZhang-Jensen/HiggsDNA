import os
import sys
import json
import copy

import logging
logger = logging.getLogger(__name__)

def load_config(config):
    """
    Load a dictionary or path to json file into a dictionary.

    :param config: dictionary or path to json file with config options
    :type config: dict/str
    :return: config options as dictionary
    :rtype: dict
    """
    if isinstance(config, dict):
        return config

    elif isinstance(config, str):
        if not config.endswith(".json"):
            logger.warning("[load_config] : You specified a string, meaning we load config options from a json file, but %s does not look like a json file!" % config)

        with open(expand_path(config), "r") as f_in:
            options = json.load(f_in)
            return options

    else:
        message = "[load_config] : You tried to load a config file with '%s' which has type <%s>, but only loading from json file (<str>) or <dict> are supported." % (str(config), str(type(config)))
        logger.exception(message)
        raise TypeError(message)


def get_module(module_name):
    """
    Check if module has been imported.
    Raises ImportError if not, otherwise returns the module.

    :param module_name: name of module in format "higgs_dna.module", e.g. "higgs_dna.systematics.lepton_systematics"
    :type module_name: str
    :return: corresponding module in higgs_dna
    :rtype: module
    """
    if module_name not in sys.modules:
        logger.exception("[misc_utils : get_module] Module '%s' has not been imported. Please add 'import %s' to your analysis script" % (module_name, module_name))
        raise ImportError()
        
    return sys.modules[module_name]


def expand_path(relative_path):
    """
    Convert a relative path (assumed to be the path under HiggsDNA/) into an absolute path.

    :param relative_path: path under HiggsDNA/
    :type relative_path: str
    :return: absolute path
    :rtype: str
    """

    dir = os.path.dirname(__file__)
    subdirs = dir.split("/")

    base_path = ""
    for subdir in subdirs:
        if subdir == "higgs_dna":
            break
        base_path += subdir + "/"
        if subdir == "HiggsDNA":
            break

    return base_path + relative_path


def update_dict(original, new):
    """
    Update nested dictionary (dictionary possibly containing dictionaries)
    If a field is present in new and original, take the value from new.
    If a field is present in new but not original, insert this field 


    :param original: source dictionary
    :type original: dict
    :param new: dictionary to take new values from
    :type new: dict
    :return: updated dictionary
    :rtype: dict
    """

    updated = copy.deepcopy(original)

    for key, value in original.items():
        if key in new.keys():
            if isinstance(value, dict):
                updated[key] = update_dict(value, new[key])
            else:
                updated[key] = new[key]

    return updated

from higgs_dna.utils.metis_utils import do_cmd
def check_proxy():
    """
    Check if a valid grid proxy exists.

    :return: path to proxy if it exists, otherwise None
    :rtype: str
    """

    proxy = None
    bad_proxy = False
    proxy_info = do_cmd("voms-proxy-info").split("\n")
    for line in proxy_info:
        if "path" in line:
            proxy = line.split(":")[-1].strip()

        if "timeleft" in line:
            time_left = int(line.replace("timeleft", "").replace(":", "").strip())
            if not time_left > 0:
                bad_proxy = True

        if "Couldn't find a valid proxy." in line:
            bad_proxy = True

    if proxy is None or bad_proxy:
        logger.warning("[misc_utils : check_proxy] We were not able to find grid proxy or proxy was found to be expired. Output of 'voms-proxy-info': %s" % (str(proxy_info)))
        return None

    else:
        return proxy
