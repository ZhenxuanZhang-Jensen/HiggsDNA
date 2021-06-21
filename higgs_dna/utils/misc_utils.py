import json

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

        with open(config, "r") as f_in:
            options = json.load(f_in)
            return options

    else:
        message = "[load_config] : You tried to load a config file with '%s' which has type <%s>, but only loading from json file (<str>) or <dict> are supported." % (str(config), str(type(config)))
        logger.exception(message)
        raise TypeError(message)
            


