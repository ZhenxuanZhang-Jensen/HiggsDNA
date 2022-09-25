"""
Various utility functions for sample management and job submission.
Much of the code is taken from:
    ProjectMetis metis/Utils.py: https://github.com/aminnj/ProjectMetis/blob/f9e71556cb84496731fa71dcab2dfc82b6e3022f/metis/Utils.py
    Author: Nick Amin
"""

import subprocess
import datetime

import logging
logger = logging.getLogger(__name__)

def do_cmd(cmd, returnStatus=False, dryRun=False):
    """

    """
    if dryRun:
        print("dry run: {}".format(cmd))
        status, out = 1, ""
    else:
        status, out = subprocess.getstatusoutput(cmd)
    if returnStatus:
        return status, out
    else:
        return out

def do_cmd_timeout(cmd, timeout, loop = False):
    """
    Run a command ``cmd`` and kill the process if it does not finish before timeout.

    :param loop: whether to indefinitely loop through this until it finishes under the timeout
    :type loop: bool
    :returns: True if completed, False if timed out
    :rtype: bool
    """
    logger.debug("[misc_utils : do_cmd_timeout] Running command '%s' with timeout of %.0f seconds with indefinite looping until completion set to '%s'." % (cmd, timeout, loop))

    p = subprocess.Popen(cmd.split(" "))
    try:
        p.wait(timeout)
        return True
    except:
        p.kill()
        if loop:
            return do_cmd_timeout(cmd, timeout, loop)
        else:
            return False

def get_proxy_file():
    return "/tmp/x509up_u{0}".format(os.getuid())

def get_timestamp():
    # return current time as a unix timestamp
    return int(datetime.datetime.now().strftime("%s"))

def from_timestamp(timestamp):
    # return datetime object from unix timestamp
    return datetime.datetime.fromtimestamp(int(timestamp))

def timedelta_to_human(td):
    if td.days >= 2:
        return "{} days".format(td.days)
    else:
        return "{} hours".format(int(td.total_seconds()//3600))

def num_to_ordinal_string(n):
    # https://stackoverflow.com/questions/3644417/python-format-datetime-with-st-nd-rd-th-english-ordinal-suffix-like
    return str(n)+("th" if 4<=n%100<=20 else {1:"st",2:"nd",3:"rd"}.get(n%10, "th"))


