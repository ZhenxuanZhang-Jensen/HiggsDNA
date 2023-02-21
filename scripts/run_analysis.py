import argparse
import json
import contextlib
import sys 

from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.utils.misc_utils import expand_path
from higgs_dna.analysis import AnalysisManager

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run taggers with specified configuration")

    # Optional arguments
    parser.add_argument(
        "--config",
        required=False,
        default=None,
        type=str,
        help="json config file to run analysis")

    parser.add_argument(
        "--log-level",
        required=False,
        default="INFO",
        type=str,
        help="Level of information printed by the logger")

    parser.add_argument(
        "--log-file",
        required=False,
        type=str,
        help="Name of the log file")

    parser.add_argument(
        "--output_dir",
        required=False,
        default="output",
        type=str,
        help="output directory to store output files")

    parser.add_argument(
        "--batch_system",
        required=False,
        default=None, # will default to "local" inside AnalysisManager
        type=str,
        help="batch system to run jobs on")

    parser.add_argument(
        "--merge_outputs",
        required=False,
        action="store_true",
        help="merge output files all into a single file")

    parser.add_argument(
        "--unretire_jobs",
        required=False,
        action="store_true",
        help="resubmit jobs that failed more than the max number of tries and were retired. Only applicable if you are re-running an existing analysis and you had jobs that failed multiple times (presumably due to corruptions/transient xrootd errors).")

    parser.add_argument(
        "--retire_jobs",
        required=False,
        action="store_true",
        help="retire all unfinished jobs and allow HiggsDNA to finish: calculate and apply scale1fb (for MC), merge outputs (if selected) and write summary json.")

    parser.add_argument(
        "--reconfigure_jobs",
        required=False,
        action="store_true",
        help="force the JobsManager to reconfigure each job. This option could be useful in the scenario that you want to update some detail in the job scripts (e.g. change the requested list of sites, change the requested memory, etc) but don't want to run a whole new analysis. Likely only useful for non-local job submission")

    parser.add_argument(
        "--short",
        required=False,
        action="store_true",
        help="just run 1 job for each sample/year to test workflow")

    parser.add_argument(
        "--years",
        required=False,
        default=None,
        type=str,
        help="csv list of years to process (overrides setting in config json")

    parser.add_argument(
        "--sample_list",
        required=False,
        default=None,
        type=str,
        help="csv list of samples to process (overrides setting in config json")

    parser.add_argument(
        "--fpo",
        required=False,
        default=None,
        type=int,
        help="number of input files per each job. Overrides any fpo specified for individual samples. If rerunning from a partial run of this script, does not override previous fpo.")

    parser.add_argument(
        "--n_cores",
        required=False,
        default=6,
        type=int,
        help="number of cores to use for running jobs in parallel. Only applicable if running locally.")

    parser.add_argument(
        "--use_xrdcp",
        required=False,
        default=False,
        type=bool,
        help="use xrdcp to copy to nanoAOD file to local or not ")

    return parser.parse_args()


def main(args):
    logger = setup_logger(args.log_level, args.log_file)

    if args.config is not None:
        with open(expand_path(args.config), "r") as f_in:
            args.config = json.load(f_in)

    logger.debug("Running HiggsDNA analysis with config:")

    args = {k:v for k,v in vars(args).items() if v is not None} # throw away None-value args
    analysis = AnalysisManager(**args)
    analysis.run()


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
