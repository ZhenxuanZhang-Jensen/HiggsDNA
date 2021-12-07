import argparse
import json

from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.utils.misc_utils import expand_path
from higgs_dna.analysis import AnalysisManager

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run taggers with specified configuration")

    # Required arguments
    parser.add_argument(
        "--config",
        required=True,
        type=str,
        help="json config file to run analysis")

    # Optional arguments
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
        default="local",
        type=str,
        help="batch system to run jobs on")

    parser.add_argument(
        "--merge_outputs",
        required=False,
        action="store_true",
        help="merge output files all into a single file")

    parser.add_argument(
        "--resubmit_retired",
        required=False,
        action="store_true",
        help="resubmit jobs that failed more than the max number of tries and were retired. Only applicable if you are re-running an existing analysis and you had jobs that failed multiple times (presumably due to corruptions/transient xrootd errors).")

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

    return parser.parse_args()


def main(args):
    logger = setup_logger(args.log_level, args.log_file)

    with open(expand_path(args.config), "r") as f_in:
        args.config = json.load(f_in)

    logger.debug("Running HiggsDNA analysis with config:")

    analysis = AnalysisManager(**vars(args))
    if args.short:
        analysis.run(max_jobs = 1)
    else:
        analysis.run()


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
