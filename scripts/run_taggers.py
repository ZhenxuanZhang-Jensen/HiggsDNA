import argparse
from higgs_dna.utils import setup_logger


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run taggers with specified configuration")

    # Required arguments


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

    return parser.parse_args()


def main(args):
    logger = setup_logger(args.log_level, args.log_file)
    logger.debug("Running script to apply taggers")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
