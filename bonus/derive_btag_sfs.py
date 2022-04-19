import awkward
import glob
import argparse
import json
import numpy
import os

from higgs_dna.utils.logger_utils import setup_logger

def parse_arguments():
    parser = argparse.ArgumentParser(
            description="Derive btag normalization sfs from HiggsDNA parquet output files")

    parser.add_argument(
        "--input_dir",
        required=True,
        type=str,
        help="path to directory with HiggsDNA parquet output files")

    parser.add_argument(
        "--output_dir",
        required=False,
        default="output",
        type=str,
        help="output dir to store json in")

    return parser.parse_args()


def get_proc_and_year(f):
    subdir = f.split("/")[-2]
    return subdir


def get_sfs(f):
    result = {}
    events = awkward.from_parquet(f)

    weights = [x for x in events.fields if "weight_btag_deepjet_sf_SelectedJet_" in x]
    for w in weights:
        variation = w.replace("weight_btag_deepjet_sf_SelectedJet_", "")
        result[variation] = 1. / awkward.mean(events[w])

    return result


def main(args):
    logger = setup_logger("DEBUG")

    inputs = glob.glob(args.input_dir + "/*/merged_nominal.parquet")
    logger.debug("[HiggsDNABonusTool] Found %d input files in directory '%s'" % (len(inputs), args.input_dir))

    os.system("mkdir -p %s" % args.output_dir)

    results = {}

    for f in inputs:
        name = get_proc_and_year(f)
        sfs = get_sfs(f)
        results[name] = sfs

    with open(args.output_dir + "/btag_sfs.json", "w") as f_out:
        json.dump(results, f_out, indent = 4, sort_keys = True)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
