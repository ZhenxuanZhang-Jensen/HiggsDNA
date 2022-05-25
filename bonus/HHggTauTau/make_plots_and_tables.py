import os
import argparse

from parallel_utils import submit_jobs

def parse_arguments():
    parser = argparse.ArgumentParser(
            description="Assess HiggsDNA parquet output files")

    parser.add_argument(
        "--input_dir",
        required=True,
        type=str,
        help="path to directory with HiggsDNA parquet output files")

    

    return parser.parse_args()

BASE = """python assess.py --input_dir "INPUT_DIR" --output_dir OUTPUT_DIR --group_procs "GJets:GJets_HT-100To200,GJets_HT-200To400,GJets_HT-400To600,GJets_HT-40To100,GJets_HT-600ToInf|tt+X:TTGG,TTGamma,TTJets|VGamma:ZGamma,WGamma|SingleH:VBFH_M125,VH_M125,ggH_M125,ttH_M125|HH:HHggTauTau,HHggWW_dileptonic,HHggWW_semileptonic" --signals "HH,SingleH" --blind --backgrounds "tt+X,VGamma,DiPhoton,GJets" --plots "HHggTauTau/plots.json" --make_webpage --assess_systematics"""

BASE_TABLES = """python assess.py --input_dir "INPUT_DIR" --output_dir OUTPUT_DIR --make_tables --group_procs "GJets:GJets_HT-100To200,GJets_HT-200To400,GJets_HT-400To600,GJets_HT-40To100,GJets_HT-600ToInf|tt+X:TTGG,TTGamma,TTJets" --signals "HHggTauTau,HHggWW_dileptonic,HHggWW_semileptonic,VBFH_M125,VH_M125,ggH_M125,ttH_M125" --blind --backgrounds 'tt+X,WGamma,ZGamma,DiPhoton,GJets'"""

cuts = [
    ("inclusive", None),
    ("1tau_0lep", "category:[8,8]"),
    ("1tau_1IsoTrack", "category:[7,7]"),
    ("2tau_0lep", "category:[3,3]"),
    ("1tau_1lep", "category:[1,2]"),
    ("0tau_2lep", "category:[4,6]")
]

def main(args):
    os.chdir("../")

    base_command = BASE.replace("INPUT_DIR", args.input_dir)
    table_command = BASE_TABLES.replace("INPUT_DIR", args.input_dir)

    command_list = []
    for name, cut in cuts:
        command = base_command.replace("OUTPUT_DIR", "~/public_html/HiggsDNA/HHggTauTau_25May2022/" + name)
        tab_command = table_command.replace("OUTPUT_DIR", "~/public_html/HiggsDNA/HHggTauTau_25May2022/" + name)
        if cut is not None:
            command += ' --cuts "%s"' % cut
            tab_command += ' --cuts "%s"' % cut
            continue
        command_list.append(command)
        #command_list.append(tab_command)

    submit_jobs(command_list, 12)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
