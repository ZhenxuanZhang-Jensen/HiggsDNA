import awkward
import xgboost
import argparse
import glob
import numpy
import os
import json
import itertools
import vector

from sklearn import metrics

from higgs_dna.utils.logger_utils import setup_logger

def parse_arguments():
    parser = argparse.ArgumentParser(
            description="Train a BDT from HiggsDNA parquet output files")

    parser.add_argument(
        "--input_dir",
        required=True,
        type=str,
        help="path to directory with HiggsDNA parquet output files")

    parser.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="path to output directory")

    parser.add_argument(
        "--config",
        required=True,
        type=str,
        help="json file with BDT config options")

    parser.add_argument(
        "--cuts",
        required=False,
        default=None,
        type=str,
        help="|-separated list of cuts, e.g. 'n_jets:[4,999]|LeadPhoton_mvaID:[-0.7,1.0]'") 

    return parser.parse_args()


def do_cut(events, field, range):
    cut = (events[field] >= range[0]) & (events[field] <= range[1])
    return events[cut]


def parse_cuts(cuts):
    if cuts is None:
        return None

    cuts_parsed = []

    cuts = cuts.split("|")
    for cut in cuts:
        field = cut.split(":")[0]
        range = [float(x.replace("[","").replace("]","")) for x in cut.split(":")[1].split(",")]
        cuts_parsed.append((field, range))

    return cuts_parsed


def main(args):
    logger = setup_logger("DEBUG")
    

    nominal = args.input_dir + "/merged_nominal.parquet"
    logger.debug("[HiggsDNABonusTool] Loading events from file '%s' to train BDT." % nominal)

    events = awkward.from_parquet(nominal) 
    
    cuts = parse_cuts(args.cuts)
    if cuts is not None:
        for field, range in cuts:
            events = do_cut(events, field, range)    

    with open(args.config, "r") as f_in:
        bdt_config = json.load(f_in)

    logger.debug("[HiggsDNABonusTool] Training BDT with the following options:")
    for key, val in bdt_config.items():
        logger.debug("\t %s : %s" % (str(key), str(val)))

    with open(args.input_dir + "/summary.json", "r") as f_in:
        process_map = json.load(f_in)["sample_id_map"]

    # Anchor by diphoton_eta
    eta_fields = [x for x in bdt_config["features"] if "eta" in x]
    for field in eta_fields:
        events[field] = events[field] * numpy.sign(events.Diphoton_eta)

    # Calculate all delta etas
    #combinations = itertools.combinations(eta_fields, 2)
    #for x in combinations:
    #    name = "delta_eta_%s_%s" % (x[0], x[1])
    #    events[name] = abs(events[x[0]] - events[x[1]])
    #    bdt_config["features"].append(name)

    """
    jets = []
    for i in range(1,9):
        jet = vector.awk({
            "pt" : events["jet_%d_pt"],
            "eta" : events["jet_%d_eta"],
            "phi" : events["jet_%d_phi"],
            "mass" : events["jet_%d_mass"]
        })
        jet = awkward.unflatten(jet, counts=1)
        jets.append(jet)

    jets = awkward.concatentate(jets, axis =1)

    dijets = awkward.combinations(jets, 2, fields = ["l", "sl"])
    dijets = dijets[awkward.argsort(dijets.l.pt + dijets.sl.pt, ascending=False, axis=1)]
    dijets["comp"] = dijets.l + dijets.sl
    trijets = awkward.combinations(jets, 3, fields = ["a", "b", "c"])
    trijets = trijets[awkward.argsort(trijets.a.pt + trijets.b.pt + trijets.c.pt)]
    """


    # Add label
    events["train_label"] = awkward.ones_like(events.event) * -1

    logger.debug("[HiggsDNABonusTool] Found processes by the following proc ids: ")
    for proc, id in process_map.items():
        if proc == "Data":
            cat = proc
            label = -1
        elif proc in bdt_config["signal"]:
            cat = "Signal"
            label = 1
        elif proc in bdt_config["background"]:
            cat = "Background"
            label = 0
        else:
            continue
        logger.debug("\t %s : %s (%s)" % (proc, str(id), cat))

        events["train_label"] = awkward.where(
                events.process_id == id,
                awkward.ones_like(events.train_label) * label,
                events.train_label
        )

    logger.debug("[HiggsDNABonusTool] Processes marked as 'Signal' or 'Background' will be used to train the BDT as signal and background, respectively.")

    logger.debug("[HiggsDNABonusTool] Out of %d total events, found %d signal and %d background events." % (len(events), len(events[events.train_label == 1]), len(events[events.train_label == 0])))

    # Add test/train/val split     
    split = numpy.random.randint(3, size = len(events))
    events["train_split"] = awkward.from_numpy(split) # 0 = train, 1 = test, 2 = val

    events = awkward.nan_to_num(events, nan=-999., posinf=-999., neginf=-999.)

    # Select only signal/background events
    events_bdt = events[events.train_label >= 0]
    events_bdt_train = awkward.values_astype(events_bdt[events_bdt.train_split == 0], numpy.float64)
    events_bdt_test = awkward.values_astype(events_bdt[events_bdt.train_split == 1], numpy.float64)

    features_train = awkward.to_numpy(events_bdt_train[bdt_config["features"]])
    features_train = features_train.view((float, len(features_train.dtype.names)))

    features_test = awkward.to_numpy(events_bdt_test[bdt_config["features"]])
    features_test = features_test.view((float, len(features_test.dtype.names))) 

    # Make dmatrix for xgboost
    d_train = xgboost.DMatrix(
            features_train,
            awkward.to_numpy(events_bdt_train["train_label"]),
            weight = awkward.to_numpy(abs(events_bdt_train["weight_central"]))
    )

    d_test = xgboost.DMatrix(
            features_test,
            awkward.to_numpy(events_bdt_test["train_label"]),
            weight = awkward.to_numpy(abs(events_bdt_test["weight_central"]))
    )
 
    eval_list = [(d_train, "train"), (d_test, "test")]
    progress = {}
    bdt_config["mva"]["param"]["scale_pos_weight"] = awkward.sum(events_bdt[events_bdt.train_label == 0]["weight_central"]) / awkward.sum(events_bdt[events_bdt.train_label == 1]["weight_central"])
    bdt = xgboost.train(
            bdt_config["mva"]["param"],
            d_train,
            bdt_config["mva"]["n_trees"],
            eval_list, evals_result = progress,
            early_stopping_rounds = bdt_config["mva"]["early_stopping_rounds"]
    )

    os.system("mkdir -p %s" % (args.output_dir))
    bdt.save_model(args.output_dir + "/weights.xgb")

    # Predict
    evts = awkward.values_astype(events, numpy.float64)
    features = awkward.to_numpy(evts[bdt_config["features"]])
    features = features.view((float, len(features.dtype.names)))
    events["mva_score"] = bdt.predict(xgboost.DMatrix(features))

    awkward.to_parquet(events, args.output_dir + "/merged_nominal.parquet")

    # Copy config files to output dir
    os.system("cp %s %s" % (args.input_dir + "/summary.json", args.output_dir))
    os.system("cp %s %s" % (args.config, args.output_dir))


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
