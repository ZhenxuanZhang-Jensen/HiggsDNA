import awkward
import glob
import json
import os
import numpy
import argparse
import math

from yahist import Hist1D
from yahist.utils import plot_stack
import matplotlib.pyplot as plt
from matplotlib import ticker

import mplhep as hep
plt.style.use(hep.style.CMS)

from higgs_dna.utils.logger_utils import setup_logger

def parse_arguments():
    parser = argparse.ArgumentParser(
            description="Assess HiggsDNA parquet output files")

    parser.add_argument(
        "--input_dir",
        required=True,
        type=str,
        help="path to directory with HiggsDNA parquet output files")

    parser.add_argument(
        "--make_tables",
        required=False,
        action="store_true",
        help="make data/MC tables")

    parser.add_argument(
        "--assess_systematics",
        required=False,
        action="store_true",
        help="make summary plots of each systematic")

    parser.add_argument(
        "--plots",
        required=False,
        default=None,
        type=str,
        help="path to json file with plot options")

    parser.add_argument(
        "--signals",
        required=False,
        default=None,
        type=str,
        help="csv list of signal processes (to be excluded from data/MC comparison)")

    parser.add_argument(
        "--group_procs",
        required=False,
        default=None,
        type=str,
        help="csv list of processes to group, with groups separated by '|', e.g. 'SMHiggs:ggH_M125,VBFH_M125|tt+X:ttGG,ttG,tt'")

    parser.add_argument(
        "--cuts",
        required=False,
        default=None,
        type=str,
        help="|-separated list of cuts, e.g. 'n_jets:[4,999]|LeadPhoton_mvaID:[-0.7,1.0]'")

    parser.add_argument(
        "--output_dir",
        required=False,
        default="output",
        type=str,
        help="output dir to store tables and plots in")

    parser.add_argument(
        "--make_webpage",
        required=False,
        action="store_true",
        help="make a nicely browsable web page")

    return parser.parse_args()


def regroup_processes(process_map, groups):
    for proc, id in process_map.items():
        process_map[proc] = [id]

    if groups is not None:
        groups = args.group_procs.split("|")
        for group in groups:
            new_name = group.split(":")[0]
            procs = group.split(":")[1].split(",")
            process_map[new_name] = []
            for proc in procs:
                process_map[new_name] += [process_map[proc]]
                process_map.pop(proc)

    return process_map

def infer_systematics(inputs, cuts):
    # Weight systs
    systs = {
        "weights" : {},
        "ics" : {}
    }
    for x in inputs:
        if "merged_nominal" in x:
            systs["ics"]["nominal"] = { "input" : x , "events" : awkward.from_parquet(x, lazy=True)}

        else:
            syst = x.split("/")[-1].replace("merged_", "").replace(".parquet", "").replace("_down", "").replace("_up", "")
            var = "up" if "up" in x else "down"
            if syst not in systs["ics"].keys():
                systs["ics"][syst] = {}
            systs["ics"][syst][var] = { "input" : x , "events" : awkward.from_parquet(x, lazy=True)}

    if cuts is not None:
        for syst, info in systs["ics"].items():
            if syst == "nominal":
                x = info["events"]
                for field, cut_range in cuts:
                    x = do_cut(x, field, cut_range)
                systs["ics"][syst]["events"] = x
                continue
            for var, var_info in info.items():
                x = var_info["events"]
                for field, cut_range in cuts:
                    x = do_cut(x, field, cut_range)

                systs["ics"][syst][var]["events"] = x
    
    events = systs["ics"]["nominal"]["events"] 

    weights = [x for x in events.fields if "weight_" in x]
    for x in weights:
        if "weight_central" in x:
            continue
        if "central" in x:
            var = "central"
            syst = x.replace("weight_", "").replace("_central", "")
        elif "up" in x or "down" in x:
            if "up" in x:
                var = x[x.find("up"):]
            elif "down" in x:
                var = x[x.find("down"):]
            syst = x.replace("weight_", "").replace("_"+var, "")

        if syst not in systs["weights"].keys():
            systs["weights"][syst] = {}

        systs["weights"][syst][var] = { "branch" : x }

    return systs


def events_with_ids(events, ids):
    cut = events["process_id"] == -999
    for id in ids:
        cut = (cut) | (events["process_id"] == id)
    return events[cut]


def do_cut(events, field, cut_range):
    cut = (events[field] >= cut_range[0]) & (events[field] <= cut_range[1])
    return events[cut]


def parse_cuts(cuts):
    if cuts is None:
        return None

    cuts_parsed = []

    cuts = cuts.split("|")
    for cut in cuts:   
        field = cut.split(":")[0]
        cut_range = [float(x.replace("[","").replace("]","")) for x in cut.split(":")[1].split(",")]
        cuts_parsed.append((field, cut_range))

    return cuts_parsed 


def get_yield_and_unc(events):
    return awkward.sum(events["weight_central"]) 


def table_line(proc, results, n_bkg):
    if n_bkg > 0:
        bkg_frac = results["n"] / n_bkg
    else:
        bkg_frac = 0.

    if proc == "total_bkg":
        proc = "Total MC bkg"
    line = "\t \t %s & %.3f & $\\pm \\text{%.3f}$  & $~^{+\\text{%.3f}}_{-\\text{%.3f}}$  & %.4f \\\\ \n" % (proc.replace("_", "-"), results["n"], results["stat_unc"], results["syst_unc_up"], results["syst_unc_down"], bkg_frac)
    return line

def make_tables(events, process_map, signals):
    yields = {}
   
    for proc, ids in process_map.items():
        yields[proc] = {
                "n" : 0,
                "stat_unc" : 0,
                "syst_unc_up" : 0,
                "syst_unc_down" : 0
        }
        nominal = events["ics"]["nominal"]["events"] 
        nominal = events_with_ids(nominal, ids)

        yields[proc]["n"] = awkward.sum(nominal["weight_central"])
        yields[proc]["stat_unc"] = awkward.sum(nominal["weight_central"] * nominal["weight_central"]) ** (0.5)

        if proc == "Data":
            continue

        for syst, info in events["ics"].items():
            if syst == "nominal":
                continue
            for var, var_info in info.items(): 
                syst_events = var_info["events"]
                syst_events = events_with_ids(syst_events, ids)

                n_syst = awkward.sum(syst_events["weight_central"])
                if n_syst >= yields[proc]["n"]:
                    yields[proc]["syst_unc_up"] += (yields[proc]["n"] - n_syst) ** 2
                elif n_syst < yields[proc]["n"]:
                    yields[proc]["syst_unc_down"] += (yields[proc]["n"] - n_syst) ** 2

        for syst, info in events["weights"].items():
            for var, var_info in info.items():
                if var == "central":
                    continue
                if "central" in info.keys():
                    syst_central = nominal[info["central"]["branch"]]
                else:
                    syst_central = awkward.ones_like(nominal["weight_central"])
                n_syst = awkward.sum(nominal["weight_central"] * (nominal[var_info["branch"]] / syst_central))

                if n_syst >= yields[proc]["n"]:
                    yields[proc]["syst_unc_up"] += (yields[proc]["n"] - n_syst) ** 2
                elif n_syst < yields[proc]["n"]:
                    yields[proc]["syst_unc_down"] += (yields[proc]["n"] - n_syst) ** 2

        yields[proc]["syst_unc_up"] = yields[proc]["syst_unc_up"] ** (0.5)
        yields[proc]["syst_unc_down"] = yields[proc]["syst_unc_down"] ** (0.5)

    table = ""
    table += "\\begin{center}\n"
    table += "\t \\begin{tabular}{ l | r r r | r}\n"
    table += "\t \t \\multicolumn{5}{c}{\\texttt{HiggsDNAPhysTool} : Data/MC Yield Table} \\\\ \\hline \\hline\n"
    table += "\t \t Process & Yield & Stat. unc. & Syst. unc. & $\\mathcal F$ of bkg \\\\ \\hline\n"

    yields["total_bkg"] = {
        "n" : 0,
        "stat_unc" : 0,
        "syst_unc_up" : 0,
        "syst_unc_down" : 0
    }
    for proc, result in yields.items():
        if proc == "Data" or proc in signals or proc == "total_bkg":
            continue
        yields["total_bkg"]["n"] += result["n"]
        for unc in ["stat_unc", "syst_unc_up", "syst_unc_down"]:
            yields["total_bkg"][unc] += result[unc] ** 2

    for unc in ["stat_unc", "syst_unc_up", "syst_unc_down"]:
        yields["total_bkg"][unc] = yields["total_bkg"][unc] ** (0.5)

   
    sort_idx = awkward.argsort([yields[x]["n"] for x in signals], ascending = True)
    signals = numpy.array(signals)[sort_idx] 
    for proc in signals:
        table += table_line(proc, yields[proc], yields["total_bkg"]["n"])

    table += "\t \t \\hline \n"
    
    bkgs = [x for x in yields.keys() if (x not in signals and x not in ["Data", "total_bkg"])]
    sort_idx = awkward.argsort([yields[x]["n"] for x in bkgs], ascending = True)
    bkgs = numpy.array(bkgs)[sort_idx]
    for proc in bkgs:
        table += table_line(proc, yields[proc], yields["total_bkg"]["n"])

    table += "\t \t \\hline \n"
    table += table_line("total_bkg", yields["total_bkg"], yields["total_bkg"]["n"])

    table += "\t \t \\hline \n"
        
    if "Data" in yields.keys():
        table += table_line("Data", yields["Data"], yields["total_bkg"]["n"])

    table += "\t \t \\hline \\hline \n"

    table += "\t \\end{tabular}\n" 
    table += "\\end{center}\n"

    print(table)
    return table


def make_data_mc_plot(data, bkg, sig, savename, **kwargs):
    normalize = kwargs.get("normalize", False)
    x_label = kwargs.get("x_label", None)
    y_label = kwargs.get("y_label", "Events" if not normalize else "Fraction of events")
    rat_label = kwargs.get("rat_label", "Data/MC")
    title = kwargs.get("title", None)
    y_lim = kwargs.get("y_lim", None)
    x_lim = kwargs.get("x_lim", None)
    rat_lim = kwargs.get("rat_lim", [0.0, 2.0])
    overflow = kwargs.get("overflow", False)
    log_y = kwargs.get("log_y", False)

    bins = kwargs.get("bins")

    h_data = Hist1D(data, bins=bins, overflow=overflow, label="Data")
    h_data = h_data.to_poisson_errors() # asymmetric poisson errors

    h_bkg = []
    h_bkg_syst = {}
    for proc, plot_data in bkg.items():
        h = Hist1D(plot_data["array"], weights = plot_data["weights"], bins = bins, overflow = overflow, label=proc)
        h = h.to_poisson_errors()
        h_bkg.append(h)
        for syst, syst_array, syst_weights in plot_data["syst_arrays"]:
            h_syst = Hist1D(syst_array, weights = syst_weights, bins = bins, overflow = overflow)
            h_syst = h_syst.to_poisson_errors()
            if syst not in h_bkg_syst.keys():
                h_bkg_syst[syst] = h_syst.copy().to_poisson_errors()
            else:
                h_bkg_syst[syst] += h_syst
        
        for syst, syst_weights in plot_data["syst_weights"]:
            h_syst = Hist1D(plot_data["array"], weights = syst_weights, bins = bins, overflow = overflow, label=proc)
            h_syst = h_syst.to_poisson_errors()
            if syst not in h_bkg_syst.keys():
                h_bkg_syst[syst] = h_syst.copy().to_poisson_errors()
            else:
                h_bkg_syst[syst] += h_syst


    h_bkg_total = None
    h_bkg_total_syst_up = None
    h_bkg_total_syst_down = None

    for h in h_bkg:
        if h_bkg_total is None:
            h_bkg_total = h.copy()
        else:
            h_bkg_total += h

    h_bkg_total = h_bkg_total.to_poisson_errors()

    h_bkg_total_syst_up = None
    h_bkg_total_syst_down = None


    h_sig = []
    for proc, plot_data in sig.items():
        h = Hist1D(plot_data["array"], weights = plot_data["weights"], bins = bins, overflow = overflow, label=proc)
        h = h.to_poisson_errors()
        h_sig.append(h)

    fig, (ax1,ax2) = plt.subplots(2, sharex=True, figsize=(12,9), gridspec_kw=dict(height_ratios=[3, 1]))
    plt.grid()
    h_data.plot(ax=ax1, color = "black", errors = True)
    plt.sca(ax1)
    hep.cms.label(" Preliminary",loc=0,data=True,lumi=137,fontsize=18)

    stack = sorted(h_bkg, key = lambda x : x.integral)
    plot_stack(stack, ax=ax1, histtype="stepfilled")

    for idx, h in enumerate(h_sig):
        h.plot(ax=ax1, color = "C%d" % (idx+5), errors = False, linewidth=3)

    ### Calculate errors ###
    # When plotting data errorbands in the ratio panel, show only the data error scaled by the total bkg (MC error will be plotted separately)
    ratio = h_data.divide(h_bkg_total, binomial = True)
    ratio._errors_up = h_data._errors_up / h_bkg_total.counts
    ratio._errors_down = h_data._errors_down / h_bkg_total.counts


    mc_stat_err_up = h_bkg_total._errors_up 
    mc_stat_err_down = h_bkg_total._errors_down 

    # Calculate syst error
    mc_syst_err_up = numpy.zeros(h_data.nbins)
    mc_syst_err_down = numpy.zeros(h_data.nbins)

    # Sum all sources in quadrature
    for syst, hist in h_bkg_syst.items():
        for i in range(hist.nbins):
            if hist.counts[i] >= h_bkg_total.counts[i]:
                mc_syst_err_up[i] += (hist.counts[i] - h_bkg_total.counts[i]) ** 2
            elif hist.counts[i] < h_bkg_total.counts[i]:
                mc_syst_err_down[i] += (hist.counts[i] - h_bkg_total.counts[i]) ** 2

    mc_syst_err_up = mc_syst_err_up ** (0.5)
    mc_syst_err_down = mc_syst_err_down ** (0.5)

    # Calculate total uncertainty by adding stat and syst in quadrature
    mc_stat_syst_err_up = ((mc_stat_err_up ** 2) + (mc_syst_err_up **2))**(0.5) 
    mc_stat_syst_err_down = ((mc_stat_err_down ** 2) + (mc_syst_err_down **2))**(0.5) 

    # Ratio pad errors
    mc_ratio_stat_err_up = 1 + (mc_stat_err_up / h_bkg_total.counts)
    mc_ratio_stat_err_down = 1 - (mc_stat_err_down / h_bkg_total.counts)

    mc_ratio_stat_syst_err_up = 1 + (mc_stat_syst_err_up / h_bkg_total.counts)
    mc_ratio_stat_syst_err_down = 1 - (mc_stat_syst_err_down / h_bkg_total.counts)

    ratio.metadata["label"] = None
    ratio.plot(ax=ax2, errors=True, color="black")

    if x_label is not None:
        ax2.set_xlabel(x_label)

    if y_label is not None:
        ax1.set_ylabel(y_label)

    if rat_label is not None:
        ax2.set_ylabel(rat_label)

    if title is not None:
        ax1.set_title(title)

    if y_lim is not None:
        ax1.set_ylim(y_lim)

    if rat_lim is not None:
        ax2.set_ylim(rat_lim)

    if x_lim is not None:
        ax1.set_xlim(x_lim)

    if log_y:
        ax1.set_yscale("log")

    for i in range(h_data.nbins):
        if h_bkg_total.counts[i] == 0:
            continue
        # Stat err 
        ax1.axhspan(h_bkg_total.counts[i] - mc_stat_err_down[i], h_bkg_total.counts[i] + mc_stat_err_up[i], float(i) / float(h_data.nbins), float(i+1) / float(h_data.nbins), color = "black", alpha = 0.25, linewidth = 0.)
        # Stat err ratio pad
        ax2.axhspan(mc_ratio_stat_err_down[i], mc_ratio_stat_err_up[i], float(i) / float(h_data.nbins), float(i+1) / float(h_data.nbins), color = "black", alpha = 0.25, linewidth = 0.)
        if h_bkg_syst:
            # Stat \oplus syst error
            ax1.axhspan(h_bkg_total.counts[i] + mc_stat_err_up[i], h_bkg_total.counts[i] + mc_stat_syst_err_up[i], float(i) / float(h_data.nbins), float(i+1) / float(h_data.nbins), color = "red", alpha = 0.25, linewidth = 0.) 
            ax1.axhspan(h_bkg_total.counts[i] - mc_stat_syst_err_down[i], h_bkg_total.counts[i] - mc_stat_err_down[i], float(i) / float(h_data.nbins), float(i+1) / float(h_data.nbins), color = "red", alpha = 0.25, linewidth = 0.) 
            # Stat \oplus syst error ratio pad
            ax2.axhspan(mc_ratio_stat_err_up[i], mc_ratio_stat_syst_err_up[i], float(i) / float(h_data.nbins), float(i+1) / float(h_data.nbins), color = "red", alpha = 0.25, linewidth = 0.)
            ax2.axhspan(mc_ratio_stat_syst_err_down[i], mc_ratio_stat_err_down[i], float(i) / float(h_data.nbins), float(i+1) / float(h_data.nbins), color = "red", alpha = 0.25, linewidth = 0.)

    plt.savefig(savename)
    plt.clf()


def make_shape_comparisons(plot_config, output_dir, events, process_map, signals):
    for field, info in plot_config.items():
        arrays = []
        weights = []
        names = []

        for proc, ids in process_map.items():
            if proc not in signals:
                continue
            evts = events_with_ids(events["ics"]["nominal"]["events"], ids)
            arrays.append(evts[field])
            weights.append(evts["weight_central"])
            names.append(proc)

        plot_shapes(arrays, weights, names, output_dir + "/%s_shape.pdf" % field, **info)


def make_plots(plot_config, output_dir, events, process_map, signals, bkgs):
    if "Data" not in process_map.keys():
        return

    for field, info in plot_config.items():
        data = events_with_ids(events["ics"]["nominal"]["events"], process_map["Data"])[field]
        bkg = {}
        sig = {}

        for b in bkgs:
            evts = events_with_ids(events["ics"]["nominal"]["events"], process_map[b])
            array = evts[field]
            weights = events_with_ids(events["ics"]["nominal"]["events"], process_map[b])["weight_central"]
            bkg[b] = { "array" : array, "weights" : weights, "syst_weights" : [], "syst_arrays" : [] }
            
            for syst, syst_info in events["ics"].items():
                if syst == "nominal":
                    continue
                for var, var_info in syst_info.items():
                    syst_array = events_with_ids(var_info["events"], process_map[b])[field]
                    syst_weights = events_with_ids(var_info["events"], process_map[b])["weight_central"]
                    bkg[b]["syst_arrays"].append((syst+"_"+var, syst_array, syst_weights))

            for syst, syst_info in events["weights"].items():
                for var, var_info in syst_info.items():
                    if var == "central":
                        continue
                    if "central" in syst_info.keys():
                        syst_central = evts[syst_info["central"]["branch"]]
                    else:
                        syst_central = awkward.ones_like(evts["weight_central"])

                    syst_weights = weights * (evts[var_info["branch"]] / syst_central)
                    bkg[b]["syst_weights"].append((syst+"_"+var, syst_weights)) 


        for s in signals:
            evts = events_with_ids(events["ics"]["nominal"]["events"], process_map[s])
            array = evts[field]
            weights = events_with_ids(events["ics"]["nominal"]["events"], process_map[s])["weight_central"]
            sig[s] = { "array" : array, "weights" : weights, "syst_weights" : [], "syst_arrays" : [] }

        make_data_mc_plot(data, bkg, sig, savename = "%s/%s_dataMC.pdf" % (output_dir, field), **info)


def get_quantile_range(array, ci):
    sorted = numpy.flip(numpy.sort(array))


def summarize_weights(systs, process_map, output_dir):
    #events = awkward.from_parquet(systs["ics"]["nominal"]["input"])
    events = systs["ics"]["nominal"]["events"]

    for proc, ids in process_map.items():
        evts_proc = events_with_ids(events, ids)
        for weight, vars in systs["weights"].items():
            for var, info in vars.items():
                branch = info["branch"]
                mean = awkward.mean(evts_proc[branch])
                std = awkward.std(evts_proc[branch])

                systs["weights"][weight][var][proc] = {
                    "mean" : mean,
                    "median" : numpy.quantile(evts_proc[branch], 0.5),
                    "p1sigma" : numpy.quantile(evts_proc[branch], 0.84),
                    "p2sigma" : numpy.quantile(evts_proc[branch], 0.975),
                    "m1sigma" : numpy.quantile(evts_proc[branch], 0.16),
                    "m2sigma" : numpy.quantile(evts_proc[branch], 0.025),
                }

    with open(output_dir + "/weight_syst_summary.json", "w") as f_out:
        json.dump(systs["weights"], f_out, indent = 4, sort_keys = True)


def summarize_ics(systs, process_map, output_dir):
    nominal = systs["ics"]["nominal"]["events"]

    for proc, ids in process_map.items():
        if proc == "Data":
            continue
        for syst, syst_info in systs["ics"].items():
            if syst == "nominal":
                continue

            nominal_proc = events_with_ids(nominal, ids)
            arrays = [nominal_proc["Diphoton_mass"]]
            weights = [nominal_proc["weight_central"]]
            names = ["Nominal"]
            
            for var, var_info in syst_info.items():
                syst_events = var_info["events"]
                syst_events = events_with_ids(syst_events, ids)

                arrays.append(syst_events["Diphoton_mass"])
                weights.append(syst_events["weight_central"])
                names.append(syst + " " + var)

            plot_shapes(arrays, weights, names, output_dir + "/%s_mgg_%s.pdf" % (syst, proc), x_label = "m_{gg} [GeV]", bins = "20, 120, 130")


def plot_shapes(arrays, weights, names, savename, **kwargs):
    x_label = kwargs.get("x_label", None)
    y_label = kwargs.get("y_label", "Fraction of events")
    rat_label = kwargs.get("rat_label", "Ratio to %s" % names[0])
    title = kwargs.get("title", None)
    y_lim = kwargs.get("y_lim", None)
    x_lim = kwargs.get("x_lim", None)
    rat_lim = kwargs.get("rat_lim", [0.0, 2.0])
    overflow = kwargs.get("overflow", False)
    log_y = kwargs.get("log_y", False)
    bins = kwargs.get("bins")

    hists = []
    norm = []
    mean = []
    std = []
    for a, w, n in zip(arrays, weights, names):
        h = Hist1D(a, weights = w, bins = bins, overflow = overflow) 
        norm.append(h.integral)
        mean.append(h.mean())
        std.append(h.std())
        h = h.normalize()
        hists.append(h)

    fig, (ax1,ax2) = plt.subplots(2, sharex=True, figsize=(12,9), gridspec_kw=dict(height_ratios=[3, 1]))
    plt.grid()

    for idx, h in enumerate(hists):
        h.plot(ax=ax1, color="C%d" % idx, errors = False, linewidth=3, label = "%s [N : %.3f, Mean : %.3f, Std : %.3f]" % (names[idx], norm[idx], mean[idx], std[idx]))

    plt.sca(ax1)
    hep.cms.label(" Preliminary",loc=0,data=True,lumi=137,fontsize=18)        
    
    if len(hists) >= 2:
        for i in range(1, len(hists)):
            ratio = hists[i].divide(hists[0])
            ratio.plot(ax=ax2, color="C%d" % i, errors = False, linewidth=3)

    if x_label is not None:
        ax2.set_xlabel(x_label)

    if y_label is not None:
        ax1.set_ylabel(y_label)

    if rat_label is not None:
        ax2.set_ylabel(rat_label)

    if title is not None:
        ax1.set_title(title)

    if y_lim is not None:
        ax1.set_ylim(y_lim)

    else:
        if log_y:
            ax1.set_yscale("log")
        y_min, y_max = ax1.get_ylim()
        ax1.set_ylim(y_min, y_max * (1 + (0.2 * len(hists))))

    if rat_lim is not None:
        ax2.set_ylim(rat_lim)

    if x_lim is not None:
        ax1.set_xlim(x_lim)

    if log_y:
        ax1.set_yscale("log")
 

    plt.savefig(savename)
    plt.clf()


def plot_weights(systs, process_map, output_dir):
    for proc, ids in process_map.items():
        if proc == "Data":
            continue
        means = []
        p1sigma = []
        p2sigma = []
        m1sigma = []
        m2sigma = []
        y = []
        labels = []
        
        for weight, vars in systs["weights"].items():
            for var, info in vars.items():
                means.append(info[proc]["mean"])
                labels.append(weight + "_" + var + " [mu = %.3f]" % info[proc]["mean"])
                p1sigma.append(info[proc]["p1sigma"])
                p2sigma.append(info[proc]["p2sigma"])
                m1sigma.append(info[proc]["m1sigma"])
                m2sigma.append(info[proc]["m2sigma"])


        sort_idx = awkward.argsort(means, ascending = False)

        means = numpy.array(means)[sort_idx]
        labels = numpy.array(labels)[sort_idx]
        p1sigma = numpy.array(p1sigma)[sort_idx]
        p2sigma = numpy.array(p2sigma)[sort_idx]
        m1sigma = numpy.array(m1sigma)[sort_idx]
        m2sigma = numpy.array(m2sigma)[sort_idx]

        fig, ax = plt.subplots()
        ax.set_xscale("log", base=10)

        plt.plot(means, numpy.arange(len(labels)), marker="o", color = "black", linewidth=0)
        x_min = 1. / 3.
        x_max = 3.0
        assert (1. / x_min) == (x_max / 1.) # need to have symmetric axes for the plotting to work
        x_range = x_max - x_min
        for i in range(len(means)):
            ax.axvspan(m2sigma[i], p2sigma[i], (i + 0.25) / float(len(means)), (i + 0.75) / float(len(means)), color = "yellow")
            ax.axvspan(m1sigma[i], p1sigma[i], (i + 0.25) / float(len(means)), (i + 0.75) / float(len(means)), color = "green")

            #ax.axhspan(i-0.25,i+0.25,(m2sigma[i] - x_min)/x_range,(p2sigma[i]-x_min)/x_range,color = "yellow")
            #ax.axhspan(i-0.25,i+0.25,(m1sigma[i] - x_min)/x_range,(p1sigma[i]-x_min)/x_range, color = "green")

        plt.xlim([x_min, x_max])
        plt.ylim([-0.5, len(means) - 0.5])
        plt.gca().xaxis.grid(True)
        ax.set_xticks(numpy.logspace(math.log(x_min, 10), math.log(x_max, 10), 7))
        ax.get_xaxis().set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        ax.tick_params(axis='x',which='minor',bottom=False,top=False,labelbottom=False)
        plt.yticks(numpy.arange(len(labels)), labels, fontsize = 10)

        plt.tight_layout()
        plt.savefig(output_dir + "/weight_syst_%s.pdf" % proc)
        plt.clf()


def main(args):
    logger = setup_logger("DEBUG")
    
    inputs = glob.glob(args.input_dir + "/merged_*.parquet")
    logger.debug("[HiggsDNABonusTool] Found %d input files in directory '%s'" % (len(inputs), args.input_dir))

    with open(args.input_dir + "/summary.json", "r") as f_in:
        process_map = json.load(f_in)["sample_id_map"]

    os.system("mkdir -p %s" % args.output_dir)

    signals = []
    if args.signals is not None:
        signals = args.signals.split(",")

    process_map = regroup_processes(process_map, args.group_procs)
    bkgs = [x for x in process_map.keys() if (x not in signals and x != "Data")]

    logger.debug("[HiggsDNABonusTool] Grouping processes by the following proc ids: ")
    for proc, ids in process_map.items():
        if proc == "Data":
            cat = proc
        elif proc in signals:
            cat = "Signal"
        else:
            cat = "Background"
        logger.debug("\t %s : %s (%s)" % (proc, str(ids), cat)) 
    logger.debug("[HiggsDNABonusTool] Processes marked as 'Signal' will not be included in data/MC comparisons and will be plotted as individual lines.")
    logger.debug("[HiggsDNABonusTool] Processes marked as 'Background' will be summed up to calculate the total background estimate and will be plotted in a stacked histogram in data/MC plots.")

    cuts = parse_cuts(args.cuts)
    events = infer_systematics(inputs, cuts)

    if args.make_tables:        
        logger.debug("[HiggsDNABonusTool] Making data/MC yield tables.")
        table = make_tables(events, process_map, signals)

        with open(args.output_dir + "/data_mc_yield_table.txt", "w") as f_out:
            f_out.write(table)

    if args.plots is not None:
        logger.debug("[HiggsDNABonusTool] Making data/MC plots.")
        with open(args.plots, "r") as f_in:
            plot_config = json.load(f_in)
            
        make_plots(plot_config, args.output_dir, events, process_map, signals, bkgs)
        logger.debug("[HiggsDNABonusTool] Making shape comparisons (no systematics).")
        make_shape_comparisons(plot_config, args.output_dir, events, process_map, signals)

    if args.assess_systematics:
        logger.debug("[HiggsDNABonusTool] Summarizing weight systematics.")
        summarize_weights(events, process_map, args.output_dir)
        plot_weights(events, process_map, args.output_dir)
        logger.debug("[HiggsDNABonusTool] Summarizing systematics with independent collections.")
        summarize_ics(events, process_map, args.output_dir)

    logger.debug("[HiggsDNABonusTool] Plots and tables written to directory '%s'" % (args.output_dir))

    if args.make_webpage:
        os.system("cp web/index.php %s" % args.output_dir)
        os.system("chmod 755 %s" % args.output_dir)
        os.system("chmod 755 %s/*" % args.output_dir)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)

