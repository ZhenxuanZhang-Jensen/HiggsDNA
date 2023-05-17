import uproot
import awkward as ak
import numpy as np
import json

def get_events_in_cut_range(events, weights, lower_cut, upper_cut, output_file_name):
    selected_events = events[(events["fatjet_H_Hqqqq_qqlv_vsQCDTop"] > lower_cut) & (events["fatjet_H_Hqqqq_qqlv_vsQCDTop"] <= upper_cut)]
    return selected_events

def main():
    # Load ROOT files
    data_sideband = uproot.open("/eos/user/z/zhenxuan/wwyy/root_files/boosted_cat/Data_cat1.root:cat1")
    data = data_sideband.arrays()

    # Apply mask
    mask = ((data["Diphoton_mass"] < 115) | (data["Diphoton_mass"] > 135))
    data_sideband = data[mask]
    boosted_mass = ["M1000","M1100","M1200","M1300","M1400","M1500","M1600","M1700","M1800","M1900","M2000","M2200","M2400","M2600","M2800","M3000"]
    # boosted_mass_test = ["M250","M260"]
    signal_mass_points = boosted_mass
    signal_scores_dict = {}
    signal_weights_dict = {}

    # Load signal samples
    for mp in signal_mass_points:
        signal_file = "/eos/user/z/zhenxuan/hhwwgg_parquet/FHSL_channel/UL17_R_gghh_FHSL_" + mp + "_2017_cat1.parquet"
        signal = ak.from_parquet(signal_file)

    # Sort PN scores
    data_sideband_scores = np.sort(data_sideband["fatjet_H_Hqqqq_qqlv_vsQCDTop"])[::-1]

    # Initialize variables
    idx = 0
    lower_cut = data_sideband_scores[idx]
    upper_cut = 1
    prev_limit = 0
    improvement = 1
    categories_dict = {mp: [] for mp in signal_mass_points}

    for mp in signal_mass_points:
        idx = 0
        lower_cut = data_sideband_scores[idx]
        upper_cut = 1
        prev_limit = 0
        improvement = 1

        while lower_cut >= 0.2 and idx < len(data_sideband_scores):
            idx += 10
            if idx < len(data_sideband_scores):
                lower_cut = data_sideband_scores[idx]
            else:
                lower_cut = 0.2

            total_limit = 0
            for cat_lower_cut, cat_upper_cut in categories_dict[mp] + [(lower_cut, upper_cut)]:
                signal_name_for_scan = "Signal_"+'tmp'+"_FHSL_2017_1jets_"+"_for_scan_limit"+".root"
                data_sideband_name_for_scan = "Data_FHSL_2017_cat_1jets_"+'_for_scan_limit'+ "_"+ 'tmp' + ".root"
                signal_selcted_events = get_events_in_cut_range(signal, lower_cut, upper_cut)
                data_sideband_selcted_events = get_events_in_cut_range(data_sideband, lower_cut, upper_cut)
            # Compute limit for this category
            limit = get_limit(data_weights_in_cut, signal_weights_in_cut)

            total_limit += limit

        improvement = (total_limit - prev_limit) / prev_limit if prev_limit != 0 else 1

        if improvement > 0.01:
            categories_dict[mp].append((lower_cut, upper_cut))
            upper_cut = lower_cut
            prev_limit = total_limit

    # Modify the last category to start from the lowest cut
    lowest_cut = 0.2
    categories_dict[mp][-1] = (lowest_cut, categories_dict[mp][-1][1])

    # Print the final list of categories and the total limit
    for mp in signal_mass_points:
        print(f"Categories for {mp}: {categories_dict[mp]}")

    print(f"Total Limit: {prev_limit}")

    # Write the categories to a JSON file
    with open("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/categories_all_mass_point.json", "w") as f:
        json.dump(categories_dict, f, indent=1)
