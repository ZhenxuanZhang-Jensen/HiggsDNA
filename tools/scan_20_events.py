import uproot
import awkward as ak
import numpy as np
import json

def compute_significance(signal_weight_sum, bkg_weight_sum):
    s = signal_weight_sum
    b = bkg_weight_sum
    return s / np.sqrt(s+b)

def get_weighted_events_in_cut_range(events, weights, lower_cut, upper_cut):
    return weights[(events > lower_cut) & (events <= upper_cut)]

def main():
    # Load ROOT files
    data_sideband = uproot.open("/eos/user/z/zhenxuan/wwyy/root_files/boosted_cat/Data_cat1.root:cat1")
    data = data_sideband.arrays()

    # Apply mask
    mask = ((data["Diphoton_mass"] < 115) | (data["Diphoton_mass"] > 135))
    data_sideband = data[mask]
    boosted_mass = ["M250","M260","M270","M280","M300","M320","M350","M400","M450","M550","M600","M650","M700","M750","M800","M850","M900","M1000","M1100","M1200","M1300","M1400","M1500","M1600","M1700","M1800","M1900","M2000","M2200","M2400","M2600","M2800","M3000"]
    signal_mass_points = boosted_mass
    signal_scores_dict = {}
    signal_weights_dict = {}

    # Load signal samples
    for mp in signal_mass_points:
        signal_file = "/eos/user/z/zhenxuan/hhwwgg_parquet/FHSL_channel/UL17_R_gghh_FHSL_" + mp + "_2017_cat1.parquet"
        signal = ak.from_parquet(signal_file)
        signal_scores_dict[mp] = ak.Array(signal["fatjet_H_Hqqqq_qqlv_vsQCDTop"])
        signal_weights_dict[mp] = ak.Array(signal["weight_central"])

    # Load background samples
    mc_bkg1 = uproot.open("/eos/user/z/zhenxuan/wwyy/root_files/boosted_cat/DataDriven_data_cat1.root:cat1")
    mc_bkg2 = uproot.open("/eos/user/z/zhenxuan/wwyy/root_files/boosted_cat/UL17_DiphotonJetsBox.root:cat1")
    mc_bkg1_scores = ak.Array(mc_bkg1["fatjet_H_Hqqqq_qqlv_vsQCDTop"].array())
    mc_bkg2_scores = ak.Array(mc_bkg2["fatjet_H_Hqqqq_qqlv_vsQCDTop"].array())
    mc_bkg1_weights = ak.Array(mc_bkg1["weight_central"].array())
    mc_bkg2_weights = ak.Array(mc_bkg2["weight_central"].array())

    # Sort PN scores
    data_sideband_scores = np.sort(data_sideband["fatjet_H_Hqqqq_qqlv_vsQCDTop"])[::-1]

    # Initialize variables
    idx = 0
    lower_cut = data_sideband_scores[idx]
    upper_cut = 1
    prev_total_significance = 0
    improvement = 1
    categories_dict = {mp: [] for mp in signal_mass_points}

    for mp in signal_mass_points:
        idx = 0
        lower_cut = data_sideband_scores[idx]
        upper_cut = 1
        prev_total_significance = 0
        improvement = 1

        while lower_cut >= 0.2 and idx < len(data_sideband_scores):
            idx += 10
            if idx < len(data_sideband_scores):
                lower_cut = data_sideband_scores[idx]
            else:
                lower_cut = 0.2

            total_significance = 0
            for cat_lower_cut, cat_upper_cut in categories_dict[mp] + [(lower_cut, upper_cut)]:
                signal_weights_in_cut = get_weighted_events_in_cut_range(signal_scores_dict[mp], signal_weights_dict[mp], cat_lower_cut, cat_upper_cut)
                mc_bkg1_weights_in_cut = get_weighted_events_in_cut_range(mc_bkg1_scores, mc_bkg1_weights, cat_lower_cut, cat_upper_cut)
                mc_bkg2_weights_in_cut = get_weighted_events_in_cut_range(mc_bkg2_scores, mc_bkg2_weights, cat_lower_cut, cat_upper_cut)

                # Compute significance for this category
                significance = compute_significance(np.sum(signal_weights_in_cut), np.sum(mc_bkg1_weights_in_cut) + np.sum(mc_bkg2_weights_in_cut))
                total_significance = np.sqrt(total_significance**2 + significance**2)

            improvement = (total_significance - prev_total_significance) / prev_total_significance if prev_total_significance != 0 else 1

            if improvement > 0.01:
                categories_dict[mp].append((lower_cut, upper_cut))
                upper_cut = lower_cut
                prev_total_significance = total_significance

        # Modify the last category to start from the lowest cut
        lowest_cut = 0.2
        categories_dict[mp][-1] = (lowest_cut, categories_dict[mp][-1][1])

        # Print the final list of categories and the total significance
        for mp in signal_mass_points:
            print(f"Categories for {mp}: {categories_dict[mp]}")

        # Write the categories to a JSON file
        with open("categories.json", "w") as f:
            json.dump(categories_dict, f)

        print(f"Total Significance: {prev_total_significance}")
if __name__ == "__main__":
    main()
