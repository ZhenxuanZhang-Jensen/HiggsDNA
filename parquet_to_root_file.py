from parquet_to_root import parquet_to_root
import subprocess
import ROOT
import awkward as ak
import numpy as np
import uproot
import matplotlib.pyplot as plt
def parquet_to_root_forDataDriven(input_file, output_file, is_sig):
    # input_path = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/"
    
    # read the parquet file
    awkward_array = ak.from_parquet(input_file)
    # print fields
    print("input events fields", awkward_array.fields)
    # only keep the fields we want
    modified_array = awkward_array[['fatjetTop_btagDeepB','fatjetTop_eta','fatjetTop_lsf3','fatjetTop_mass','fatjetTop_msoftdrop','fatjetTop_particleNetMD_QCD','fatjetTop_particleNetMD_Xbb','fatjetTop_particleNet_TvsQCD','fatjetTop_particleNet_WvsQCD','fatjetTop_phi','fatjetTop_pt','fatjetTop_subJetIdx1','fatjetTop_subJetIdx2','fatjetTop_inclParTMDV1_HWWlvqqvsQCDTop','fatjetTop_inclParTMDV1_HWW4q3qvsQCD','fatjetTop_Hqqqq_vsQCDTop','fatjetTop_leadSubjet_btagDeepB','fatjetTop_subleadSubjet_btagDeepB','category','event','weight_central_initial','weight_central','weight_central_no_lumi','LP_weight','match_uncertainty','LP_weight_sys_up','LP_weight_sys_down']]
    # 保存成parquet，否则无法使用parquet_to_root
    ak.to_parquet(modified_array, "tmp.parquet")

    # parquet to root    
    parquet_to_root("tmp.parquet", output_file, treename='cat1', verbose=False)
    print("output events fields", modified_array.fields)
    return None


if __name__ == "__main__":    
    input_file = "/eos/user/z/zhenxuan/SWAN_projects/HH/calibration/eventsTTInclusive_matched.parquet"    
    output_file = "/eos/user/z/zhenxuan/SWAN_projects/HH/calibration/eventsTTInclusive_matched.root"

    parquet_to_root_forDataDriven(input_file= input_file, output_file= output_file,  is_sig=True)


    #####data
    # list_of_data = ["UL17_dataB_2017/merged_nominal","UL17_dataC_2017/merged_nominal","UL17_dataD_2017/merged_nominal","UL17_dataE_2017/merged_nominal","UL17_dataF_2017/merged_nominal"]
    # for i in range(len(list_of_data)):
    #     input_file_sig = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/" + list_of_data[i] + ".parquet"
    #     output_file_sig = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/" + list_of_data[i] + "_cat1.parquet"
    #     input_file = input_file_sig
    #     output_file = output_file_sig
    #     field_to_remove = "nothing"
        
    #     parquet_to_root_forDataDriven(input_file= input_file, output_file= output_file, field_to_remove=field_to_remove, is_sig=True)