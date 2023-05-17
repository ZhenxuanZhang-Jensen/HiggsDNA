from parquet_to_root import parquet_to_root
import subprocess
import ROOT
import awkward as ak
import numpy as np
import uproot
import matplotlib.pyplot as plt
def parquet_to_root_forDataDriven(input_file, output_file, field_to_remove, is_sig):
    # input_path = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/"
    
    # read the parquet file
    awkward_array = ak.from_parquet(input_file)
    # print fields
    print("input events fields", awkward_array.fields)
    # modified the input parquet file to datadriven
    # 定义要删除的字段名称
    field_to_remove = field_to_remove
    # 使用 列表推导式 删除指定字段
    modified_array = awkward_array[[x for x in ak.fields(awkward_array) if x != field_to_remove]]
    # 加字段 minphotonID, maxphotonID
    if is_sig == True:
        modified_array = modified_array[modified_array['Diphoton_minID']>-0.9]
    # 选择 cat
    modified_array = modified_array[modified_array['category'] == 1]
    # 保存成parquet，否则无法使用parquet_to_root
    ak.to_parquet(modified_array, "tmp.parquet")

    # parquet to root    
    parquet_to_root("tmp.parquet", output_file, treename='cat1', verbose=False)
    print("output events fields", modified_array.fields)
    return None


if __name__ == "__main__":
    # list_of_all_FHSL_sig = ["UL17_R_gghh_FHSL_M-400_2017","UL17_R_gghh_FHSL_M-2600_2017","UL17_R_gghh_FHSL_M-3000_2017","UL17_R_gghh_FHSL_M-700_2017","UL17_R_gghh_FHSL_M-1900_2017","UL17_R_gghh_FHSL_M-800_2017","UL17_R_gghh_FHSL_M-550_2017","UL17_R_gghh_FHSL_M-280_2017","UL17_R_gghh_FHSL_M-260_2017","UL17_R_gghh_FHSL_M-2800_2017","UL17_R_gghh_FHSL_M-320_2017","UL17_R_gghh_FHSL_M-650_2017","UL17_R_gghh_FHSL_M-1200_2017","UL17_R_gghh_FHSL_M-2400_2017","UL17_R_gghh_FHSL_M-1600_2017","UL17_R_gghh_FHSL_M-1000_2017","UL17_R_gghh_FHSL_M-1700_2017","UL17_R_gghh_FHSL_M-1100_2017","UL17_R_gghh_FHSL_M-900_2017","UL17_R_gghh_FHSL_M-1300_2017","UL17_R_gghh_FHSL_M-850_2017","UL17_R_gghh_FHSL_M-1800_2017","UL17_R_gghh_FHSL_M-270_2017","UL17_R_gghh_FHSL_M-350_2017","UL17_R_gghh_FHSL_M-600_2017","UL17_R_gghh_FHSL_M-300_2017","UL17_R_gghh_FHSL_M-750_2017","UL17_R_gghh_FHSL_M-250_2017","UL17_R_gghh_FHSL_M-1400_2017","UL17_R_gghh_FHSL_M-1500_2017","UL17_R_gghh_FHSL_M-2200_2017","UL17_R_gghh_FHSL_M-450_2017","UL17_R_gghh_FHSL_M-2000_2017"]
    # list_of_all_FHSL_sig = ["UL17_R_gghh_FHSL_M-2800_2017","UL17_R_gghh_FHSL_M-2400_2017","UL17_R_gghh_FHSL_M-2000_2017","UL17_R_gghh_FHSL_M-1400_2017","UL17_R_gghh_FHSL_M400_2017","UL17_R_gghh_FHSL_M450_2017","UL17_R_gghh_FHSL_M550_2017","UL17_R_gghh_FHSL_M650_2017","UL17_R_gghh_FHSL_M700_2017","UL17_R_gghh_FHSL_M750_2017","UL17_R_gghh_FHSL_M800_2017","UL17_R_gghh_FHSL_M850_2017","UL17_R_gghh_FHSL_M900_2017","UL17_R_gghh_FHSL_M600_2017","UL17_R_gghh_FHSL_M1600_2017","UL17_R_gghh_FHSL_M1200_2017","UL17_R_gghh_FHSL_M1900_2017","UL17_R_gghh_FHSL_M1500_2017","UL17_R_gghh_FHSL_M1100_2017","UL17_R_gghh_FHSL_M270_2017","UL17_R_gghh_FHSL_M1300_2017","UL17_R_gghh_FHSL_M280_2017","UL17_R_gghh_FHSL_M260_2017","UL17_R_gghh_FHSL_M350_2017","UL17_R_gghh_FHSL_M3000_2017","UL17_R_gghh_FHSL_M2600_2017","UL17_R_gghh_FHSL_M320_2017","UL17_R_gghh_FHSL_M300_2017","UL17_R_gghh_FHSL_M1700_2017","UL17_R_gghh_FHSL_M1000_2017","UL17_R_gghh_FHSL_M250_2017","UL17_R_gghh_FHSL_M2200_2017","UL17_R_gghh_FHSL_M1800_2017"]
    list_of_all_FHSL_sig = ["UL17_R_gghh_FHSL_M3000_2017"]
    for i in range(len(list_of_all_FHSL_sig)):
        input_file_sig = "/eos/user/z/zhenxuan/hhwwgg_parquet/FHSL_channel/" + list_of_all_FHSL_sig[i] + "_cat1.parquet"
        output_file_sig = "/eos/user/z/zhenxuan/wwyy/root_files/" + list_of_all_FHSL_sig[i] + "_cat1.root"
        input_file = input_file_sig
        output_file = output_file_sig
        field_to_remove = "nothing"        
        parquet_to_root_forDataDriven(input_file= input_file, output_file= output_file, field_to_remove=field_to_remove, is_sig=True)


    #####data
    # list_of_data = ["UL17_dataB_2017/merged_nominal","UL17_dataC_2017/merged_nominal","UL17_dataD_2017/merged_nominal","UL17_dataE_2017/merged_nominal","UL17_dataF_2017/merged_nominal"]
    # for i in range(len(list_of_data)):
    #     input_file_sig = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/" + list_of_data[i] + ".parquet"
    #     output_file_sig = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/" + list_of_data[i] + "_cat1.parquet"
    #     input_file = input_file_sig
    #     output_file = output_file_sig
    #     field_to_remove = "nothing"
        
    #     parquet_to_root_forDataDriven(input_file= input_file, output_file= output_file, field_to_remove=field_to_remove, is_sig=True)