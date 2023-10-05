import subprocess
import ROOT
import awkward as ak
import numpy as np
import uproot
import matplotlib.pyplot as plt
def parquet_to_root_forDataDriven(input_file, output_file, field_to_remove):
    # input_path = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/"
    
    # read the parquet file
    awkward_array = ak.from_parquet(input_file)
    # modified the input parquet file to datadriven
    # 定义要删除的字段名称
    field_to_remove = field_to_remove
    # 使用 列表推导式 删除指定字段
    modified_array = awkward_array[[x for x in ak.fields(awkward_array) if x != field_to_remove]]
    # 选择 cat
    modified_array = modified_array[modified_array['category'] == 1]
    # 保存成parquet，否则无法使用parquet_to_root
    ak.to_parquet(modified_array, "tmp.parquet")

    # parquet to root
    from parquet_to_root import parquet_to_root
    parquet_to_root("tmp.parquet", output_file, treename='cat1', verbose=False)
    return None
if __name__ == "__main__":
    input_file_diphoton_jets = "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1_wocut/UL17_DiPhotonJetsBox_MGG_80toInf_2017/merged_nominal.parquet"
    output_file_diphoton_jets = "/eos/user/z/zhenxuan/wwyy/root_files/UL17_DiphotonJetsBox_wocut.root"

    input_file_data_wocut = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL_cat1_wocut/merged_nominal.parquet"
    output_file_data_wocut = "/eos/user/z/zhenxuan/wwyy/root_files/Data_cat1_wocut.root"

    input_file_data = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL_cat1/merged_nominal.parquet"
    output_file_data = "/eos/user/z/zhenxuan/wwyy/root_files/Data_cat1.root"

    input_file_qcd = "/eos/user/z/zhenxuan/hhwwgg_parquet/qcd/hhwwgg_qcd_FHSL_cat1/merged_nominal.parquet"
    output_file_qcd = "/eos/user/z/zhenxuan/wwyy/root_files/Data_cat1.root"

    output_file_data_for_DataDriven = "/eos/user/z/zhenxuan/wwyy/root_files/Data_cat1_wocut_forDataDriven.root"

    input_file = input_file_data
    output_file = output_file_data

    field_to_remove = "nothing"
    parquet_to_root_forDataDriven(input_file= input_file, output_file= output_file, field_to_remove=field_to_remove)
    # datadriven need to remove the weight_central

    # field_to_remove = "weight_central"
    # parquet_to_root_forDataDriven(input_file= input_file_data, output_file= output_file_data_for_DataDriven, field_to_remove=field_to_remove)