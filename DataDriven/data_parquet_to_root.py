import subprocess
import ROOT
import awkward as ak
import numpy as np
import uproot
import matplotlib.pyplot as plt
input_path = "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_FHSL/"
# input_path = "/eos/user/s/shsong/combined_WWgg/parquet/data/"
output_path = "/eos/user/z/zhenxuan/wwyy/root_files/"
# read the parquet file
awkward_array = ak.from_parquet(input_path + 'merged_nominal.parquet')
# modified the input parquet file to datadriven
# 定义要删除的字段名称
field_to_remove = "weight_central"
# 使用 列表推导式 删除指定字段
modified_array = awkward_array[[x for x in ak.fields(awkward_array) if x != field_to_remove]]
# 选择 cat
modified_array = modified_array[modified_array['category'] == 1]
# 保存成parquet，否则无法使用parquet_to_root
ak.to_parquet(modified_array, output_path+"tmp.parquet")

# parquet to root
output_root_name = 'Data_pre_datadriven.root'
from parquet_to_root import parquet_to_root
parquet_to_root(output_path+"tmp.parquet", output_path + output_root_name, treename='cat1', verbose=False)
