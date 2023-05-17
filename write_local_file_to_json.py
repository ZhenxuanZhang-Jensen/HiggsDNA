import os
import json

# 指定输入文件夹路径
input_folder = '/eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/_das_samples_top_jets_calibration'

# 创建一个字典来存储数据
data = {}

# 遍历文件夹中的所有子文件夹
for subdir in os.listdir(input_folder):
    # 创建一个字典来存储每个子文件夹的数据
    subdir_data = {}
    
    # 按照给定的格式填写每个子文件夹的数据
    subdir_data['xs'] = 4.958
    subdir_data['bf'] = 1.0
    subdir_data['fpo'] = 5
    
    # 遍历子文件夹中的所有文件
    files = []
    if os.path.isdir(os.path.join(input_folder, subdir)):
        for file in os.listdir(os.path.join(input_folder, subdir)):
            file_path = os.path.join(input_folder, subdir, file)
            if os.path.isfile(file_path):
                files.append(os.path.join(input_folder, subdir) + "/" + file.strip()) # 去除文件名前后的空格
    else:
        print(f"Warning: {os.path.join(input_folder, subdir)} is not a directory and will be skipped.")
    
    # 将文件列表添加到子文件夹的数据中
    subdir_data['files'] = {'2017': files}
    
    # 将子文件夹数据添加到总数据中
    data[subdir] = subdir_data

# 创建JSON文件并将数据写入其中
with open('output.json', 'w') as f:
    json.dump(data, f, indent=4)
