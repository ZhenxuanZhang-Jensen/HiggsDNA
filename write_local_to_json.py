import os
import json

# Specify the input folder path
# input_folder = '/eos/cms/store/group/phys_b2g/jodervan/NMSSM_XToYHTo2G2WTo2G4Q_UL2017/'
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2018_signal_YH_SL/UL2018'
UL17_HH_FH = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2017_signal/UL2017/'
UL17_HH_bbgg = "/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2017_bbgg/UL2017"
UL17_YH_bbgg = "/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2017_bbgg/UL2017"
# UL17_HH_SL = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2017_signal_SL/UL2017'
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2017/UL2017'#17data
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2018/UL2018'#17data
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2016APV/UL2016APV/'#16Predata
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2016_signal/UL2016/'#16Predata
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2017_bkg/UL2017/'#2017
# input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/nanoAOD_24Oct2023/UL2018_bkg/UL2018/'#2017
input_folder_list = []
input_folder_list.append(UL17_HH_bbgg)
# input_folder_list.append(UL17_YH_bbgg)
isdata = False
year = '2017'
is_slimed = True
# Create a dictionary to store the data
data = {}
subdir_data = {}
for input_folder in input_folder_list:
    # Traverse through all subfolders in the folder
    for subdir in os.listdir(input_folder):
        # Create a dictionary to store the data for each subfolder
        
        # Fill in the data for each subfolder according to the given format
        if isdata:
            subdir_data['process_id'] = 0
            subdir_data['fpo'] = 3
        else:    
            subdir_data['xs'] = 0.001 # cross section : -> signal:0.001, background: check in das
            subdir_data['bf'] = 1.0
            subdir_data['process_id'] = 0
            subdir_data['fpo'] = 3
            
        # Traverse through all files in the subfolder
        files = []
        if os.path.isdir(os.path.join(input_folder, subdir)):
            count = 0
            for file in os.listdir(os.path.join(input_folder, subdir)):
                if is_slimed:
                    count +=1
                if count >= 10:
                    break
                file_path = os.path.join(input_folder, subdir, file)
                if os.path.isfile(file_path):
                    files.append(os.path.join(input_folder, subdir) + "/" + file.strip())  # Remove leading and trailing whitespaces from the file name
        else:
            print(f"Warning: {os.path.join(input_folder, subdir)} is not a directory and will be skipped.")
        
        # Add the list of files to the data for the subfolder
        subdir_data['files'] = {year: files}
        
        # Add the subfolder data to the overall data
        data[subdir] = subdir_data

# Create a JSON file and write the data to it
with open('output.json', 'w') as f:
    json.dump(data, f, indent=4)