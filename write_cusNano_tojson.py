import os
import json

# Specify the input folder path
input_folder = '/eos/cms/store/group/phys_b2g/zhenxuan/custom_nanoAOD/YH_afterSKIM/UL2017_ttgg/GluGluToRadionToHHTo2G2Tau/'

# Create a dictionary to store the data
data = {}

# Traverse through all subfolders in the folder
for subdir in os.listdir(input_folder):
    # Create a dictionary to store the data for each subfolder
    subdir_data = {}
    
    # Fill in the data for each subfolder according to the given format
    subdir_data['process_id'] = 0 # process ID
    subdir_data['xs'] = 0.001 # cross section : -> signal:0.001, background: check in das
    subdir_data['bf'] = 0.000284658
    subdir_data['fpo'] = 5
    
    # Traverse through all files in the subfolder
    files = []
    if os.path.isdir(os.path.join(input_folder, subdir)):
        for file in os.listdir(os.path.join(input_folder, subdir)):
            file_path = os.path.join(input_folder, subdir, file)
            if os.path.isfile(file_path):
                files.append(os.path.join(input_folder, subdir) + "/" + file.strip())  # Remove leading and trailing whitespaces from the file name
    else:
        print(f"Warning: {os.path.join(input_folder, subdir)} is not a directory and will be skipped.")
    
    # Add the list of files to the data for the subfolder
    subdir_data['files'] = {'2017': files}
    
    # Add the subfolder data to the overall data
    data[subdir] = subdir_data

# Create a JSON file and write the data to it
with open('output.json', 'w') as f:
    json.dump(data, f, indent=4)