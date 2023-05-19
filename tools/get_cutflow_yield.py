import json
import pandas as pd
import os
###########################################################################
folder_path = "/eos/user/s/shsong/combined_WWgg/parquet/test"

# List to store the file paths
file_paths = []

# Recursive function to traverse the directory
def search_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if "combined_eff.json" in file:
                file_path = os.path.join(root, file)
                file_paths.append(file_path)

# Call the recursive function to search for files
search_files(folder_path)

# Print the file paths
file_list=[]
for file_path in file_paths:
    file_list.append(file_path)

###########################################################################
i=1
dfs = {}
for file in file_list:
    # List to store the dictionaries from each file
    data_list = []

    # File names
    with open(file, 'r') as f:
            # Load the JSON data from the file
        data = json.load(f)

            # Extend the data_list with the dictionaries from each file
        data_list.extend(data)

    # Initialize empty lists to store the values
    eff_field_names = []
    eff_data_dict = {}
    field_names = []
    data_dict = {}
    # Process the data
    for data in data_list:
        # Iterate over the keys in the dictionary
        for key, value in data.items():
            # Extract field name and data
            if "efficiency" in key:
                eff_field_name = key.split("-")[1]
                eff_field_data = value
            if "event number" in key:
                field_name = key.split("-")[1]
                field_data = value
            if eff_field_name not in eff_field_names:
                eff_field_names.append(eff_field_name)
            if field_name not in field_names:
                field_names.append(field_name)
            # Add the data to the corresponding field in the dictionary
            
            if eff_field_name not in eff_data_dict:
                eff_data_dict[eff_field_name] = [eff_field_data]
                
            if field_name not in data_dict:
                data_dict[field_name] = [field_data]
    dfs[f"dfeff{i}"]=pd.DataFrame(eff_data_dict)
    dfs[f"df{i}"]=pd.DataFrame(data_dict)
    i=i+1

###########################################################################
combined_effdf = pd.DataFrame(0, index=dfs[f"dfeff{1}"].index, columns=dfs[f"dfeff{1}"].columns)
for j in range(1, i):
    key = f"dfeff{j}"
    if key in dfs:
        combined_effdf = combined_effdf+dfs[key]
    df_eff=combined_effdf/(i-1)
eventdf = pd.DataFrame(0, index=dfs[f"df{1}"].index, columns=dfs[f"df{1}"].columns)
for j in range(1, i):
    key = f"df{j}"
    if key in dfs:
        eventdf = eventdf+dfs[key]
output="/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/test/"

eventdf.to_parquet(output+"event_yield.parquet")
df_eff.to_parquet(output+"cutflow_eff.parquet")