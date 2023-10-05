import os
import argparse

def check_parquet_files(directory):
    parquet_dirs = []
    non_parquet_dirs = []
    
    for root, dirs, files in os.walk(directory):
        if root == directory:
            for subdir in dirs:
                subdir_path = os.path.join(root, subdir)
                if any(file.endswith('.parquet') for file in os.listdir(subdir_path)):
                    parquet_dirs.append(subdir_path)
                else:
                    non_parquet_dirs.append(subdir_path)
            break
    
    return parquet_dirs, non_parquet_dirs

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', type=str, default='/eos/user/z/zhenxuan',
                        required=False, help='')
    args = parser.parse_args()
    
    parquet_dirs, non_parquet_dirs = check_parquet_files(args.dir)
    
    # print("Subdirectories with Parquet files:")
    # for directory in parquet_dirs:
    #     print(directory)
    
    print("Subdirectories without Parquet files:")
    for directory in non_parquet_dirs:
        # only print the subdirectory name
        print(directory.split('/')[-1])
