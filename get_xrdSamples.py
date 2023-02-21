import os 
cmd = 'dasgoclient --query="file dataset=' + "/DoubleEG/Run2017F-UL2017_MiniAODv2-v2/MINIAOD"+ '"'
root_files = os.popen(str(cmd)).readlines()
print(root_files)