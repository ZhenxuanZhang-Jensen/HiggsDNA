import json
jsonFile = "/eos/user/z/zhenxuan/brilws/lumi_cal.json"
com_list = json.load(jsonFile)
print(com_list.keys())