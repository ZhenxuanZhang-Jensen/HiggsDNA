import uproot
import json
with open("/afs/cern.ch/user/z/zhenxuan/HiggsDNA/metadata/samples/hhww_cusNANO_local.json") as fout:
    data = json.load(fout)
keys = ['UL17_dataB','UL17_dataC','UL17_dataD','UL17_dataE','UL17_dataF']
# keys = ["DoubleEG_Run2017B","DoubleEG_Run2017C","DoubleEG_Run2017D","DoubleEG_Run2017E","DoubleEG_Run2017F"]
# print(keys)
# read uproot files one by one
for key in range(len(keys)):
    files = data[keys[key]]['files']['2017']
    print("the key is: ", keys[key])
    print("files number is:",len(files) )
    for file in files:
        with uproot.open(file, timeout = 1800) as f:
            #attention new block to read lumi and save in the txt file to read whold data lumi to make sure we have enough lumi
            lumi = f['LuminosityBlocks'].arrays(['run','luminosityBlock'])
            Runs = f['Runs'].arrays(['run'])
            import itertools
            for i in range(len(Runs.run)):
                list_lumi = lumi.luminosityBlock[lumi.run == Runs.run[i]].to_list()
                range_list = [[t[0][1], t[-1][1]] for t in (tuple(g[1]) for g in itertools.groupby(enumerate(list_lumi), lambda list_lumi: list_lumi[1]-list_lumi[0]))]
                with open("/eos/user/z/zhenxuan/brilws/UL17_data_customized_v1.txt","b") as ftxt:
                    ftxt.write("\n")
                    # ftxt.write(str(file))
                    ftxt.write('"'+ str(Runs.run[i]) + '"')
                    ftxt.write(":")
                    ftxt.write(str(range_list))
                    ftxt.write(",")

# merge same keys

def merge_keys(input_json, output_json):
    def myhook(pairs):
        d = {}
        tmp_k = 0
        for k, v in pairs:
            if k not in d:
                d[k] = v
            else:
                d[k] += v
        return d
    with open(input_json,'r',encoding='utf8') as fp:
        json_data = json.load(fp,object_pairs_hook=myhook)


    with open(output_json,"w", encoding='utf-8') as f: ## 设置'utf-8'编码
        f.write("{")
        f.write('\n')
        for i in range(len(list(json_data.keys()))):
            f.write( '"' + str(list(json_data.keys())[i]) + '"' + ':'+ str(json_data[list(json_data.keys())[i]])+','  )  
            f.write('\n')
        f.write('\n')
        f.write("}")