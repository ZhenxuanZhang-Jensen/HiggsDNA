"""
to write the higgsdna json file with reading the files from local dir
"""
import os
import json
import sys
import argparse
import json
parser = argparse.ArgumentParser()
parser.add_argument('--EOSdir', type=str,
                    default='/eos/user/z/zhenxuan',
                    required=False,
                    help='''
                    ''')
args = parser.parse_args()
##############json file template
json_file_template_1 = '''
"{tag_id}":
    {{ "xs":0.001,
    "bf":1.0,
    "files":{{
    "2017":
     '''
json_file_template_2 = '''
    }
    },
'''
#########read dir files info
import subprocess
with open('/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_25/src/das_samples_UL17_R_gghh.json') as f:
# with open('/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_25/src/das_samples_UL17_data.json') as f:
	data = json.load(f)
json_file_name = "hhww_cusNANO_local"
with open(json_file_name + ".json","a") as fout:
    # fout.write("{")
    for i in range(len(data)):
        fout.write(json_file_template_1.format(   tag_id = str(list(data.keys())[i])
                                                ))
        # fout.write(json_file_template_1)
        fout.write("[")
        root_file_list = []
        samples_find = [j.split(".root")[0].split('-')[-1] for j in data[list(data.keys())[i]]]
        # print(list(data.keys())[i])
        for j in range(len(samples_find)):
            # print(samples_find[j])
            shell_cmd = 'ls /eos/user/z/zhenxuan/customized_NanoAOD/UL17/' + list(data.keys())[i] + '/*' + samples_find[j] + '*'
            # shell_cmd = 'ls /eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/UL17/' + list(data.keys())[i] + '/*' + samples_find[j] + '*'
            return_cmd = subprocess.run(shell_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8',shell=True)
            if return_cmd.returncode == 0:
                if (j != len(samples_find)-1):
                    fout.write('"' + str(return_cmd.stdout).split('\n')[0] + '"' + ',' +'\n' )
                else:
                    fout.write('"' + str(return_cmd.stdout).split('\n')[0] + '"'  +'\n' ) # not add , in the final one

                # print(str(return_cmd.stdout))
                # root_file_list.append(str(return_cmd.stdout))
            else:
                print("file can't find: %s, tag is : %s" %(return_cmd.stdout,samples_find[j]))
        fout.write("]")
        fout.write(json_file_template_2)
#########write to json
    # with open(json_file_name + ".json","a") as fout:
    #     print("dir name:",str(list(data.keys())[i]))
    #     fout.write(json_file_template.format(   tag_id = str(list(data.keys())[i]),
    #                                             samples_list = str(root_file_list).replace("'",'"')
    #                                             ))

# with open(json_file_name + ".json","a") as fout:
#     fout.write("}")
                                                