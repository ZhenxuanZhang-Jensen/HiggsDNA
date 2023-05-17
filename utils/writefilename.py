import os
def file_name_walk(file_dir):
    dic_path="/eos/user/z/zhenxuan/wwyy/files/"
    path=dic_path+'json.txt'
    file=open(path,'w')
    for root,dirs,files in os.walk(file_dir):
        if "/eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/UL17/UL17_" in root:
            if "/eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/UL17/UL17_data" in root:
                continue
            else:
                sample_list_name=root.split("/")[-1]
                
                file.write('"'+sample_list_name+'":\n{"xs":0.001,\n"bf":1.0,\n"fpo":5,\n"files":{\n"2017":\n[')
                for root_file in files:
                    file.write('"/eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/UL17/'+sample_list_name+'/'+root_file+'",\n')
                file.write(']\n}\n},\n')
        else:
            continue

    


file_name_walk("/eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/UL17/")