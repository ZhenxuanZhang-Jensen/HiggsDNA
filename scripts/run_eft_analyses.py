import json
import os

eft_config = "/home/users/smay/HiggsDNA/metadata/analysis/hhggtautau/hh_ggtautau_sr_ul_eft.json"
benchmarks = [1,2,3,4,5,6,7,8,9,10,11,12]

bdts = [
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_1.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_2.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_3.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_4.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_5.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_6.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_7.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_8.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_9.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_10.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_11.xgb",
    "/home/users/fsetti/HHggTauTau/HggAnalysisDev/MVAs/output/eft_node_12.xgb" 
]

srs = [
    [  0.912695,  0.9830 ],
    [  0.994044,  0.9995 ],
    [  0.894935,  0.9857 ],
    [  0.909370,  0.9887],
    [  0.861141,  0.9816],
    [  0.828323,  0.9826],
    [  0.947518,  0.9875],
    [  0.479675,  0.9865],
    [  0.920841,  0.9844],
    [  0.924548,  0.9860],
    [  0.945699,  0.9772],
    [  0.957525,  0.9861]
]

def write_config(eft_config, n):
    with open(eft_config, "r") as f_in:
        lines = f_in.readlines()

    for idx, x in enumerate(lines):
        if 'EFT_BDT_FILE' in x:
            lines[idx] = x.replace('EFT_BDT_FILE', '"%s"' % bdts[n-1])
        if 'EFT_BDT_CUTS' in x:
            lines[idx] = x.replace('EFT_BDT_CUTS', str(srs[n-1]))
        if 'EFT_SAMPLE' in x:
            lines[idx] = x.replace('EFT_SAMPLE', '"HH_ggTauTau_EFT_node_%d"' % n)

    f_eft = eft_config.replace(".json", "_%d.json" % n)
    with open(f_eft, "w") as f_out:
        for x in lines:
            f_out.write(x)         

    return f_eft


for n in benchmarks:
    config = write_config(eft_config, n)
    command = "/bin/nice -n 19 python scripts/run_analysis.py --config %s --output_dir 'ggtautau_sr_eft_%d_21Apr2022' --log-level 'INFO' --batch_system 'local' --merge_outputs --n_cores 12" % (config.replace("/home/users/smay/HiggsDNA/", ""), n)
    print(command)
    os.system(command)
