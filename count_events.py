import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
def read_files(file_name):
    events = ak.from_parquet(file_name)
    return events
def read_files_uproot(file_name):
    events = uproot.open(file_name)
    return events
def count_events(events_name, var_name, mask):
    counts, bins, patches = plt.hist(events_name[var_name][mask], weights = events_name['weight_central'][mask])
    total_events = np.sum(counts)
    return total_events
def plot_var(*var_arr,mask, weight, legend_name,bins):
    for i in range(len(var_arr)):
        plt.hist(var_arr[i][mask], weights = weight[mask], histtype="step", density=True, label=legend_name,bins=bins)
        plt.legend()
def delta_R_O(phi1, phi2):
    delta_phi = phi2 - phi1
    delta_eta = 0
    return np.sqrt(delta_phi**2 + delta_eta**2)

def delta_R(phi1, eta1, phi2, eta2):
    delta_phi = phi2 - phi1
    delta_eta = eta2 - eta1
    return np.sqrt(delta_phi**2 + delta_eta**2)

import multiprocessing
import awkward as ak

def read_and_count_parquet(filename, queue):
    events = read_files(filename)
    # 读取parquet文件
    cat1_data = count_events(events_name = events, var_name='weight_central', mask=((events['weight_central']>0) & (events['category']==1)))
    cat2_data = count_events(events_name = events, var_name='weight_central', mask=((events['weight_central']>0) & (events['category']==2)))
    cat3_data = count_events(events_name = events, var_name='weight_central', mask=((events['weight_central']>0) & (events['category']==3)))

    # 将统计结果放入进程间通信队列中
    queue.put((cat1_data,cat2_data,cat3_data))

if __name__ == '__main__':
    filenames = ["/eos/user/s/shsong/combined_WWgg/bkg/UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_DiPhotonJetsBox_M40_80_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_GluGluHToGG_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_TTGG_0Jets_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_VHToGG_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_WWG_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_ttHJetToGG_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_TTGJets_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_VBFHToGG_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_W3JetsToLNu_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_WWTo1L1Nu2Q_4f_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_ttWJets_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_DiPhotonJetsBox_MGG_80toInf_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_QCD_Pt-30to40_MGG-80toInf_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_QCD_Pt-30toInf_MGG-40to80_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_QCD_Pt-40ToInf_MGG-80ToInf_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_TTJets_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_TTToHadronic_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_W1JetsToLNu_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_W2JetsToLNu_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_W4JetsToLNu_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/bkg/UL17_WGJJToLNu_EWKnotop_QCD_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/UL17_dataB_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/UL17_dataC_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/UL17_dataD_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/UL17_dataE_2017/merged_nominal.parquet","/eos/user/s/shsong/combined_WWgg/UL17_dataF_2017/merged_nominal.parquet"]
    queues = [multiprocessing.Queue() for _ in filenames]
    processes = [multiprocessing.Process(target=read_and_count_parquet, args=(filename, queue)) for filename, queue in zip(filenames, queues)]

    # 启动进程
    for process in processes:
        process.start()

    # 获取统计结果
    results = []
    for queue in queues:
        result = queue.get()
        results.append(result)

    # 等待进程结束
    for process in processes:
        process.join()

    # 打印统计结果
    for filename, cat1_data,cat2_data,cat3_data in zip(filenames, results):
        print(f'{filename} 包含 {cat1_data} 个cat1, {cat2_data} 个cat2, {cat3_data} 个 cat3 ')
