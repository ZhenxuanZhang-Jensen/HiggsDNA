import awkward as ak
import numpy as np
import uproot 

signal_file=ak.from_parquet("/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/WWgg/SL/SL3000/UL17_R_gghh_SL_M-3000_2017/merged_nominal.parquet")
data_file=ak.from_parquet("/eos/user/s/shsong/combined_WWgg/parquet/data_latest/merged_nominal.parquet")
sigcut=(data_file.Diphoton_mass>135.)|(data_file.Diphoton_mass<115.)
datasideband=data_file[sigcut]
datasideband_fatjet1_score=ak.unflatten(datasideband["fatjet_1_Hqqqq_vsQCDTop"],counts=1)
datasideband_fatjet2_score=ak.unflatten(datasideband["fatjet_2_Hqqqq_vsQCDTop"],counts=1)
datasideband_fatjet3_score=ak.unflatten(datasideband["fatjet_3_Hqqqq_vsQCDTop"],counts=1)
datasideband_fatjet_score=ak.concatenate([datasideband_fatjet1_score,datasideband_fatjet2_score,datasideband_fatjet3_score],axis=1)
datasideband_fatjet_score = datasideband_fatjet_score[ak.argsort(datasideband_fatjet_score,ascending=False,axis=-1)]
signal_fatjet1_score=ak.unflatten(signal_file["fatjet_1_Hqqqq_vsQCDTop"],counts=1)
signal_fatjet2_score=ak.unflatten(signal_file["fatjet_2_Hqqqq_vsQCDTop"],counts=1)
signal_fatjet3_score=ak.unflatten(signal_file["fatjet_3_Hqqqq_vsQCDTop"],counts=1)
signal_fatjet_score=ak.concatenate([signal_fatjet1_score,signal_fatjet2_score,signal_fatjet3_score],axis=1)
signal=ak.zip({"H_3q4q":signal_fatjet_score,"weight":signal_file["weight_central"] })
signal_fatjet_score = signal_fatjet_score[ak.argsort(signal_fatjet_score,axis=-1,ascending=False)]
signal=signal[ak.argsort(signal.H_3q4q,axis=-1,ascending=False)]
def calculate_significance(cut,signal,datasideband):
    s=ak.sum(signal.weight[signal.H_3q4q>cut])
    b=len(datasideband_fatjet_score[ak.num(datasideband_fatjet_score[datasideband_fatjet_score>cut])>0])
    return s / np.sqrt(s+b)
significance=[]
for i in range(1, 11):
    num = i / 10.0
    sob=calculate_significance(num,signal,datasideband_fatjet_score)
    significance.append(sob)

#plot the scatter plot
def plot_scatter_significance():
    size = 20  
    color = 'blue'  

    plt.scatter(cut, significance, s=size, c=color)
    plt.xticks(fontsize=12)  
    plt.yticks(fontsize=12)  
    plt.xlabel('H_3q4q score',fontsize=12)
    plt.ylabel('s/b',fontsize=12)
    plt.title('Significance', fontsize=14)
    plt.grid(True)
    plt.show()
