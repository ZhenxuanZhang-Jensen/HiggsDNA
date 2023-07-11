import uproot
import awkward as ak
import numpy as np
import json

def get_weighted_events_in_cut_range(Hscores, weights, lower_cut, upper_cut):
    return weights[(Hscores > lower_cut) & (Hscores <= upper_cut)]
def compute_significance(signal_weight_sum, bkg_weight_sum):
    s = signal_weight_sum
    b = bkg_weight_sum
    significance=np.sqrt(2*((s+b)*np.log(1+s/b)-s))

    return significance
def get_W_events_out_of_cut_range(Hscores,Wscores, weights, lower_cut, W_workingpoints):
    return weights [(Hscores < lower_cut) & (Wscores > W_workingpoints)]

def calculate_3_cat_significance(sig_fatjet_W_cat1,sig_fatjet_W_cat2,bkg_fatjet_W_cat1,bkg_fatjet_W_cat2,total_significance):
    s_Wcat1=ak.sum(sig_fatjet_W_cat1.weight[:,0])
    b_Wcat1=len(bkg_fatjet_W_cat1)
    if b_Wcat1!=0:
        significance_Wcat1=np.sqrt(2*((s_Wcat1+b_Wcat1)*np.log(1+s_Wcat1/b_Wcat1)-s_Wcat1))
    else:
        significance_Wcat1=999
    s_Wcat2=ak.sum(sig_fatjet_W_cat2.weight)
    b_Wcat2=len(bkg_fatjet_W_cat2)    
    if b_Wcat2!=0:
        
        significance_Wcat2=np.sqrt(2*((s_Wcat2+b_Wcat2)*np.log(1+s_Wcat2/b_Wcat2)-s_Wcat2))
    else:
        significance_Wcat2=999
    W_significance=np.sqrt(significance_Wcat1*significance_Wcat1+significance_Wcat2*significance_Wcat2)
    significance_Hcat=total_significance
    sum_significance=np.sqrt(significance_Wcat1*significance_Wcat1+significance_Wcat2*significance_Wcat2+significance_Hcat*significance_Hcat)
    return sum_significance,W_significance

signal=ak.from_parquet("/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggFH/UL17_R_gghh_M-3000_2017/merged_nominal.parquet")
signal_H_scores_dict = {}
signal_H_weights_dict = {}
signal_H_weights_dict["M3000"]=signal["weight_central"][signal["fatjet_H_Hqqqq_vsQCDTop"]!=-999]
signal_H_scores_dict['M3000']=signal["fatjet_H_Hqqqq_vsQCDTop"][signal["fatjet_H_Hqqqq_vsQCDTop"]!=-999]
data=ak.from_parquet("/eos/user/s/shsong/combined_WWgg/parquet/data_v1/merged_nominal.parquet")
mask = ((data["Diphoton_mass"] < 115) | (data["Diphoton_mass"] > 135))
data_sideband = data[mask]
data_sideband_H_scores = np.sort(data_sideband["fatjet_H_Hqqqq_qqlv_vsQCDTop"])[::-1][np.sort(data_sideband["fatjet_H_Hqqqq_qqlv_vsQCDTop"])[::-1]!=-999]


idx = 0
lower_cut = data_sideband_H_scores[idx]
upper_cut = 1
prev_total_significance = 0
improvement = 1
signal_mass_points = ["M3000"]
categories_dict = {mp: [] for mp in signal_mass_points}
total_significance = 0
H_total_significance=[]
H_significance=[]
significance_3cats=[]
W_total_significance=[]
while lower_cut >= 0. and idx < len(data_sideband_H_scores):
    idx += 10
    if idx < len(data_sideband_H_scores):
        lower_cut = data_sideband_H_scores[idx]
    else:
        lower_cut = 0.
    
    for cat_lower_cut, cat_upper_cut in categories_dict['M3000'] + [(lower_cut, upper_cut)]:
        signal_weights_in_cut = get_weighted_events_in_cut_range(signal_H_scores_dict['M3000'], signal_H_weights_dict['M3000'], cat_lower_cut, cat_upper_cut)
        no_H_fatjetsig=signal[(signal["nGoodAK8jets"]>0)&(signal['fatjet_H_Hqqqq_vsQCDTop']<cat_lower_cut)]
        datasideband=data_sideband_H_scores[(data_sideband_H_scores>cat_lower_cut)&(data_sideband_H_scores<cat_upper_cut)]
        significance=compute_significance(ak.sum(signal_weights_in_cut),len(datasideband))
        total_significance = np.sqrt(total_significance**2 + significance**2)
    improvement = (total_significance - prev_total_significance) / prev_total_significance if prev_total_significance != 0 else 1
    if improvement > 0.01:
        H_significance.append(significance)
        categories_dict['M3000'].append((lower_cut, upper_cut))
        upper_cut = lower_cut
        prev_total_significance = total_significance
        H_total_significance.append(total_significance)
        
        no_H_fatjet_sig=signal[(signal["nGoodAK8jets"]>0)&(signal['fatjet_H_Hqqqq_vsQCDTop']<cat_lower_cut)]
        sig_fatjet1_score = ak.unflatten(no_H_fatjet_sig["fatjet_1_particleNet_WvsQCD"],counts=1)
        sig_fatjet2_score = ak.unflatten(no_H_fatjet_sig["fatjet_2_particleNet_WvsQCD"],counts=1)
        sig_fatjet3_score = ak.unflatten(no_H_fatjet_sig["fatjet_3_particleNet_WvsQCD"],counts=1)
        sig_fatjets_W_score=ak.concatenate([sig_fatjet1_score,sig_fatjet2_score,sig_fatjet3_score],axis=1)
        sig_fatjets_W_score=sig_fatjets_W_score[ak.argsort(sig_fatjets_W_score,ascending = False, axis = 1)][:,:2]
        sig_fatjet_W=ak.zip({"W":sig_fatjets_W_score,"nAK4":no_H_fatjet_sig['nGoodAK4jets'],"weight":no_H_fatjet_sig['weight_central']})
        sig_fatjet_W_cat1=sig_fatjet_W[ak.num(sig_fatjet_W[sig_fatjet_W.W>0.71])==2]
        sig_fatjet_W_cat2=sig_fatjet_W[(ak.num(sig_fatjet_W[sig_fatjet_W.W>0.71])==1)][:,0]
        sig_fatjet_W_cat2=sig_fatjet_W_cat2[sig_fatjet_W_cat2.nAK4>=2]
        
        no_H_fatjet_bkg=data_sideband[(data_sideband["nGoodAK8jets"]>0)&(data_sideband['fatjet_H_Hqqqq_vsQCDTop']<cat_lower_cut)]
        bkg_fatjet1_score = ak.unflatten(no_H_fatjet_bkg["fatjet_1_particleNet_WvsQCD"],counts=1)
        bkg_fatjet2_score = ak.unflatten(no_H_fatjet_bkg["fatjet_2_particleNet_WvsQCD"],counts=1)
        bkg_fatjet3_score = ak.unflatten(no_H_fatjet_bkg["fatjet_3_particleNet_WvsQCD"],counts=1)
        bkg_fatjets_W_score=ak.concatenate([bkg_fatjet1_score,bkg_fatjet2_score,bkg_fatjet3_score],axis=1)
        bkg_fatjets_W_score=bkg_fatjets_W_score[ak.argsort(bkg_fatjets_W_score,ascending = False, axis = 1)][:,:2]
        bkg_fatjet_W=ak.zip({"W":bkg_fatjets_W_score,"nAK4":no_H_fatjet_bkg['nGoodAK4jets'],"weight":no_H_fatjet_bkg['weight_central']})
        bkg_fatjet_W_cat1=bkg_fatjet_W[ak.num(bkg_fatjet_W[bkg_fatjet_W.W>0.71])==2]
        bkg_fatjet_W_cat2=bkg_fatjet_W[(ak.num(bkg_fatjet_W[bkg_fatjet_W.W>0.71])==1)][:,0]
        bkg_fatjet_W_cat2=bkg_fatjet_W_cat2[bkg_fatjet_W_cat2.nAK4>=2]
        sum_significance,W_significance=calculate_3_cat_significance(sig_fatjet_W_cat1,sig_fatjet_W_cat2,bkg_fatjet_W_cat1,bkg_fatjet_W_cat2,total_significance)
        significance_3cats.append(sum_significance)
        W_total_significance.append(W_significance)