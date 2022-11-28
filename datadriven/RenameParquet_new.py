from email.utils import decode_rfc2231
import sys
# import sys
import awkward as ak
# sys.path.append("/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/py2-xgboost/0.72/lib/python2.7/site-packages/")
# sys.path.append("/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/py2-scikit-learn/0.19.1/lib/python2.7/site-packages/")
# sys.path.append("/afs/cern.ch/work/c/chuw/.higgsenv/lib/python3.9/site-packages/")
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# import xgboost as xgb
# from sklearn.model_selection import train_test_split
from random import choice
import parquet_to_root
from math import *
 
def get2BodyMass(data,pt1,eta1,phi1,m1,pt2,eta2,phi2,m2,name):
    
    px1 = data[pt1]*np.cos(data[phi1])
    py1 = data[pt1]*np.sin(data[phi1])
    pz1 = data[pt1]*np.sinh(data[eta1])
    px2 = data[pt2]*np.cos(data[phi2])
    py2 = data[pt2]*np.sin(data[phi2])
    pz2 = data[pt2]*np.sinh(data[eta2])
    E1=np.sqrt(px1*px1+py1*py1+pz1*pz1+data[m1]*data[m1])
    E2=np.sqrt(px2*px2+py2*py2+pz2*pz2+data[m2]*data[m2])
#     if (px2*px2+py2*py2+pz2*pz2-data[m2]*data[m2] < 0):
#         print("E1 is negative")
    mass=np.sqrt((E1+E2)*(E1+E2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2))
    data[name]=mass
    data["mStar"]=data[name]-data[m1]-data[m2]+250
    data["Dipho_pt_mStar"]=data["Diphoton_pt"]/data["mStar"]
    data["FatH4qJet_pt_mStar"]=data["fatjet_H_pt"]/data["mStar"]
#     return mass
#     E1=sqrt(px1*px1+py1*py1+pz1*pz1-m1*m1)
#     E2=sqrt(px2*px2+py2*py2+pz2*pz2-m2*m2)
# #     v1 = TLorentzVector(px1, py1, pz1, E1)
#     px2 = pt2*cos(phi2)
#     py2 = pt2*sin(phi2)
#     pz2 = pt2*sinh(eta2)
#     v2 = TLorentzVector(px2, py2, pz2, E2)
#     v3 = v1+v2
#     return v3.M()
def load_parquet(fname,isSignal:bool,Xmass,Type):
    pq=ak.from_parquet(fname)
    pq=pq[pq.category>0]
    df=ak.to_pandas(pq)
    print("before selection")
    print(df.shape[0])
    print(df.shape[1])
    if isSignal==True:
        df["label"]=1
    else:
        df["label"]=0
    if Type=="ggH":
        df["singleH"]=1
    elif Type=="ttH":
        df["singleH"]=2
    elif Type=="VH":
        df["singleH"]=3
    elif Type=="VBFH":
        df["singleH"]=4
    elif Type=="THQ":
        df["singleH"]=5
    else:
        df["singleH"]=0
    df["Xmass"]=Xmass
    df["Type"]=Type
    # df=df.query("category>0")
    if ("QCD" not in Type):
        GetMinMaxID(df)
    df["Dipho_pt_mgg"]=df.Diphoton_pt/df.Diphoton_mass
    get2BodyMass(df,"Diphoton_pt","Diphoton_eta","Diphoton_phi","Diphoton_mass","fatjet_H_pt","fatjet_H_eta","fatjet_H_phi","fatjet_H_mass","WWggMass")
    return df

def get_mass(x):
    Xmass=125
    mass_list=[1000,1250,1500,2000,2500,3000]
    if x>999:
        Xmass = x 
    else:
        Xmass=choice(mass_list)
    return Xmass
def getMinID(leadID,subleadID):
    if leadID<subleadID:
        return leadID
    else:
        return subleadID
def getMinIDFlav(leadID,subleadID,LeadPhoton_genPartFlav,SubleadPhoton_genPartFlav):
    if leadID<subleadID:
        return LeadPhoton_genPartFlav  
    else:
        return SubleadPhoton_genPartFlav
    
def getMaxID(leadID,subleadID):
    if leadID>=subleadID:
        return leadID
    else:
        return subleadID 
def reMapLeadPhoID(hasMaxLead,maxID,minID):
    if hasMaxLead:
        return maxID
    else:
        return minID
def reMapSubPhoID(hasMaxLead,maxID,minID):
    if hasMaxLead:
        return minID
    else:
        return maxID  
    
def reMapID(df):
     df["Leading_Photon_MVA"] = df.apply(lambda x: reMapLeadPhoID(x["hasMaxLead"],x["maxID"],x["minID"]),axis=1)
     df["Subleading_Photon_MVA"] = df.apply(lambda x: reMapSubPhoID(x["hasMaxLead"],x["maxID"],x["minID"]),axis=1)    
def GetMinMaxID(df):
    print("GetMinMaxID")
    # clum = df.shape[1]
    if type == "Data":
        df["minID"] = df.apply(lambda x: getMinID(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"]),axis=1)
        df["maxID"] = df.apply(lambda x: getMaxID(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"]),axis=1)
    else:
        df["minID"] = df.apply(lambda x: getMinID(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"]),axis=1)
        df["minID_genPartFlav"] = df.apply(lambda x: getMinIDFlav(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"],x["LeadPhoton_genPartFlav"],x["SubleadPhoton_genPartFlav"]),axis=1)
        # df1 = df.apply(lambda x: getMinID(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"]),axis=1)

    # df.insert(clum,'minID',df1)    
        df["maxID"] = df.apply(lambda x: getMaxID(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"]),axis=1)
    # df2 = df.apply(lambda x: getMaxID(x["LeadPhoton_mvaID"],x["SubleadPhoton_mvaID"]),axis=1)
    # df.insert(clum+1,'maxID',df2)
def RenameDf(df,type):
    if type == "Data":
        df=df.rename({"Diphoton_mass":"CMS_hgg_mass", 
                       "LeadPhoton_mvaID":"Leading_Photon_MVA",
                       "SubleadPhoton_mvaID":"Subleading_Photon_MVA",
                       "weight_central":"weight"
                      }, axis='columns')
    elif type == "BKG":
        df=df.rename({
                       "Diphoton_mass":"CMS_hgg_mass", 
                       "LeadPhoton_mvaID":"Leading_Photon_MVA",
                       "SubleadPhoton_mvaID":"Subleading_Photon_MVA",
                       "weight_central_no_lumi":"weight"
                      }, axis='columns')
#     print(df.columns)
    return df
def RenameDfQCD(df):
     df=df.rename({"Diphoton_mass":"CMS_hgg_mass", 
                       "LeadPhoton_mvaID":"Leading_Photon_MVA_old",
                        "SubleadPhoton_mvaID":"Subleading_Photon_MVA_old"
                       }, axis='columns')
     print(df.columns)
     return df
fname=sys.argv[1]
type = sys.argv[2]
print("loading files")
df=load_parquet(fname,False,0,type)
df_new=RenameDf(df,type)
df_new.to_parquet("/eos/user/s/shsong/hhwwggFH_parquet/Datadriven"+sys.argv[3]+".parquet")
# df_new.to_parquet(fname.split(".")[0]+"_new_2017.parquet")
parquet_to_root("/eos/user/s/shsong/hhwwggFH_parquet/Datadriven"+sys.argv[3],"/eos/user/s/shsong//eos/user/s/shsong/hhwwggFH_root/Datadriven"+sys.argv[3]+".root",treename=sys.argv[4],verbose=False)