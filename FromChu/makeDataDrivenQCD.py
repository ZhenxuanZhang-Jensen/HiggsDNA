import ROOT
from array import array
GJets=ROOT.TFile.Open("Gjet_fake.root")
GJets_tree=GJets.Get("tree")
cuts="(Diphoton_lead_pt_mgg>0.35)*(Diphoton_sublead_pt_mgg>0.25)*(CMS_hgg_mass<115 || CMS_hgg_mass>135)"
h_minphotonID=ROOT.TH1F("h_minphotonID_gjet","h_minphotonID_gjet",19,-0.9,1);
GJets_tree.Project("h_minphotonID_gjet","minID","weight*41.5*"+cuts)
photonIDPDF_fake=ROOT.TF1("photonIDPDF_fake","pol7",-0.9,1.)
h_minphotonID.Fit(photonIDPDF_fake,"R")
c1=ROOT.TCanvas("c1","c1",600,800)
h_minphotonID.Draw("E1")
c1.SaveAs("test.png")

Data=ROOT.TFile.Open("Data.root")
Data_tree=Data.Get("tree")
nevents=Data_tree.GetEntries()
# new_weight = array('f', [0])
# print(new_weight)
# out = ROOT.TFile.Open("QCD.root", "RECREATE")
new_weight=-999
weights=[]
# newtree = Data_tree.CopyTree("1")
# _newWeight = newtree.Branch(
                        # "new_weight", new_weight, "new_weight/F")
minID=[]
maxID=[]
hasMaxLead=[]
originalminID=[]
print(nevents)
for i in range(0,nevents):
    Data_tree.GetEntry(i)
    # weights.append(1)
    if(Data_tree.Leading_Photon_MVA < Data_tree.Subleading_Photon_MVA):
        hasleadIDmin=True
        original_Photon_MVA_min = Data_tree.Leading_Photon_MVA
        Photon_MVA_max = Data_tree.Subleading_Photon_MVA
    else:
        hasleadIDmin=False
        original_Photon_MVA_min = Data_tree.Subleading_Photon_MVA
        Photon_MVA_max = Data_tree.Leading_Photon_MVA
    originalminID.append(original_Photon_MVA_min)
    maxID.append(Photon_MVA_max)
    # weights.append(1)
    if(not (original_Photon_MVA_min<-0.7 and Photon_MVA_max>-0.7)):
        new_weight=-999
        minID.append(-999)
        hasMaxLead.append(-999)
    else:
        if(hasleadIDmin):
            hasMaxLead.append(0)
            Leading_Photon_MVA=photonIDPDF_fake.GetRandom(-0.7,Photon_MVA_max)
            PhotonID_min=Leading_Photon_MVA
        else:
            Subleading_Photon_MVA=photonIDPDF_fake.GetRandom(-0.7,Photon_MVA_max)
            PhotonID_min=Subleading_Photon_MVA
            hasMaxLead.append(1)
        minID.append(PhotonID_min)
        new_weight = photonIDPDF_fake.Integral(-0.7,Photon_MVA_max) / photonIDPDF_fake.Integral(-0.9,-0.7);
    weights.append(new_weight)
    if(i%100000==0):
        print("Read entry:",i,new_weight)
    # _newWeight.Fill()
    # newtree.Write()
# out.Close()
# from root_numpy import tree2array
# newArr=tree2array(newtree)
print(sum(weights))
d={"new_weight":weights,"minID":minID,"maxID":maxID,"originalminID":originalminID,"hasMaxLead":hasMaxLead}
import pandas
dataframe=pandas.DataFrame(d) 
print(dataframe)
print(dataframe.shape[0])
dataframe.to_csv("data_weight.csv")