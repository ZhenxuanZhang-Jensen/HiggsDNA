#include "TMath.h"
//#include "TROOT.h"
#include <TH1D.h>
#include <TH1D.h>
#include <TH1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeafF.h>
#include <TChain.h>
#include <TFile.h>
#include "TSystem.h"
#include <TChain.h>
#include "TSystem.h"
#include <TString.h>
#include <iostream>
#include <vector>
#include <TPostScript.h>
#include <iostream>
#include <iomanip>  //for precision
#include <TTree.h>
#include <TBranch.h>
//#include "TLorentzVector.h"
#include "TMVA/Reader.h"

using namespace std;

const string InputFileName = "/afs/cern.ch/user/z/zhenxuan/HiggsDNA/DataDriven/New_cat1.root";
const string InputTree = "cat1";
const string OutputFileName = "/eos/user/z/zhenxuan/wwyy/root_files/DataDriven_data_cat1.root";

// const string InputFakePdfFileName = "/eos/user/s/shsong/combined_WWgg/datadriven/category2/GJet.root";
const string InputFakePdfFileName = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2016/preVFP/pf_FakeIdMvaPdf.root";
const string FakePdfProb = "f1";

const float minIDCut = -0.9; 
const float minIDSideband = -0.95; 
void DataDrivenQCD(){
  //input file
  TChain *fChain=new TChain(InputTree.c_str());
  fChain->Add(InputFileName.c_str());

  //fake prob.
  TFile *myfile = new TFile(InputFakePdfFileName.c_str());
  TF1 *fake_prob = (TF1*)myfile->Get("f1");
  float IntProb_sideband = fake_prob->Integral(minIDSideband, minIDCut);
  cout<<" JTao : Integral of tf1 = "<<IntProb_sideband<<" from "<<minIDSideband<<" to "<<minIDCut<<endl; 
  
  //new output file
  TFile *newfile = new TFile(OutputFileName.c_str(),"recreate");
  //TTree *newtree = fChain->CloneTree();
  TTree *newtree = fChain->CloneTree(0);
  cout<<" JTao: Nin="<<fChain->GetEntries()<<" Nclone="<<newtree->GetEntries()<<endl;
  //return;
  //----
  float minIDmvaRaw;
  float maxIDmvaRaw;
  float minIDmva;
  float leadmvaRaw;
  float subleadmvaRaw; 
  //float diphoMVANew;
  float diphoMVARaw;  
  float weight;
  newtree->Branch("minIDmvaRaw", &minIDmvaRaw, "minIDmvaRaw/F");
  newtree->Branch("maxIDmvaRaw", &maxIDmvaRaw, "maxIDmvaRaw/F");
  newtree->Branch("minIDmva", &minIDmva, "minIDmva/F");
  newtree->Branch("leadmvaRaw", &leadmvaRaw, "leadmvaRaw/F");
  newtree->Branch("subleadmvaRaw", &subleadmvaRaw, "subleadmvaRaw/F");



  //----
  Float_t         LeadPhoton_mvaID;
  Float_t         SubleadPhoton_mvaID;
  Float_t         Diphoton_minID;
  Float_t         Diphoton_maxID;
  fChain->SetBranchAddress("LeadPhoton_mvaID", &LeadPhoton_mvaID);
  fChain->SetBranchAddress("SubleadPhoton_mvaID", &SubleadPhoton_mvaID);
  fChain->SetBranchAddress("Diphoton_minID", &Diphoton_minID);
  fChain->SetBranchAddress("Diphoton_maxID", &Diphoton_maxID);

  float leadptom;
  float subleadptom;
  float leadeta;
  float subleadeta;
  float Diphoton_mass;
  float sigmawv;
  float CosPhi;
  float vtxprob;

  fChain->SetBranchAddress("Diphoton_mass", &Diphoton_mass);
  //---
  //----diphootn BDT
  float Diphoton_mass_=0.0;
  float leadmva_=0.0;
  float subleadmva_=0.0;
  


  //---
  Long64_t nentries = fChain->GetEntries();
  float weight_central;
  newtree->Branch("weight_central", &weight_central, "weight_central/F");
  for (Long64_t i=0;i<nentries;i++) {
//  for (Long64_t i=0;i<100;i++) {
    fChain->GetEntry(i);
    newtree->GetEntry(i);
    maxIDmvaRaw = max(LeadPhoton_mvaID,SubleadPhoton_mvaID);
    minIDmvaRaw = min(LeadPhoton_mvaID,SubleadPhoton_mvaID);
//    if(i%1000==0) cout<<" JTao: maxIDmvaRaw="<<maxIDmvaRaw<<" and minIDmvaRaw="<<minIDmvaRaw<<endl;
    if(maxIDmvaRaw <= minIDCut){
      cout<<"JTao: Warning maxIDmvaRaw "<<maxIDmvaRaw<<" <= minIDCut "<<minIDCut<<" with i = "<<i<<endl;    
      maxIDmvaRaw = minIDCut + 1.0e-6;
    }
    fake_prob->SetRange(minIDCut, maxIDmvaRaw);
    minIDmva = fake_prob->GetRandom(); 
/*
    if(maxIDmvaRaw <= minIDCut) {
      cout<<"JTao: Warning maxIDmvaRaw "<<maxIDmvaRaw<<" <= minIDCut "<<minIDCut<<" with i = "<<i<<endl;
      maxIDmvaRaw = minIDCut + 1.0e-6;
    }
*/
    weight_central = fake_prob->Integral(minIDCut, maxIDmvaRaw)/IntProb_sideband;
    if(i%100==0 ) cout<<" JTao:  minIDmva="<<minIDmva<<", weight_central="<<weight_central<<" with i = "<<i<<endl;
    //===
    leadmvaRaw = LeadPhoton_mvaID; subleadmvaRaw = SubleadPhoton_mvaID;
    LeadPhoton_mvaID = LeadPhoton_mvaID>SubleadPhoton_mvaID?LeadPhoton_mvaID:minIDmva;
    SubleadPhoton_mvaID = SubleadPhoton_mvaID>LeadPhoton_mvaID?SubleadPhoton_mvaID:minIDmva;
    Diphoton_maxID = LeadPhoton_mvaID>SubleadPhoton_mvaID?LeadPhoton_mvaID:SubleadPhoton_mvaID;
    Diphoton_minID = LeadPhoton_mvaID>SubleadPhoton_mvaID?SubleadPhoton_mvaID:LeadPhoton_mvaID;
    //==diphoton MVA===
    Diphoton_mass_=Diphoton_mass;
    leadmva_=LeadPhoton_mvaID;
    subleadmva_=SubleadPhoton_mvaID;   
    newtree->Fill();
  }
  cout<<" JTao: Nnew="<<newtree->GetEntries()<<endl;
  newtree->Write();
  delete newfile;
  // // Open input file
  // TFile* input_file = new TFile(OutputFileName.c_str(), "UPDATE");
  // if (!input_file || !input_file->IsOpen() || input_file->IsZombie()) {
  //     std::cerr << "Failed to open input file." << std::endl;
  //     return 1;
  // }
  // Replace weight_central with weight in the "tree" TTree
  // replace_weight_central(input_file, "cat1");

  // Close the input file
  // input_file->Close();



}
