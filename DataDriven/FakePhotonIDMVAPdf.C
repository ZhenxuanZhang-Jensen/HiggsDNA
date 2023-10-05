//#include "setTDRStyle.C"
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
//#include "TLorentzVector.h"

using namespace std;
// This is the input file name, which is the outputGJet bkg file produced by flashgg framework
const string InputFileName = "/eos/user/z/zhenxuan/hhwwgg_root/bkgs_UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf_2017cat1.root";
// This is the input tree name
const string InputTree = "cat1";
// This is the output root file name
const string OutputFileName = "/eos/user/z/zhenxuan/wwyy/files/pf_FakeIdMvaPdf_cat1.root";
// Here have the definition of global values
const float BinXLow = -0.95; //minIDCut
const float BinXHig = 1.0;
const int BinTotal=10;
// Here begin the main function
void FakePhotonIDMVAPdf(){
	// Here use the Hgg paper style for plotting
  gROOT->ProcessLine(".x hggPaperStyle.C");
  gStyle->SetOptStat("iemr");
  //gStyle->SetOptFit(111);
  gStyle->SetOptFit(0);
	// Here define the new chain for read and write the tree vars
  TChain *fChain=new TChain(InputTree.c_str());
  fChain->Add(InputFileName.c_str());
  TFile *newfile = new TFile(OutputFileName.c_str(),"recreate");
  TH1D *fakeIDMVA = new TH1D("fakeIDMVA", "", BinTotal, BinXLow, BinXHig);
  Double_t         weight_central;
  Float_t         LeadPhoton_mvaID;
  Float_t         SubleadPhoton_mvaID;
  Float_t         lead_genMatchType;
  Float_t         sublead_genMatchType;
  fChain->SetBranchAddress("weight_central", &weight_central);
  fChain->SetBranchAddress("LeadPhoton_mvaID", &LeadPhoton_mvaID);
  fChain->SetBranchAddress("SubleadPhoton_mvaID", &SubleadPhoton_mvaID);
  Float_t nentries = fChain->GetEntries();
  for (Float_t i=0;i<nentries;i++) {
    fChain->GetEntry(i);
  	{
      if(SubleadPhoton_mvaID>=BinXLow) fakeIDMVA->Fill(SubleadPhoton_mvaID, weight_central);
      if(LeadPhoton_mvaID>=BinXLow)  fakeIDMVA->Fill(LeadPhoton_mvaID, weight_central);
    } 
  }
  cout<<"JTao: out of the loop"<<endl;
  TF1 *f1 = new TF1("f1","pol7(0)",BinXLow, BinXHig);
  //=======
  TCanvas *c1 = new TCanvas("reconstruction1","reconstruction1");
  c1->SetFillColor(0);
  
  fakeIDMVA->SetLineColor(4);
  fakeIDMVA->SetFillColor(0);
  fakeIDMVA->SetLineStyle(1);
  fakeIDMVA->SetLineWidth(2);
  fakeIDMVA->GetXaxis()->SetTitle("ID MVA score of fake photons");
  double WidthBin=(BinXHig-BinXLow)/BinTotal;
  string PreTitleY( Form("Events / %.2g ",WidthBin) );
  string TitleY =  PreTitleY + "";
  fakeIDMVA->GetYaxis()->SetTitle(TitleY.c_str());
  //fakeIDMVA->GetYaxis()->SetTitle("Events / 0.01");
  fakeIDMVA->SetTitleSize(0.05,"XY");
  fakeIDMVA->SetTitleOffset(1.3, "Y");
  fakeIDMVA->SetMarkerColor(4);
  fakeIDMVA->SetMarkerSize(0.6);
  fakeIDMVA->SetMarkerStyle(20);
  
  //-------------
  c1->cd();
  //gPad->SetLogy(1);
  fakeIDMVA->Draw("PE1");
  f1->SetLineColor(2);
  f1->SetLineWidth(2);
  fakeIDMVA->Fit("f1");
  cout<<"JTao: N_hist = "<<fakeIDMVA->Integral()<<", integral from tf1 = "<<f1->Integral(BinXLow, BinXHig)<<endl;
  c1->Print("pf_FakeIdMva_FitPol7.pdf");
  c1->Clear();
  //=============
  fakeIDMVA->Write();
  f1->Write();
  newfile->Write();
  newfile->Close();
}