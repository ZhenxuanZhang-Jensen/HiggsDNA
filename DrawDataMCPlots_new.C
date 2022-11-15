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

//==============
const int debug=1;

//signal only shape for comparison
// const TString InputSignalFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/XGboost_diphoMVA_UL18_newSig_newModel.root";
// const TString InputSignalFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/output_sig125_IncludeLumi.root";
// const TString TreeNameSig = "Sig125";
// const TString InputDataFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/New_UL2017data/MassUL2018_ETSS.root";
const TString InputDataFile = "/eos/user/s/shsong/hhwwggSL_root/data/2017data.root";
// const TString InputDataFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/NewDiphotonBDT_UL2018data_lessoverfitting.root";
// const TString InputDataFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/NewDiphotonBDT_UL2018data_lessoverfitting.root";
const TString TreeNameData = "Data_13TeV_2017";
//---
// const TString InputPPFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/New_pp.root";
// const TString InputPPFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/XGboost_diphoMVA_UL18_newSig.root";
const TString InputggFile = "/eos/user/s/shsong/hhwwggSL_root/bkg/DiphotonJetBox.root";
const TString TreeNamegg = "DiPhotonJetsBox";
const TString InputGJetsFile = "/eos/user/s/shsong/hhwwggSL_root/bkg/GJets.root";
const TString TreeNameGJets = "Gjet";
const TString InputWJets = "/eos/user/s/shsong/hhwwggSL_root/bkg/WJets.root";
const TString TreeNameWJets = "WJetsToLNu";
const TString InputTTFile = "/eos/user/s/shsong/hhwwggSL_root/bkg/TT.root";
const TString TreeNameTT = "TT";
const TString InputWGJJFile = "/eos/user/s/shsong/hhwwggSL_root/bkg/WGJJToLNu_2017.root";
const TString TreeNameWGJJ = "WGJJToLNu";
const TString InputttWFile = "/eos/user/s/shsong/hhwwggSL_root/bkg/ttW_2017.root";
const TString TreeNamettW = "ttW";

// const TString InputPPFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/XGboost_diphoMVA_UL18_newSig_newModel.root";
// const TString InputPPFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/.root";
// const TString TreeNamePP = "pp";
// const TString InputQCDFile = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/UL2018/MassUL2018_ETSS/XGboost_diphoMVA_UL18_newSig_newModel.root";//SF + sigmaE/E weights + Tao 2D pT weights
// onst TString TreeNameQCD = "DataDriven_QCD";
// const string OutputPlotDir = "DataMCComparisonPlots_mvaCuts";
const string OutputPlotDir = "/eos/user/s/shsong/DataMCComparisonPlots";
// const string OutputPlotDir = "DataMCComparisonPlots_old";
//const string OutputPlotDir = "DataMCComparisonPlotsBothEB";
//const string OutputPlotDir = "DataMCComparisonPlotsNotEBEB";

// const string Preselections="";
const string Preselections="(CMS_hgg_mass <= 115. || CMS_hgg_mass >= 135.)";
//const string Preselections="(CMS_hgg_mass <= 115. || CMS_hgg_mass >= 135.) && fabs(leadeta)<1.5 && fabs(subleadeta)<1.5";
//const string Preselections="(CMS_hgg_mass <= 115. || CMS_hgg_mass >= 135.) && (fabs(leadeta)>1.5 || fabs(subleadeta)>1.5)";
// const string MCWeight = "weight";
const string MCWeight = "weight";// "weight*41.5"; // add a new weight to scale diphoton bdt score >-0.85

// */
//====================



///////////////////////////////////////

void DrawMyPlots(string Object, string Selections,  string XTitle, string YUnit, string PlotName, int BinTotal, double BinXLow, double BinXHig, int LegendLR=0, int IfLogY=0, int IfStatErr=0){


  TCanvas *c1 = new TCanvas("reconstruction1","reconstruction1");
  c1->SetFillColor(0);

  TLegend *legend;

  if(IfStatErr==1){
    if(LegendLR==0) legend = new TLegend(0.65,0.67,0.89,0.9); 
    else if (LegendLR==-1) legend = new TLegend(0.45,0.67,0.69,0.9); //Middle
    else legend = new TLegend(0.2,0.67,0.44,0.9);
  }else{
    if(LegendLR==0) legend = new TLegend(0.65,0.7,0.89,0.9); //(0.7,0.75,0.9,0.9);
    else if (LegendLR==-1) legend = new TLegend(0.45,0.7,0.69,0.9); //Middle
    else legend = new TLegend(0.2,0.7,0.44,0.9);  
  }
  legend->SetFillColor(0);

  //====add root file
  TChain *Data_Tree=new TChain(TreeNameData);
  TChain *MCgg_Tree=new TChain(TreeNamegg);
  TChain *MCGJets_Tree=new TChain(TreeNameGJets);
  TChain *MCWJets_Tree=new TChain(TreeNameWJets);
  TChain *MCTT_Tree=new TChain(TreeNameTT);
  TChain *MCWGJJ_Tree=new TChain(TreeNameWGJJ);
  TChain *MCttW_Tree=new TChain(TreeNamettW);  


  //====================
  Data_Tree->Add(InputDataFile);
  MCgg_Tree->Add(InputggFile);
  MCGJets_Tree->Add(InputGJetsFile);
  MCWJets_Tree->Add(InputWJets);
  MCTT_Tree->Add(InputTTFile);
  MCWGJJ_Tree->Add(InputWGJJFile);
  MCttW_Tree->Add(InputttWFile);


  //=========entries================
  int entries_Data = Data_Tree->GetEntries();
  if(debug==1) cout <<"JTao: nEntries_Data = "<<entries_Data<<endl;

  //int entries_MC = MC_Tree->GetEntries();
  //if(debug==1) cout <<"JTao: nEntries_MC = "<<entries_MC<<endl;

  c1->cd();

  //=================
  char *myLimits= new char[100];
  sprintf(myLimits,"(%d,%f,%f)",BinTotal,BinXLow,BinXHig);
  TString Taolimits(myLimits);

  cout<<"JTao : selections -- "<<Selections<<endl;
  //====data=======
  TString variable_Data = Object + ">>Histo_Data_temp" + Taolimits;
  Data_Tree->Draw(variable_Data, Selections.c_str());
  TH1D *h_data = (TH1D*)gDirectory->Get("Histo_Data_temp");
  h_data->SetTitle("");
  c1->Clear();

  double Ntot_Data=h_data->Integral();
  if( debug==1 ) cout<<"JTao: N_Data= "<<Ntot_Data<<endl;

  Double_t scale_Data = 1.0/Ntot_Data;
  h_data->Sumw2();
  //h_data->Scale(scale_Data);

  //---MC---
  string MCSelections = MCWeight + "*(" + Selections + ")";
  const string Sig_weight = "weight"; // add a new weight to scale diphoton bdt score >-0.85
  cout<<"JTao : MC selections -- "<<MCSelections<<endl;
  TString variable_MCgg = Object + ">>Histo_MCgg_temp" + Taolimits;
  MCgg_Tree->Draw(variable_MCgg,MCSelections.c_str());
  TH1D *h_MCgg = (TH1D*)gDirectory->Get("Histo_MCgg_temp");
  c1->Clear();
  TString variable_MCGJets = Object + ">>Histo_MCGJets_temp" + Taolimits;
  MCGJets_Tree->Draw(variable_MCGJets,MCSelections.c_str());
  TH1D *h_MCGJets = (TH1D*)gDirectory->Get("Histo_MCGJets_temp");
  c1->Clear();
  TString variable_MCWJets = Object + ">>Histo_MCWJets_temp" + Taolimits;
  MCWJets_Tree->Draw(variable_MCWJets,MCSelections.c_str());
  TH1D *h_MCWJets = (TH1D*)gDirectory->Get("Histo_MCWJets_temp");
  c1->Clear();
  TString variable_MCTT = Object + ">>Histo_MCTT_temp" + Taolimits;
  MCTT_Tree->Draw(variable_MCTT,MCSelections.c_str());
  TH1D *h_MCTT = (TH1D*)gDirectory->Get("Histo_MCTT_temp");
  c1->Clear();
  TString variable_MCWGJJ = Object + ">>Histo_MCWGJJ_temp" + Taolimits;
  MCWGJJ_Tree->Draw(variable_MCWGJJ,MCSelections.c_str());
  TH1D *h_MCWGJJ = (TH1D*)gDirectory->Get("Histo_MCWGJJ_temp");
  c1->Clear();
  TString variable_MCttW = Object + ">>Histo_MCttW_temp" + Taolimits;
  MCttW_Tree->Draw(variable_MCttW,MCSelections.c_str());
  TH1D *h_MCttW = (TH1D*)gDirectory->Get("Histo_MCttW_temp");
  c1->Clear();
  // TString variable_MCpp = Object + ">>Histo_MCpp_temp" + Taolimits;
  // MCpp_Tree->Draw(variable_MCpp,MCSelections.c_str());
  // TH1D *h_MCpp = (TH1D*)gDirectory->Get("Histo_MCpp_temp");
  
  // signal shape only for comparsion:
//  TString variable_MCsig = Object + ">>Histo_MCsig_temp" + Taolimits;
  // MCsig_Tree->Draw(variable_MCsig, Sig_weight.c_str());
  // TH1D *h_MCsig = (TH1D*)gDirectory->Get("Histo_MCsig_temp");
  // h_MCsig->Scale(1000);
  // std::cout << "MCsig entries=" << h_MCsig->Integral() << std::endl;
  // c1->Clear();


  // TString variable_MCqcd = Object + ">>Histo_MCqcd_temp" + Taolimits;
  // MCqcd_Tree->Draw(variable_MCqcd, MCSelections.c_str());
  // TH1D *h_MCqcd = (TH1D*)gDirectory->Get("Histo_MCqcd_temp");
  // c1->Clear();
  std::cout << "MCgg entries=" << h_MCgg->Integral() << std::endl;
  std::cout << "MCGJets entries=" << h_MCGJets->Integral() << std::endl;
  std::cout << "MCWJets entries=" << h_MCWJets->Integral() << std::endl;
  std::cout << "MCTT entries=" << h_MCTT->Integral() << std::endl;
  std::cout << "MCWGJJ entries=" << h_MCWGJJ->Integral() << std::endl;
  std::cout << "MCttW entries=" << h_MCttW->Integral() << std::endl;
  h_MCgg->SetLineColor(50);
  h_MCgg->SetFillColor(50);
  // h_MCgg->SetFillStyle(3004);
  h_MCgg->SetLineStyle(1);
  h_MCgg->SetLineWidth(2);

  h_MCGJets->SetLineColor(38);
  h_MCGJets->SetFillColor(38);
  //h_MCGJets->SetFillStyle(3004);
  h_MCGJets->SetLineStyle(1);
  h_MCGJets->SetLineWidth(2);

  h_MCWJets->SetLineColor(30);
  h_MCWJets->SetFillColor(30);
  //hpf->SetFillStyle(3005);
  h_MCWJets->SetLineStyle(1);
  h_MCWJets->SetLineWidth(2);

  h_MCTT->SetLineColor(34);
  h_MCTT->SetFillColor(34);
  //hpf->SetFillStyle(3005);
  h_MCTT->SetLineStyle(1);
  h_MCTT->SetLineWidth(2);  

  h_MCWGJJ->SetLineColor(40);
  h_MCWGJJ->SetFillColor(40);
  //hpf->SetFillStyle(3005);
  h_MCWGJJ->SetLineStyle(1);
  h_MCWGJJ->SetLineWidth(2);

  h_MCttW->SetLineColor(60);
  h_MCttW->SetFillColor(60);
  //hpf->SetFillStyle(3005);
  h_MCttW->SetLineStyle(1);
  h_MCttW->SetLineWidth(2);
  // cout << "scale factor:" << scale_MC << endl;
  // hs->Scale(scale_MC);
  //h_MC->Scale(MCXSweight);
  TH1D *h_MC=new TH1D("h_MC","",BinTotal,BinXLow,BinXHig);
  h_MC->Sumw2();
  h_MC->Add(h_MCgg,1.0);
  h_MC->Add(h_MCGJets,1.0);
  h_MC->Add(h_MCWJets,1.0);
  h_MC->Add(h_MCTT,1.0);
  h_MC->Add(h_MCWGJJ,1.0);
  h_MC->Add(h_MCttW,1.0);
  double Ntot_MC=h_MC->Integral();
  if( debug==1 ) cout<<"JTao: N_MC= "<<Ntot_MC<<endl;
  Double_t scale_MC = Ntot_Data*1.0/Ntot_MC;
  cout << "nomalization scale factor = " <<scale_MC<<endl;
  // scale_MC = 1;  
  h_MC->Sumw2();
  h_MC->Scale(scale_MC);  
  h_MCgg->Scale(scale_MC);  
  h_MCGJets->Scale(scale_MC);
  h_MCWJets->Scale(scale_MC);  
  h_MCTT->Scale(scale_MC);  
  h_MCWGJJ->Scale(scale_MC);  
  h_MCttW->Scale(scale_MC); 

  THStack *hs = new THStack("hs","");
  hs->Add(h_MCgg);
  hs->Add(h_MCGJets);
  hs->Add(h_MCWJets);
  hs->Add(h_MCTT);
  hs->Add(h_MCWGJJ);
  hs->Add(h_MCttW);

  double Chi2=0.;
  for(int ibin=0; ibin<BinTotal; ibin++){
    double Nd = h_data->GetBinContent(ibin+1);
    double Nm = h_MC->GetBinContent(ibin+1);
    double NmErr = h_MC->GetBinError(ibin+1);
    
    Chi2 += fabs(NmErr)>1e-9?(Nm-Nd)*(Nm-Nd)*1.0/(NmErr*NmErr):0.0;
  }

  cout<<"JTao: chi2 = "<<Chi2<<endl;

  //Stat Err
  TH1D *htot=new TH1D("htot","",BinTotal,BinXLow,BinXHig);
  TH1D *htot_Norm=new TH1D("htot_Norm","",BinTotal,BinXLow,BinXHig);
  // if(IfStatErr==1){
  //   htot->Add(h_MC,1.0);
  //   htot->Sumw2();
  //   //double Ntot_MC=htot->Integral();
  //   //float scale = Ntot_Data*1.0/Ntot_MC;
  //   //htot->Scale(scale);

  //   }
  // }
  htot->Sumw2();
  for(int ibin=0; ibin<BinTotal; ibin++){
    double Nm = h_MC->GetBinContent(ibin+1);
    double NmErr = h_MC->GetBinError(ibin+1);
    double RelErr=Nm>0?NmErr*1.0/Nm:0.0;
    htot_Norm->SetBinContent(ibin+1, 1);
    //    htot_Norm->SetBinError(ibin+1, RelErr);
    double Nd = h_data->GetBinContent(ibin+1);
    double ScaledRelErr = Nm>0.?Nd/Nm*RelErr:0.0;
    htot_Norm->SetBinError(ibin+1, ScaledRelErr);
  }

  double maxY=max(h_data->GetMaximum(),h_MC->GetMaximum());
  double minY=min(h_data->GetMinimum(),h_MC->GetMinimum());
  //  h_data->GetYaxis()->SetRangeUser(0.95*minY, 1.2*maxY);
  h_data->GetYaxis()->SetRangeUser(0.95*minY, 1.05*maxY);
  //h_data->SetMaximum(1.2*maxY);
  if(IfLogY==1) h_data->GetYaxis()->SetRangeUser(0.1, 1.5*maxY);

  h_data->SetLineColor(1);
  h_data->SetFillColor(0);
  h_data->SetLineStyle(1);
  h_data->SetLineWidth(2);
  //h_data->GetXaxis()->SetTitle(XTitle.c_str());
  double WidthBin=(BinXHig-BinXLow)/BinTotal;
  //TString TitleY( Form("A.U. / %.2g GeV",WidthBin) );
  //TString TitleY( Form("No. of Entries in data / %.2g GeV",WidthBin) );
  //TString TitleY = "A.U";
  string PreTitleY( Form("Events / %.2g ",WidthBin) );
  //  string PreTitleY( Form("No. of Entries / %.2g ",WidthBin) );
  string TitleY =  PreTitleY + YUnit;
  h_data->GetYaxis()->SetTitle(TitleY.c_str());

  h_data->SetTitleSize(0.06,"X");
  h_data->SetTitleSize(0.06,"Y");
  //h_data->SetTitleOffset(1.3, "Y");
  h_data->SetTitleOffset(1.1, "Y");

  h_data->SetMarkerColor(kBlack);
  //h_data->SetMarkerSize(1.0);
  h_data->SetMarkerSize(0.8);
  h_data->SetMarkerStyle(20);

  h_MC->SetFillColor(0);
  h_MC->SetMarkerStyle(0);
  h_MC->SetLineColor(1);
  h_MC->SetLineStyle(1);
  h_MC->SetLineWidth(2);

  //legend->AddEntry(h_data,"data","pe");
  //legend->AddEntry(h_MC,"MC","f");

  legend->AddEntry(h_data,"Data","pe");
  legend->AddEntry(h_MCgg,"#gamma#gamma","f");
  legend->AddEntry(h_MCGJets,"#gammaJets","f");
  legend->AddEntry(h_MCWJets,"WJets","f");
  legend->AddEntry(h_MCTT,"tt","f");
  legend->AddEntry(h_MCWGJJ,"W#gammaJJ","f");
  legend->AddEntry(h_MCttW,"ttW","f");
  // if(IfStatErr==1) 
  legend->AddEntry(htot_Norm, " MC Stat. Err.","f");

  //prepare 2 pads
  const Int_t nx=1;
  const Int_t ny=2;
  const Double_t x1[2] = {0.0,0.0};
  const Double_t x2[2] = {1.0,1.0};
  //const Double_t y1[] = {1.0,0.3};
  //const Double_t y2[] = {0.3,0.00};
  const Double_t y1[2] = {0.3,0.0};
  const Double_t y2[2] = {1.0,0.3};
  Double_t psize[2];
  TPad *pad;
  const char *myname = "c";
  char *name2 = new char [strlen(myname)+6];
  Int_t n = 0;
  for (int iy=0;iy<ny;iy++) {
    for (int ix=0;ix<nx;ix++) {
      n++;
      sprintf(name2,"%s_%d",myname,n);
      if(ix==0){
        gStyle->SetPadLeftMargin(.166);
      }else{
        gStyle->SetPadLeftMargin(.002);
        gStyle->SetPadTopMargin(.002);
      }

      if(iy==0){//upper
        gStyle->SetPadTopMargin(0.05*(1./0.7)); // 0.05
        gStyle->SetPadBottomMargin(.02);
      }
      if(iy==(ny-1)){//lower pad
        gStyle->SetPadTopMargin(.05);
        //gStyle->SetPadBottomMargin(.13*(1./0.3));
        gStyle->SetPadBottomMargin(.40);


      }
      pad = new TPad(name2,name2,x1[ix],y1[iy],x2[ix],y2[iy]);
      pad->SetNumber(n);
      pad->Draw();
      psize[iy]=y1[iy]-y2[iy];
      //if(iy>0 )pad->SetGrid(kTRUE);
    }// end of loop over x
  }// end of loop over y
  delete [] name2;

  //===Drawing====
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);

  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  h_data->GetXaxis()->SetLabelColor(0);
  h_data->SetNdivisions(510 ,"X");

  c1->cd(1);
  gPad->SetLogy(IfLogY);
  //gPad->SetTickx(1);
  //gPad->SetTicky(1);
  //=========
  h_data->Draw("PE1");
  legend->Draw("same");
  hs->Draw("hist,same");
  // h_MCsig->Draw("hist,same");
  //h_MC->Draw("hist,same");
  // if(IfStatErr==1){
  // }
  
  h_data->Draw("samePE1");
  h_data->Draw("Axissame");
  // check different bin content
  cout << "JTao:  1 bin:" << h_MC->GetBinContent(1) << endl;

  /*
    TLatex a;
    a.SetNDC();
    a.SetTextSize(0.05);
    a.DrawLatex(0.2,0.94, PrintInfor);
  */
  const TString PrintInfor1="#bf{CMS} #it{} #it{Preliminary}";
  const TString PrintInfor2="41.5 fb^{-1} (13TeV)"; 
  //tex = new TLatex(0.129,0.93, PrintInfor1);
  //TLatex *tex1 = new TLatex(0.16,0.94, PrintInfor1);
  TLatex *tex1 = new TLatex(0.25,0.94,PrintInfor1);
  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.045);
  tex1->SetLineWidth(2);
  tex1->Draw();

  TLatex *tex2 = new TLatex(0.56,0.933, PrintInfor2);
  tex2 = new TLatex(0.70,0.94, PrintInfor2);
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.045);
  tex2->SetLineWidth(2);
  tex2->Draw();

  ///====
  TLine *Line1 = new TLine(h_data->GetBinLowEdge(1),1,h_data->GetBinLowEdge(h_data->GetNbinsX())+ h_data->GetBinWidth(h_data->GetNbinsX()),1);
  Line1->SetLineColor(1);
  Line1->SetLineWidth(2);
  Line1->SetLineStyle(4);

  TH1D *histoRatio = new TH1D(*h_data);
  histoRatio->Divide(h_data, h_MC, 1., 1.);
    // if(IfStatErr==1){
    //   }
    // }
  for(int ibin=0; ibin<BinTotal; ibin++){
    double Nd = h_data->GetBinContent(ibin+1);
    double NdErr = h_data->GetBinError(ibin+1);
    double Nm = h_MC->GetBinContent(ibin+1);
    histoRatio->SetBinError(ibin+1, Nd/Nm * NdErr/Nd);
    // histoRatio->SetBinError(ibin+1, Nd/Nm * NdErr/Nd);
  }
  for(int ibin=0; ibin<BinTotal; ibin++){
    double Nd = h_MC->GetBinContent(ibin+1);
    double NdErr = h_MC->GetBinError(ibin+1);
    h_MC->SetBinError(ibin+1, NdErr/Nd);
  }
  histoRatio->SetLineColor(1);
  histoRatio->SetLineStyle(1);
  histoRatio->SetLineWidth(2);
  histoRatio->SetMarkerColor(1);
  histoRatio->SetMarkerStyle(20);

  c1->cd(2);
  gPad->SetLogy(0);
  histoRatio->SetTitleOffset(1,"X");
  histoRatio->SetTitleSize(0.12,"X");
  histoRatio->SetLabelSize(0.1,"X");
  histoRatio->GetXaxis()->SetTitle(XTitle.c_str());
  histoRatio->GetYaxis()->SetTitle("Data / MC");
  //  histoRatio->GetYaxis()->SetTitle("data/MC");
  //histoRatio->SetTitleOffset(0.5,"Y");
  histoRatio->SetTitleOffset(0.4,"Y");
  //histoRatio->SetTitleSize(0.12,"Y");
  histoRatio->SetTitleSize(0.14,"Y");
  histoRatio->SetLabelSize(0.1,"Y");
  histoRatio->SetLabelColor(1,"X");
  histoRatio->GetYaxis()->CenterTitle();
  //histoRatio->SetNdivisions(505 ,"Y");
  histoRatio->SetNdivisions(510 ,"X");

  gPad->SetTickx(1);
  gPad->SetTicky(1);
 
  histoRatio->GetXaxis()->SetTickLength(0.08);
  histoRatio->GetYaxis()->SetTickLength(0.06);
  histoRatio->GetYaxis()->SetNdivisions(503);

  //histoRatio->SetMinimum(0.8); //0.5
  //histoRatio->SetMaximum(1.2);  //1.5
  histoRatio->SetMinimum(0.5);
  histoRatio->SetMaximum(1.5);
  histoRatio->Draw();
  // htot->SetFillColorAlpha(kRed,0.35);
  // htot->SetMarkerStyle(1);
  // htot->Draw("same:E2");
  htot_Norm->SetFillColorAlpha(kRed,0.35);
  htot_Norm->SetMarkerStyle(0);
  htot_Norm->Draw("same:E2");
    // if(IfStatErr==1){
    // }
  Line1->Draw("same");

  //===================================
  string nameplots=OutputPlotDir + "/DataMC_"+PlotName+".png";
  c1->Print(nameplots.c_str());

  // string nameplotspdf=OutputPlotDir + "/DataMC_"+PlotName+".pdf";
  // c1->Print(nameplotspdf.c_str());

  //c1->Clear();
  //legend->Clear();
 
}


void DrawDataMCPlots_new(){

  gROOT->ProcessLine(".x hggPaperStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);	

  
  //All
  // string Preselections="((CMS_hgg_mass <= 115. || CMS_hgg_mass >= 135.))";
  // string Preselections="(CMS_hgg_mass >= -100000.)";

  // string PreselectionsMVA0=Preselections + " && diphoMVA>0.";
  // string PreselectionsMVA05=Preselections + " && diphoMVA>0.5";
  // string PreselectionsMVA08=Preselections + " && diphoMVA>0.8";
  // string PreselectionsMVA95=Preselections + " && diphoMVA>0.95";
  // string PreselectionsMVA99=Preselections + " && diphoMVA>0.99";
  // DrawMyPlots("DiphotonMVA_self", Preselections, "New Diphoton BDT score", "GeV", "NewDiphotonBDT", 40, -1., 1., 1, 0);
  // return;
  DrawMyPlots("CMS_hgg_mass", Preselections, "m_{#gamma#gamma} (GeV)", "GeV", "DiphotonMass", 50, 100., 180., 0, 0);
  // DrawMyPlots("sclaed_subleadphoton_pt", Preselections, "sclaed_subleadphoton_pt (GeV)", "GeV", "sclaed_subleadphoton_pt", 50, 0., 2., 0, 0);

  // DrawMyPlots("diphoMVA", Preselections, "Diphoton BDT score", "GeV", "TransformedDiphotonBDT", 40, -1., 1., 1, 0);
  // return;
  // DrawMyPlots("1. / (1. + exp(0.5 * log(2./(DiphotonMVA_self + 1.)-1.)))", Preselections, "Transformed new Diphoton BDT score", "GeV", "Transformed_New_Diphoton_BDT_score", 50, 0., 1., 1, 0);
  // return;
  // //VBF
  // // DrawMyPlots("dijet_prob_VBF_value", Preselections, "dijet_prob_VBF_value ", "", "dijet_prob_VBF_value", 50, 0., 1., 0, 0);
  // // DrawMyPlots("dijet_prob_VBF_value+diphoMVA", Preselections, "dijet_prob_VBF_value+Diphoton BDT Score", "", "dijet_prob_VBF_valuePlusDiphotonBDTScore", 60, -1.0, 2.0, 1, 0);
  // //return;

  // DrawMyPlots("CMS_hgg_mass", PreselectionsMVA0, "m_{#gamma#gamma} (GeV)", "GeV", "DiphotonMass_MVAG0", 80, 100., 180., 0, 0);
  // DrawMyPlots("CMS_hgg_mass", PreselectionsMVA05, "m_{#gamma#gamma} (GeV)", "GeV", "DiphotonMass_MVAG05", 80, 100., 180., 0, 0);
  // DrawMyPlots("CMS_hgg_mass", PreselectionsMVA08, "m_{#gamma#gamma} (GeV)", "GeV", "DiphotonMass_MVAG08", 80, 100., 180., 0, 0);
  // DrawMyPlots("CMS_hgg_mass", PreselectionsMVA95, "m_{#gamma#gamma} (GeV)", "GeV", "DiphotonMass_MVAG95", 80, 100., 180., 0, 0);
  // DrawMyPlots("CMS_hgg_mass", PreselectionsMVA99, "m_{#gamma#gamma} (GeV)", "GeV", "DiphotonMass_MVAG99", 80, 100., 180., 0, 0);

  // //------------------------
    // DrawMyPlots("leadR9", Preselections, "Leading #gamma R_{9}", "", "leadR9", 40, 0.4, 1.2, 1, 0);
    // DrawMyPlots("subleadR9", Preselections, "Subleading #gamma R_{9}", "", "subleadR9", 40, 0.4, 1.2, 1, 0);

    // DrawMyPlots("leadptom", Preselections, "p^{leading #gamma}_{T}/m_{#gamma#gamma}", "", "leadptom", 100, 0., 2., 0, 0);
    // DrawMyPlots("subleadptom", Preselections, "p^{subleading #gamma}_{T}/m_{#gamma#gamma}", "", "subleadptom", 100, 0., 1., 0, 0);
    // DrawMyPlots("leadmva", Preselections, "Leading #gamma ID scores", "", "leadmva", 40, -1.0, 1.0, 1, 0);
    // DrawMyPlots("subleadmva", Preselections, "Subleading #gamma ID scores", "", "subleadmva", 40, -1.0, 1.0, 1, 0);
    // DrawMyPlots("min(leadmva,subleadmva)", Preselections, "Min. #gamma ID scores", "", "minmva", 40, -1.0, 1.0, 1, 0);
    // DrawMyPlots("max(leadmva,subleadmva)", Preselections, "Max. #gamma ID scores", "", "maxmva", 40, -1.0, 1.0, 1, 0);

  // DrawMyPlots("leadeta", Preselections, "Leading #gamma #eta_{SC}", "", "leadeta", 50, -2.5, 2.5, 0, 0);
  // DrawMyPlots("subleadeta", Preselections, "Subleading #gamma #eta_{SC}", "", "subleadeta", 50, -2.5, 2.5, 0, 0);
  // DrawMyPlots("sigmarv", Preselections, "#sigma_{rv}", "", "sigmarv", 50, 0., 0.05, 0, 0);
  // DrawMyPlots("sigmawv", Preselections, "#sigma_{wv}", "", "sigmawv", 50, 0., 0.05, 0, 0);
  // DrawMyPlots("sigmarvDecorr", Preselections, "#sigma_{rv}", "", "sigmarvDecorr", 50, 0., 0.05, 0, 0);

  // DrawMyPlots("CosPhi", Preselections, "cos#phi_{#gamma#gamma}", "", "CosPhi", 50, -1., 1., 0, 0);
  // DrawMyPlots("vtxprob", Preselections, "Vertex probability", "", "vtxprob", 50, 0., 1., 0, 0);
  return;

  // DrawMyPlots("pt", Preselections, "p_{T}^{#gamma#gamma} (GeV)", "GeV", "pt", 50, 0., 500., 0, 0);
  // DrawMyPlots("pt", Preselections, "p_{T}^{#gamma#gamma} (GeV)", "GeV", "pt0to50", 50, 0., 50., 1, 0);

  //   DrawMyPlots("rho", Preselections, "#rho", "", "rho", 60, 0., 60., 0, 0);
  // DrawMyPlots("nvtx", Preselections, "#vertices", "", "nvtx", 80, 0., 80., 0, 0);

  // //  DrawMyPlots("leadmva<subleadmva?leadSigEOverE:subleadSigEOverE", Preselections, "Min. ID MVA #gamma #sigma_{E}/E", "", "minidSigEOverE", 40, 0., 0.08, 0, 0);
  /*
  DrawMyPlots("minIDSigEOverE", Preselections, "Min. ID MVA #gamma #sigma_{E}/E", "", "minidSigEOverE", 40, 0., 0.08, 0, 0);
  DrawMyPlots("maxIDSigEOverE", Preselections, "Max. ID MVA #gamma #sigma_{E}/E", "", "maxidSigEOverE", 40, 0., 0.08, 0, 0);
 
  DrawMyPlots("leadpt", Preselections, "Leading #gamma p_{T} (GeV)", "GeV", "leadpt", 60, 30., 150., 0, 0);
  DrawMyPlots("subleadpt", Preselections, "Subleading #gamma p_{T} (GeV)", "GeV", "subleadpt", 60, 20., 140., 0, 0);

  DrawMyPlots("minIDpt", Preselections, "Min. ID MVA #gamma p_{T} (GeV)", "GeV", "minIDpt", 60, 20., 140., 0, 0);
  DrawMyPlots("maxIDpt", Preselections, "Max. ID MVA #gamma p_{T} (GeV)", "GeV", "maxIDpt", 60, 20., 140., 0, 0);
  */

} 