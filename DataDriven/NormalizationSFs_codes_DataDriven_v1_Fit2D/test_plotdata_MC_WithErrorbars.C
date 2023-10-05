#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TFractionFitter.h"

#include <iostream>

using namespace std;

const double lumi_ = 41.5; // fb-1
void test_plotdata_MC_WithErrorbars(){
    TFile *fdata = TFile::Open("histos_dataUL2017.root");
    TCanvas *c1 = new TCanvas("c1", "",0,0, 800, 800);
    gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   //c1->SetGridx();
   //c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.02);
   c1->SetTopMargin(0.08);
   c1->SetBottomMargin(0.1);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   TString varname = "diphoMVA";
   TH1D *hdata = (TH1D*) fdata->Get("pp/h_" + varname + "_pp");
   TH1D *h2data = new TH1D("h2data","hist with data with error bars",20,-0.9,1);
   for (Int_t i =1 ; i<=20;i++)
   {
        Int_t bin_content = hdata->GetBinContent(i);
        h2data->SetBinContent(i,i*10000);
   }
    h2data->Draw("ep");
   c1->SaveAs("test_DataWithErrorBars.png");
}