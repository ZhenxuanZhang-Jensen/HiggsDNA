/////////////////////////////////////////////////////////////////////////////////
// This is similar to plotHistos.C, except that we are using only 2 components //
// (pp, pf+ff) instead of 3, where pp+pf is data-driven.                       //
/////////////////////////////////////////////////////////////////////////////////

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

// scale factors (data-driven pf+ff, dipho)
// double sf[2] = {0.91876, 1.72288}; // minmaxmva 4x4
// double sf[2] = {0.919982, 1.71277}; // minmaxmva 10x10
double sf[2] = {0.919982, 1.71277}; // minmaxmva 10x10

void plotHistos_DataDriven_1var(TFile *fdata, TFile *fdipho, TFile *fddqcd, TString varname, bool postfit);
void fitSFs(TString varname="minmaxmva", TString datafile="histos_dataUL2017.root", TString diphofile="histos_DiPhotonJetsBox.root", TString ddqcdfile="histos_DataDrivenQCD.root");

void plotHistos_DataDriven(bool postfit=true, TString datafile="histos_dataUL2017.root", TString diphofile="histos_DiPhotonJetsBox.root", TString ddqcdfile="histos_DataDrivenQCD.root") {
   TFile *fdata = TFile::Open(datafile);
   TFile *fdipho = TFile::Open(diphofile);
   TFile *fddqcd = TFile::Open(ddqcdfile);

   TString var[18] = {
      "mass", "pt0800", "pt0200", "pt040", "diphoMVA",
      "leadmva", "subleadmva", "leadptom", "subleadptom", "minmva", 
      "maxmva", "sigmarvDecorr", "sigmarv", "sigmawv", "leadpt",
      "subleadpt", "leadeta", "subleadeta"
   };
   for (int i=0; i<18; i++) plotHistos_DataDriven_1var(fdata, fdipho, fddqcd, var[i], postfit);
}

void plotHistos_DataDriven_1var(TFile *fdata, TFile *fdipho, TFile *fddqcd, TString varname, bool postfit) {
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

   // ------------>Primitives in pad: TopPad
   TPad *TopPad = new TPad("TopPad", "TopPad",0.01,0.28,0.99,0.99);
   TopPad->Draw();
   TopPad->cd();
   TopPad->SetFillStyle(0);
   TopPad->SetBorderMode(0);
   TopPad->SetBorderSize(2);
   TopPad->SetTickx(1);
   TopPad->SetTicky(1);
   TopPad->SetLeftMargin(0.142);
   TopPad->SetRightMargin(0.04);
   TopPad->SetTopMargin(0.08);
   TopPad->SetBottomMargin(0.02/(1.-0.26));
   TopPad->SetFrameFillStyle(0);
   TopPad->SetFrameBorderMode(0);
   TopPad->SetFrameFillStyle(0);
   TopPad->SetFrameBorderMode(0);
   if (varname.Contains("pt")) TopPad->SetLogy();

   TH1D *hdata = (TH1D*) fdata->Get("pp/h_" + varname + "_pp");
   TH1D *hdipho = (TH1D*) fdipho->Get("pp/h_" + varname + "_pp");
   TH1D *hddqcd = (TH1D*) fddqcd->Get("pp/h_" + varname + "_pp");

   // histos are normalised to 1fb-1 by default
   hdipho->Scale(lumi_);
   // hddqcd->Scale(lumi_);

   // also scale by scale factor
   if (postfit) {
      hdipho->Scale(sf[1]);
      hddqcd->Scale(sf[0]);
   }

   // stack MC
   THStack* hstack = new THStack("hstack", "hstack");
   TH1D* htotal = (TH1D*) hdipho->Clone("htotal");
   htotal->Reset();
   TH1D* hratio = (TH1D*) hdata->Clone("hratio");

   // set style
   hdipho->SetLineColor(kYellow+1);
   hdipho->SetMarkerColor(kYellow+1);
   hdipho->SetFillColor(kYellow+1);
   hdipho->SetFillStyle(1001);
   hddqcd->SetLineColor(kRed);
   hddqcd->SetMarkerColor(kRed);
   hddqcd->SetFillColor(kRed);
   hddqcd->SetFillStyle(1001);

   hstack->Add(hdipho, "hist");
   hstack->Add(hddqcd, "hist");

   htotal->Add(hddqcd);
   htotal->Add(hdipho);

   cout << "data / MC: " << hdata->Integral() << " / " << htotal->Integral() << " = " << hdata->Integral() / htotal->Integral() << endl;

   // draw
   hdata->GetXaxis()->SetTitleSize(0.);
   hdata->GetXaxis()->SetLabelSize(0.);
   hdata->GetYaxis()->SetTitleSize(0.05);
   hdata->GetYaxis()->SetLabelSize(0.05);
   hdata->SetMarkerStyle(20);
   hdata->SetMarkerSize(1.5);
   hdata->Draw("E1P");
   hstack->Draw("histsame");
   hdata->Draw("E1Psame");
   TopPad->RedrawAxis();

   // legend
   TLegend *leg = new TLegend(0.6,0.6,0.92,0.90,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hdata, "Data", "lp");
   leg->AddEntry(hdipho, "#gamma#gamma", "f");
   leg->AddEntry(hddqcd, "Data-driven QCD", "f");
   leg->Draw();
   
   // goodness-of-fit
   double chi2 = hdata->Chi2Test(htotal, "UW CHI2");
   double ndf = chi2 / hdata->Chi2Test(htotal, "UW CHI2/NDF");
   double chi2prob = TMath::Prob(chi2,ndf);
   double ksprob = hdata->KolmogorovTest(htotal);
   TLatex *t = new TLatex(.96,.935,Form("#chi^{2}/ndf = %.1f/%.0f (prob=%.2f), KS prob=%.2f",chi2,ndf,chi2prob,ksprob));
   t->SetNDC();
   t->SetTextAlign(31);
   t->SetTextSize(0.04);
   t->Draw();

   // ------------>Primitives in pad: bottomPad
   c1->cd();
   TPad *bottomPad = new TPad("bottomPad", "bottomPad",0.01,0.01,0.99,0.3);
   bottomPad->Draw();
   bottomPad->cd();
   bottomPad->SetFillStyle(0);
   bottomPad->SetBorderMode(0);
   bottomPad->SetBorderSize(2);
   bottomPad->SetTickx(1);
   bottomPad->SetTicky(1);
   bottomPad->SetLeftMargin(0.142);
   bottomPad->SetRightMargin(0.04);
   bottomPad->SetTopMargin(0.);
   bottomPad->SetBottomMargin(0.35);
   bottomPad->SetFrameFillStyle(0);
   bottomPad->SetFrameBorderMode(0);

   hratio->Divide(htotal);
   hratio->GetXaxis()->SetLabelSize(0.11);
   hratio->GetXaxis()->SetTitleSize(0.11);
   hratio->GetYaxis()->SetLabelSize(0.11);
   hratio->GetYaxis()->SetTitleSize(0.11);
   hratio->GetYaxis()->SetTitleOffset(0.39);
   hratio->GetYaxis()->SetTitle("data/MC");
   hratio->SetMinimum(0.71);
   hratio->SetMaximum(1.29);
   hratio->SetFillColor(kBlack);
   hratio->SetFillStyle(1001);
   hratio->SetMarkerStyle(20);
   hratio->SetMarkerSize(1.5);
   hratio->Draw();
   TLine gridRatio;
   gridRatio.SetLineColor(kRed);
   gridRatio.SetLineStyle(2);
   double linemin = hratio->GetBinCenter(1)-hratio->GetBinWidth(1)/2.;
   double linemax = hratio->GetBinCenter(hratio->GetNbinsX())+hratio->GetBinWidth(hratio->GetNbinsX())/2.;
   gridRatio.DrawLine(linemin,1.0,linemax,1.0);

   c1->SaveAs(varname + (postfit ? "_postfit" : "_prefit") + ".pdf");
   c1->SaveAs(varname + (postfit ? "_postfit" : "_prefit") + ".png");
}

void fitSFs(TString varname, TString datafile, TString diphofile, TString ddqcdfile) {
   TFile *fdata = TFile::Open(datafile);
   TFile *fdipho = TFile::Open(diphofile);
   TFile *fddqcd = TFile::Open(ddqcdfile);

   TH1 *hdata = (TH1*) fdata->Get("pp/h_" + varname + "_pp");
   TH1 *hdipho = (TH1*) fdipho->Get("pp/h_" + varname + "_pp");
   TH1 *hddqcd = (TH1*) fddqcd->Get("pp/h_" + varname + "_pp");

   hdipho->Scale(lumi_);
   // hddqcd->Scale(lumi_);

   double init[2], final[2], error[2];
   TH1 *hmc[2] {hddqcd, hdipho};
   TString label[2] = {
      "QCD (Data-driven)",
      "gammagamma + jets (prompt/prompt)"
   };

   cout << "Initial fractions:" << endl;
   for (int i=0; i<2; i++) {
      init[i] = hmc[i]->Integral() / hdata->Integral();
      cout << label[i] << ": " << init[i] << endl;
   }

   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   for (int i=0; i<2; i++) mc->Add(hmc[i]);
   TFractionFitter* fit = new TFractionFitter(hdata, mc); // initialise
   for (int i=0; i<2; i++) fit->Constrain(i,0.0,1.0);               // constrain fractions to be between 0 and 1
   // fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
   Int_t status = fit->Fit();               // perform the fit
   cout << "fit status: " << status << endl;
   cout << "chi2/ndf = " << fit->GetChisquare() << "/" << fit->GetNDF() << " (prob. = " << TMath::Prob(fit->GetChisquare(),fit->GetNDF()) << endl;
   if (status == 0) {                       // check on fit status
      TCanvas c1(varname+"_fractionfit");
      TH1F* result = (TH1F*) fit->GetPlot();
      hdata->Draw("Ep");
      result->Draw("same");
      c1.SaveAs(varname+"_fractionfit.pdf");

      cout << "Fitted fractions:" << endl;
      for (int i=0; i<2; i++) {
         fit->GetResult(i, final[i], error[i]);
         cout << label[i] << ": " << final[i] << " +/- " << error[i] << endl;
      }

      cout << "Scale factors:" << endl;
      for (int i=0; i<2; i++) {
         cout << label[i] << ": " << final[i] / init[i] << endl;
      }
      cout << "double sf[2] = {" << final[0] / init[0] << ", " << final[1] / init[1] << "}; // " << varname << endl;
   }
}
