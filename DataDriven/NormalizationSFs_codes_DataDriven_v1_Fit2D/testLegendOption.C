#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"


void testLegendOption()
{
  gStyle->SetOptStat(0);
  

  Double_t mhmin = 100;
  Double_t mhmax = 109;

  TF1* f1 = new TF1("f1", "gaus(0)", mhmin, mhmax);

  f1->SetParameters(1, 0.5*(mhmax+mhmin), 0.2*(mhmax-mhmin));

  TH1F* h1 = new TH1F("h1", "h1", 10*Int_t(mhmax-mhmin), mhmin, mhmax);

  h1->FillRandom("f1", 10000);


  // Draw
  //----------------------------------------------------------------------------
  TCanvas* c1 = new TCanvas();

  c1->cd();
  
  h1->SetMarkerStyle(kFullCircle);
  h1->SetTitle("");

  h1->Draw("ep");

  TLegend* legend = new TLegend(0.68, 0.72, 0.88, 0.89);

  legend->AddEntry("h1", " data", "ep");

  legend->SetBorderSize( 0);
  legend->SetFillColor ( 0);
  legend->SetTextFont  (42);

  legend->Draw("same");
  c1->SaveAs("test.png");
}