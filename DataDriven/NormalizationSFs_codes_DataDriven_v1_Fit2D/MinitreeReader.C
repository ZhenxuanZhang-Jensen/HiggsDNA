#define MinitreeReader_cxx
#include "MinitreeReader.h"
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using std::min;
using std::max;

void MinitreeReader::Loop(TString outputfile)
{
//   In a ROOT session, you can do:
//      root> .L MinitreeReader.C
//      root> MinitreeReader t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
   fChain->SetBranchStatus("*",0);  // disable all branches
   if (!fIsData) {
      fChain->SetBranchStatus("weight",1);  // activate weight
   } else {
      weight = 1;
   }
   fChain->SetBranchStatus("LeadPhoton_mvaID",1);  // activate LeadPhoton_mvaID
   fChain->SetBranchStatus("SubleadPhoton_mvaID",1);  // activate SubleadPhoton_mvaID
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile *fout = TFile::Open(outputfile, "RECREATE");
   TDirectory *dir_pp = fout->mkdir("pp");
   dir_pp->cd();
   TH1D *h_minmva_pp = new TH1D("h_minmva_pp", ";photon ID MVA (min);Entries", 19, -0.9, 1);
   TH1D *h_maxmva_pp = new TH1D("h_maxmva_pp", ";photon ID MVA (max);Entries", 20, -1, 1);
   //attention: change(-1,1) to (-0.9,1)
   TH2D *h_minmaxmva_pp = new TH2D("h_minmaxmva_pp", ";photon ID MVA (min);photon ID MVA (max)", 20, -0.9, 1, 20, -0.9, 1);
   TDirectory *dir_pf = fout->mkdir("pf");
   dir_pf->cd();
   TH1D *h_LeadPhoton_mvaID_pf = new TH1D("h_LeadPhoton_mvaID_pf", ";photon ID MVA (leading);Entries", 20, -1, 1);
   TH1D *h_SubleadPhoton_mvaID_pf = new TH1D("h_SubleadPhoton_mvaID_pf", ";photon ID MVA (subleading);Entries", 20, -1, 1);
   TH1D *h_minmva_pf = new TH1D("h_minmva_pf", ";photon ID MVA (min);Entries", 20, -0.9, 1);
   TH1D *h_maxmva_pf = new TH1D("h_maxmva_pf", ";photon ID MVA (max);Entries", 20, -1, 1);
   //attention: change(-1,1) to (-0.9,1)
   TH2D *h_minmaxmva_pf = new TH2D("h_minmaxmva_pf", ";photon ID MVA (min);photon ID MVA (max)", 19, -0.9, 1, 19, -0.9, 1);
   TDirectory *dir_ff = fout->mkdir("ff");
   dir_ff->cd();
   TH1D *h_LeadPhoton_mvaID_ff = new TH1D("h_LeadPhoton_mvaID_ff", ";photon ID MVA (leading);Entries", 20, -1, 1);
   TH1D *h_SubleadPhoton_mvaID_ff = new TH1D("h_SubleadPhoton_mvaID_ff", ";photon ID MVA (subleading);Entries", 20, -1, 1);
   TH1D *h_minmva_ff = new TH1D("h_minmva_ff", ";photon ID MVA (min);Entries", 20, -0.9, 1);
   TH1D *h_maxmva_ff = new TH1D("h_maxmva_ff", ";photon ID MVA (max);Entries", 20, -1, 1);
   //attention: change(-1,1) to (-0.9,1)
   TH2D *h_minmaxmva_ff = new TH2D("h_minmaxmva_ff", ";photon ID MVA (min);photon ID MVA (max)", 20, -0.9, 1, 20, -0.9, 1);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      // mass selection
      b_CMS_hgg_mass->GetEntry(ientry); //read only this branch
      if (CMS_hgg_mass>115 && CMS_hgg_mass<135) continue;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // pp: no matching
      h_LeadPhoton_mvaID_pp->Fill(LeadPhoton_mvaID,weight); 
      h_subLeadPhoton_mvaID_pp->Fill(SubleadPhoton_mvaID,weight); 
      h_minmva_pp->Fill(min(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
      h_maxmva_pp->Fill(max(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
      h_minmaxmva_pp->Fill(min(LeadPhoton_mvaID,SubleadPhoton_mvaID),max(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 

      // pf
      if (!fIsData
               ) {
         h_LeadPhoton_mvaID_pf->Fill(LeadPhoton_mvaID,weight); 
         h_SubleadPhoton_mvaID_pf->Fill(SubleadPhoton_mvaID,weight); 
         h_minmva_pf->Fill(min(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
         h_maxmva_pf->Fill(max(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
         h_minmaxmva_pf->Fill(min(LeadPhoton_mvaID,SubleadPhoton_mvaID),max(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
      } else if (!fIsData) {
         h_LeadPhoton_mvaID_ff->Fill(LeadPhoton_mvaID,weight); 
         h_SubleadPhoton_mvaID_ff->Fill(SubleadPhoton_mvaID,weight); 
         h_minmva_ff->Fill(min(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
         h_maxmva_ff->Fill(max(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
         h_minmaxmva_ff->Fill(min(LeadPhoton_mvaID,SubleadPhoton_mvaID),max(LeadPhoton_mvaID,SubleadPhoton_mvaID),weight); 
      }
   }

   fout->Write();
   fout->Close();
}
