//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep  3 17:54:37 2021 by ROOT version 6.22/08
// from TTree MCBkg_13TeV_UntaggedTag/MCBkg_13TeV_UntaggedTag
// found on file: /eos/cms/store/user/jtao/MassInterference/CMSSW_10_6_8/MiniTree/output_DiPhotonJetsBox.root
//////////////////////////////////////////////////////////

#ifndef MinitreeReader_h
#define MinitreeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.

class MinitreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Bool_t          fIsData;  // is it data?

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         weight;
   Float_t         LeadPhoton_mvaID;
   Float_t         SubleadPhoton_mvaID;

   // List of branches
   TBranch        *b_LeadPhoton_mvaID;   //!
   TBranch        *b_SubleadPhoton_mvaID;   //!

   MinitreeReader(TTree *tree=0, Bool_t isData=false);
   virtual ~MinitreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString outputfile);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MinitreeReader_cxx
MinitreeReader::MinitreeReader(TTree *tree, Bool_t isData) : fChain(0), fIsData(isData)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/user/jtao/MassInterference/CMSSW_10_6_8/MiniTree/output_DiPhotonJetsBox.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/user/jtao/MassInterference/CMSSW_10_6_8/MiniTree/output_DiPhotonJetsBox.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/cms/store/user/jtao/MassInterference/CMSSW_10_6_8/MiniTree/output_DiPhotonJetsBox.root:/tagsDumper/trees");
      dir->GetObject("MCBkg_13TeV_UntaggedTag",tree);

   }
   Init(tree);
}

MinitreeReader::~MinitreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MinitreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MinitreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MinitreeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

Bool_t MinitreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MinitreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MinitreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MinitreeReader_cxx
