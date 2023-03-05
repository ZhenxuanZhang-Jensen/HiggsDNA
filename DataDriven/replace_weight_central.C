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
#include "TMVA/Reader.h"

using namespace std;

void func_replace_weight_central(TFile* input_file, const char* tree_name)
{
    // Open the tree in the root file
    TTree* tree = dynamic_cast<TTree*>(input_file->Get(tree_name));
    if (!tree) {
        std::cerr << "Failed to retrieve tree from input file." << std::endl;
        return;
    }

    // Get the weight and weight_central branches
    float weight, weight_central;
    TBranch* b_weight = tree->GetBranch("weight");
    TBranch* b_weight_central = tree->GetBranch("weight_central");
    if (!b_weight || !b_weight_central) {
        std::cerr << "Failed to retrieve weight or weight_central branch from tree." << std::endl;
        return;
    }

    // Set branch addresses
    b_weight->SetAddress(&weight);
    b_weight_central->SetAddress(&weight_central);

    // Loop over entries and replace weight_central with weight
    for (int i = 0; i < tree->GetEntries(); i++) {
        b_weight->GetEntry(i);
        b_weight_central->GetEntry(i);
        weight_central = weight;
        b_weight_central->Fill();
    }

    // Write the modified tree to the input file
    input_file->cd();
    tree->Write("", TObject::kOverwrite);

    std::cout << "Successfully replaced weight_central with weight." << std::endl;
}
void replace_weight_central(){
    const string InputTree = "cat1";
    const string OutputFileName = "/eos/user/z/zhenxuan/hhwwgg_root/DataDriven_cat1.root";

    func_replace_weight_central(OutputFileName.c_str(), InputTree.c_str());
}