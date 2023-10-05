// g++ -o read_minitree `root-config --cflags --ldflags --libs` main_MinitreeReader.C
#include "MinitreeReader.C"

#include <iostream>

using namespace std;

int main(int argc, const char** argv) {
   TString inputfile(argv[1]);
   TString outputfile(argv[2]);
   bool isData = inputfile.Contains("Data") && ! inputfile.Contains("DataDriven");
   if (isData) cout << "We believe that " << inputfile.Data() << " is data." << endl;
   TFile *fin = TFile::Open(inputfile);
   TTree *tr;
   if (isData) tr = (TTree*) fin->Get("cat1");
   else {
      if (inputfile.Contains("DataDriven")) tr = (TTree*) fin->Get("cat1");
      else tr = (TTree*) fin->Get("cat1");
   }
   MinitreeReader toto(tr, isData);
   toto.Loop(outputfile);
}
