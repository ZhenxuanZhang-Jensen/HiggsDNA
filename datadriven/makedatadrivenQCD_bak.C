
void makedatadrivenQCD()
{

  // Inputs
  string common_selection = "(Leading_Photon_pt/CMS_hgg_mass>0.35)*(Subleading_Photon_pt/CMS_hgg_mass>0.25)*(passPhotonSels==1)*(CMS_hgg_mass<115 || CMS_hgg_mass>135)";
  
  vector<string> gjetfilenames = {
    //"/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/DNN_MoreVar/Backgrounds/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV.root",
    //"/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/DNN_MoreVar/Backgrounds/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf.root"
    "/eos/user/z/zhenxuan/hhwwgg_root/GJet_Pt-20to40_2017.root",
    "/eos/user/z/zhenxuan/hhwwgg_root/GJet_Pt-40toInf_2017.root"
  };

  vector<string> gjettreenames = {
    //"tagsDumper/trees/GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf_TuneCP5_13TeV_Pythia8_13TeV_HHWWggTag_1",
    //"tagsDumper/trees/GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf_TuneCP5_13TeV_Pythia8_13TeV_HHWWggTag_1"
    "GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf_TuneCP5_13TeV_Pythia8_13TeV_HHWWggTag_1",
    "GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf_TuneCP5_13TeV_Pythia8_13TeV_HHWWggTag_1"
  };

  vector<float> gjetlumis = {
    41.5,
    41.5
  };

  vector<string> datafilenames = {
    //"/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/DNN_MoreVar/Backgrounds/Data_2017.root"
    "/eos/user/a/atishelm/ntuples/HHWWgg_flashgg/WWZ_SignalTopology_Checks/2017/SingleElectron_Data_2017_hadded/SingleElectron_Data_2017.root"
  };

  vector<string> datatreenames = {
    //"tagsDumper/trees/Data_13TeV_HHWWggTag_1"
    "Data_13TeV_HHWWggTag_1"
  };

  
  //Get photon ID distribution for the fake photons of the g+jets samples
  map<string,TH1F*> h_minphotonID;
  h_minphotonID["gjet"] = new TH1F("h_minphotonID_gjet","h_minphotonID_gjet",19,-0.9,1);
  for( unsigned ifile=0; ifile<gjetfilenames.size(); ++ifile){
    TChain* ch_gjet = new TChain();
    ch_gjet->Add( Form("%s/%s",gjetfilenames[ifile].c_str(),gjettreenames[ifile].c_str()) );
    ch_gjet->Draw( "Leading_Photon_MVA >>+ h_minphotonID_gjet",
		   Form("(Leading_Photon_MVA<Subleading_Photon_MVA)*(Leading_Photon_genMatchType!=1)*(weight)*(%f)*(%s)",gjetlumis[ifile],common_selection.c_str()),
		   "goff");
    ch_gjet->Draw( "Subleading_Photon_MVA >>+ h_minphotonID_gjet",
		   Form("(Subleading_Photon_MVA<Leading_Photon_MVA)*(Subleading_Photon_genMatchType!=1)*(weight)*(%f)*(%s)",gjetlumis[ifile],common_selection.c_str()),
		   "goff");
    cout<<h_minphotonID["gjet"]->GetEntries();

  }

  //Fit the photon ID distribution for the fake photons with a pol7
  h_minphotonID["gjet"]->Scale(1./h_minphotonID["gjet"]->Integral());
  TF1* photonIDPDF_fake = new TF1("photonIDPDF_fake","pol7",-0.9,1.);
  h_minphotonID["gjet"]->Fit(photonIDPDF_fake,"R");
  TCanvas* c1 = new TCanvas();
  h_minphotonID["gjet"]->Draw("E1");
  c1->Print("gjetminphotonID.pdf");

  // Take the data events with photonID_min<-0.7 and photonID_max>-0.7, 
  // re-generate randomly the photonID_min distribution in [-0.7,photonID_max], 
  // reweight the events to match the photonID_max distribution 
  TChain* ch_data = new TChain();
  float Leading_Photon_MVA, Subleading_Photon_MVA, PhotonID_min, original_Photon_MVA_min, weight;
  for( unsigned ifile=0; ifile<datafilenames.size(); ++ifile)
    ch_data->Add( Form("%s/%s",datafilenames[ifile].c_str(),datatreenames[ifile].c_str()) );
  TFile* outfile = new TFile("datadrivenQCD_v2.root","RECREATE");
  // data-(...)=a*GJet+b*QCD 

  outfile->cd();
  outfile->mkdir("tagsDumper/trees");
  outfile->cd("tagsDumper/trees");
  ch_data->SetBranchStatus("*",1);
  ch_data->SetBranchAddress("PhotonID_min",&PhotonID_min);
  ch_data->SetBranchAddress("Leading_Photon_MVA",&Leading_Photon_MVA);
  ch_data->SetBranchAddress("Subleading_Photon_MVA",&Subleading_Photon_MVA);
  ch_data->SetBranchAddress("weight",&weight);

  TTree* outtree = ch_data->CloneTree(0);
  outtree->Branch("original_Photon_MVA_min", original_Photon_MVA_min);
  
  int nentries = ch_data->GetEntries();
  for(int ientry=0; ientry<nentries; ++ientry){
    ch_data->GetEntry(ientry);
    if(ientry%1000==0)
      cout<<"reading entry "<<ientry<<endl;
    
    bool hasleadIDmin;
    double Photon_MVA_max;
    if(Leading_Photon_MVA < Subleading_Photon_MVA) {
      hasleadIDmin=true;
      original_Photon_MVA_min = Leading_Photon_MVA;
      Photon_MVA_max = Subleading_Photon_MVA;
    }
    else {
      hasleadIDmin=false;
      original_Photon_MVA_min = Subleading_Photon_MVA;
      Photon_MVA_max = Leading_Photon_MVA;
    }      
    
    if( !(original_Photon_MVA_min<-0.7 && Photon_MVA_max>-0.7) )
      continue;

    if(hasleadIDmin){
      Leading_Photon_MVA = photonIDPDF_fake->GetRandom( -0.7, Photon_MVA_max );
      PhotonID_min = Leading_Photon_MVA;
    }
    else{
      Subleading_Photon_MVA = photonIDPDF_fake->GetRandom(-0.7, Photon_MVA_max );
      PhotonID_min = Subleading_Photon_MVA;
    }

    weight *= photonIDPDF_fake->Integral(-0.7,Photon_MVA_max) / photonIDPDF_fake->Integral(-0.9,-0.7);
    outtree->Fill();
  }
  outtree->AutoSave();
  delete ch_data;
  outfile->Close(); 
  
};
  