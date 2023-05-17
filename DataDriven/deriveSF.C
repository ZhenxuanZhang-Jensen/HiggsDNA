
class chisquare 
{

public:
  chisquare(TH1F* data1, TH1F* constMC1, TH1F* diphotonMC1, TH1F* datadrivenQCD1,
	    TH1F* data2, TH1F* constMC2, TH1F* diphotonMC2, TH1F* datadrivenQCD2)
  {
    data1_ = data1; 
    constMC1_ = constMC1; 
    diphotonMC1_ = diphotonMC1; 
    datadrivenQCD1_ = datadrivenQCD1;
    data2_ = data2; 
    constMC2_ = constMC2; 
    diphotonMC2_ = diphotonMC2; 
    datadrivenQCD2_ = datadrivenQCD2;
    SFdiphoton_ = 1.;
    SFdatadrivenQCD_ = 1.;
  }
  double operator()( double* SFs, double *p)
  {
    SFdiphoton_ = SFs[0];
    SFdatadrivenQCD_ = SFs[1];
    double chisquarevalue = 0.;
    for(int ibin=1; ibin<=data1_->GetXaxis()->GetNbins(); ++ibin) {
      double Nev_data1 = data1_->GetBinContent(ibin);
      double Nev_constMC1 = constMC1_->GetBinContent(ibin);
      double Nev_diphotonMC1 = SFdiphoton_ * diphotonMC1_->GetBinContent(ibin);
      double Nev_datadrivenQCDMC1 = SFdatadrivenQCD_ * datadrivenQCD1_->GetBinContent(ibin);
      double Nev_exp1 = Nev_constMC1+Nev_diphotonMC1+Nev_datadrivenQCDMC1; 
      chisquarevalue += (Nev_data1-Nev_exp1)*(Nev_data1-Nev_exp1) / Nev_data1;
    }
    for(int ibin=1; ibin<=data2_->GetXaxis()->GetNbins(); ++ibin) {
      double Nev_data2 = data2_->GetBinContent(ibin);
      double Nev_constMC2 = constMC2_->GetBinContent(ibin);
      double Nev_diphotonMC2 = SFdiphoton_ * diphotonMC2_->GetBinContent(ibin);
      double Nev_datadrivenQCDMC2 = SFdatadrivenQCD_ * datadrivenQCD2_->GetBinContent(ibin);
      double Nev_exp2 = Nev_constMC2+Nev_diphotonMC2+Nev_datadrivenQCDMC2; 
      chisquarevalue += (Nev_data2-Nev_exp2)*(Nev_data2-Nev_exp2) / Nev_data2;
    }

    return chisquarevalue;
  }
  
private:
  TH1F *data1_, *constMC1_, *diphotonMC1_, *datadrivenQCD1_;
  TH1F *data2_, *constMC2_, *diphotonMC2_, *datadrivenQCD2_;
  double SFdiphoton_, SFdatadrivenQCD_;
};

void deriveSF()
{

  // Inputs//(LeadPhoton_pt/CMS_hgg_mass>0.35)*(SubleadPhoton_pt/CMS_hgg_mass>0.25)*
  string common_selection = "((Diphoton_mass<115 || Diphoton_mass>135) && ( Diphoton_minID > -0.9))";
  
  map<string, vector<string> > filenames_map;
  map<string, vector<string> > treenames_map;
  map<string, vector<float> > lumis_map;

  filenames_map["diphoton"] = { vector<string>{
      "/eos/user/z/zhenxuan/wwyy/root_files/UL17_DiphotonJetsBox_wocut.root"
    }};
  filenames_map["data"] = { vector<string>{
      "/eos/user/z/zhenxuan/wwyy/root_files/Data_cat1_wocut.root"
    }};
  filenames_map["qcd"] = { vector<string>{
      "/eos/user/z/zhenxuan/wwyy/root_files/DataDriven_data_cat1.root"
    }};

  treenames_map["diphoton"] = { vector<string>{
      "cat1"
    }};
  treenames_map["data"] = { vector<string>{
      "cat1"
    }};
  treenames_map["qcd"] = { vector<string>{
      "cat1"
    }};

  lumis_map["diphoton"] = {vector<float>{
    1.
  }};
  lumis_map["data"] = {vector<float>{
    1.
  }};
  lumis_map["qcd"] = {vector<float>{
    1.
  }};

  //Get min and max photon ID distribution
  map<string,TH1F*> h_Diphoton_minID;
  map<string,TH1F*> h_Diphoton_maxID;
  for( auto filenameitr : filenames_map) {
    auto samplename = filenameitr.first;
    auto filenames = filenameitr.second;
    auto treenames = treenames_map[samplename];
    auto lumis = lumis_map[samplename];
    h_Diphoton_minID[samplename] = new TH1F(Form("h_Diphoton_minID_%s",samplename.c_str()),
					 Form("h_Diphoton_minID_%s",samplename.c_str()),
					 34,-0.9,1);
    h_Diphoton_maxID[samplename] = new TH1F(Form("h_Diphoton_maxID_%s",samplename.c_str()),
					 Form("h_Diphoton_maxID_%s",samplename.c_str()),
					 34,-0.9,1);

    for( unsigned ifile=0; ifile<filenames.size(); ++ifile) {
      TChain* ch = new TChain();
      ch->Add( Form("%s/%s",filenames[ifile].c_str(),treenames[ifile].c_str()) );
      ch->Draw( Form("Diphoton_minID >>+ h_Diphoton_minID_%s",samplename.c_str()),
      Form("(weight_central)*(%f)*(%s)",lumis[ifile],common_selection.c_str()),
      "goff");
        ch->Draw( Form("Diphoton_maxID >>+ h_Diphoton_maxID_%s",samplename.c_str()),
      Form("(weight_central)*(%f)*(%s)",lumis[ifile],common_selection.c_str()),
      "goff");

      
      cout<< "h_Diphoton_minID" << samplename <<  h_Diphoton_minID[samplename]->GetEntries()<<endl;
      cout<<"h_Diphoton_maxID" << samplename << h_Diphoton_maxID[samplename]->GetEntries()<<endl;
    }
  }
  
  // For now assume no other MC --> h_Diphoton_minID["otherMC"] and h_Diphoton_maxID["otherMC"] are left empty
  // this can be updated to further improve the data/MC agreement of few percent
  h_Diphoton_minID["otherMC"] = new TH1F("constMC1","constMC1",34,-0.9,1);
  h_Diphoton_minID["MCtot"] = new TH1F("MCtot1","MCtot1",34,-0.9,1);
  h_Diphoton_minID["scaledMCtot"] = new TH1F("MCtot1_scaled","MCtot1_scaled",34,-0.9,1);
  h_Diphoton_maxID["otherMC"] = new TH1F("constMC2","constMC2",34,-0.9,1);
  h_Diphoton_maxID["MCtot"] = new TH1F("MCtot2","MCtot2",34,-0.9,1);
  h_Diphoton_maxID["scaledMCtot"] = new TH1F("MCtot2_scaled","MCtot2_scaled",34,-0.9,1);

  TCanvas* c1 = new TCanvas();
  h_Diphoton_minID["MCtot"]->Add(h_Diphoton_minID["otherMC"]);
  h_Diphoton_minID["MCtot"]->Add(h_Diphoton_minID["diphoton"]);
  h_Diphoton_minID["MCtot"]->Add(h_Diphoton_minID["qcd"]);
  h_Diphoton_minID["MCtot"]->Draw("hist");
  h_Diphoton_minID["data"]->SetMarkerStyle(20);
  h_Diphoton_minID["data"]->Draw("E1 same");
  h_Diphoton_minID["MCtot"]->GetYaxis()->SetRangeUser(0., 1.3*TMath::Max(h_Diphoton_minID["MCtot"]->GetMaximum(),h_Diphoton_minID["data"]->GetMaximum()) );
  c1->SaveAs("Diphoton_minID.png");

  TCanvas* c12 = new TCanvas();
  h_Diphoton_maxID["MCtot"]->Add(h_Diphoton_maxID["otherMC"]);
  h_Diphoton_maxID["MCtot"]->Add(h_Diphoton_maxID["diphoton"]);
  h_Diphoton_maxID["MCtot"]->Add(h_Diphoton_maxID["qcd"]);
  h_Diphoton_maxID["MCtot"]->Draw("hist");
  h_Diphoton_maxID["data"]->SetMarkerStyle(20);
  h_Diphoton_maxID["data"]->Draw("E1 same");
  h_Diphoton_maxID["MCtot"]->GetYaxis()->SetRangeUser(0., 1.3*TMath::Max(h_Diphoton_maxID["MCtot"]->GetMaximum(),h_Diphoton_maxID["data"]->GetMaximum()) );

  chisquare chisquareobj(h_Diphoton_minID["data"], h_Diphoton_minID["otherMC"], h_Diphoton_minID["diphoton"], h_Diphoton_minID["qcd"], 
			 h_Diphoton_maxID["data"], h_Diphoton_maxID["otherMC"], h_Diphoton_maxID["diphoton"], h_Diphoton_maxID["qcd"]);
  TF2 *f = new TF2("chi2",chisquareobj,0.1,2.,0.1,2.,0);
  
  double SFdiphoton,SFdatadrivenQCD;
  double chi2 = f->GetMinimumXY(SFdiphoton,SFdatadrivenQCD);
  cout<<"observed: "<<SFdiphoton<<" "<<SFdatadrivenQCD<<" chi2="<<chi2<<endl;
  c12->SaveAs("Diphoton_maxID.png");
  TCanvas* c2 = new TCanvas();
  h_Diphoton_minID["diphoton"]->Scale(SFdiphoton);
  h_Diphoton_minID["qcd"]->Scale(SFdatadrivenQCD);
  h_Diphoton_minID["scaledMCtot"]->Add(h_Diphoton_minID["otherMC"]);
  h_Diphoton_minID["scaledMCtot"]->Add(h_Diphoton_minID["diphoton"]);
  h_Diphoton_minID["scaledMCtot"]->Add(h_Diphoton_minID["qcd"]);
  h_Diphoton_minID["scaledMCtot"]->Draw("hist");
  h_Diphoton_minID["data"]->Draw("E1 same");
  h_Diphoton_minID["scaledMCtot"]->GetYaxis()->SetRangeUser(0., 1.3*TMath::Max(h_Diphoton_minID["scaledMCtot"]->GetMaximum(),h_Diphoton_minID["data"]->GetMaximum()) );
  c2->SaveAs("Diphoton_minID_scale.png");
  TCanvas* c22 = new TCanvas();
  h_Diphoton_maxID["diphoton"]->Scale(SFdiphoton);
  h_Diphoton_maxID["qcd"]->Scale(SFdatadrivenQCD);
  h_Diphoton_maxID["scaledMCtot"]->Add(h_Diphoton_maxID["otherMC"]);
  h_Diphoton_maxID["scaledMCtot"]->Add(h_Diphoton_maxID["diphoton"]);
  h_Diphoton_maxID["scaledMCtot"]->Add(h_Diphoton_maxID["qcd"]);
  h_Diphoton_maxID["scaledMCtot"]->Draw("hist");
  h_Diphoton_maxID["data"]->Draw("E1 same");
  cout<<"QCD:"<<h_Diphoton_maxID["qcd"]->Integral()<<endl;
  cout<<"DiPhoton:"<<h_Diphoton_maxID["diphoton"]->Integral()<<endl;
  cout<<"data:"<<h_Diphoton_maxID["data"]->Integral()<<endl;
  h_Diphoton_maxID["scaledMCtot"]->GetYaxis()->SetRangeUser(0., 1.3*TMath::Max(h_Diphoton_maxID["scaledMCtot"]->GetMaximum(),h_Diphoton_maxID["data"]->GetMaximum()) );
  c22->SaveAs("Diphoton_maxID_scale.png");
}

