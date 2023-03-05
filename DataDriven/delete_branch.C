void delete_branch() {
    TFile f("/eos/user/z/zhenxuan/hhwwgg_root/combined_data_cat1.root","update");
    TTree *T = (TTree*)f.Get("cat1");
    TBranch *b = T->GetBranch("weight_central");
    T->GetListOfBranches()->Remove(b);
    T->Write();
}