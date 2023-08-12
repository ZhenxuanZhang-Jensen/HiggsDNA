import glob
import awkward

merged_events = []
files = glob.glob("/eos/user/s/shsong/HHWWgg/parquet/bkg/DiphotonJets/UL17_DiPhotonJetsBox_MGG_80toInf_2017/job_*/output*.parquet")
for f in files:
    merged_events.append(awkward.from_parquet(f))
merged_events = awkward.concatenate(merged_events)

awkward.to_parquet(merged_events, "/eos/user/s/shsong/HHWWgg/parquet/bkg/DiphotonJets/UL17_DiPhotonJetsBox_MGG_80toInf_2017/merged_nominal_all.parquet")