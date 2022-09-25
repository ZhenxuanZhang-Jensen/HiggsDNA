import glob
import awkward

merged_events = []
files = glob.glob("QCD_HT500to700_2017/job_*/output_job*.parquet")
for f in files:
    merged_events.append(awkward.from_parquet(f))
merged_events = awkward.concatenate(merged_events)

awkward.to_parquet(merged_events, "QCD_HT500to700_2017/merged_nominal.parquet")