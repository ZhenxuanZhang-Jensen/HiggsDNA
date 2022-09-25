import glob
import awkward
import numpy            
from higgs_dna.utils import awkward_utils           
           
merged_events = []
files = glob.glob("/eos/user/s/shsong/diphojetbox_mgg_passedjetcut/DiPhotonJetsBox_MGG_2017/job_*/output_job*.parquet")
for f in files:
    merged_events.append(awkward.from_parquet(f))
branches = []
for events in merged_events:
    branches +=events.fields
branches = set(branches)

for idx, events in enumerate(merged_events):
    missing_branches = [b for b in branches if b not in events.fields]
    for b in missing_branches:
        awkward_utils.add_field(
                events = merged_events[idx],
                name = b,
                data = numpy.ones(len(merged_events[idx]))
        )

logger.info("[JobsManager : merge_outputs] For syst_tag '%s' merging %d files into merged file '%s'." % (syst_tag, len(outputs), merged_file))

merged_events = awkward.concatenate(merged_events)
            
logger.debug("\t merging %d events" % len(merged_events))
awkward.to_parquet(merged_events, "/eos/user/s/shsong/diphojetbox_mgg_passedjetcut/DiPhotonJetsBox_MGG_2017/merged_nominal.parquet")

