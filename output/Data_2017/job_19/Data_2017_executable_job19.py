import json
import os
from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.analysis import run_analysis

logger = setup_logger('DEBUG')
config_file = '/afs/cern.ch/user/z/zhenxuan/HiggsDNA/output/Data_2017/job_19/Data_2017_config_job19.json'
if not os.path.exists(config_file):
    config_file = os.path.split(config_file)[-1]
with open(config_file, 'r') as f_in:
    config = json.load(f_in)

run_analysis(config)
