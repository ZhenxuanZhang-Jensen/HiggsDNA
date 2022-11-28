import json
import os
from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.analysis import run_analysis

logger = setup_logger('DEBUG')
config_file = '/afs/cern.ch/user/s/shsong/HiggsDNA/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_2017/job_3/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_2017_config_job3.json'
if not os.path.exists(config_file):
    config_file = os.path.split(config_file)[-1]
with open(config_file, 'r') as f_in:
    config = json.load(f_in)

run_analysis(config)
