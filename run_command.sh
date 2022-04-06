# python bonus/assess.py --input_dir "/eos/user/z/zhenxuan/hhwwgg" --plots "/eos/user/z/zhenxuan/hhwwgg/plots.json"
python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection.json" --merge_outputs --sample_list "GluGluToRadionToHHTo2G4Q_M300" --output_dir "/eos/user/z/zhenxuan/hhwwgg" 
# --batch_system "condor"
