# python run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/tth_presel_withSyst.json" --merge_outputs --output_dir "/eos/user/z/zhenxuan/tutorial_tth_withSyst" --batch_system "condor"
# python run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/tth_preselection.json" --merge_outputs --sample_list "ttH_M125" --output_dir "tutorial_tth" --short
python run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection.json" --merge_outputs --sample_list "ttH_M125" --output_dir "/eos/user/z/zhenxuan/hhwwgg" --batch_system "condor"
