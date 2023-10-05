#!/bin/bash
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Set flag
flag="nothing"
# Get the options
while getopts ":hf:" option; do
   case $option in
      f) # Enter a name
         flag=$OPTARG;;
   esac
done

if [ "$flag" = "nothing" ]; then
   echo "running $flag"
fi
if [ "$flag" = "all_signal_FH" ]; then
   echo "running $flag"
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_FH_custom_cat1/*.txt  -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_FH_custom_cat1/*.json  -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_FH_custom_cat1/*.pkl  -rf
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_new_FH_cat   -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_M-3000"  --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_new_FH_cat" --batch_system "condor"
   # "UL17_R_gghh_M-1000","UL17_R_gghh_M-1700","UL17_R_gghh_M-1800","UL17_R_gghh_M-2200","UL17_R_gghh_M-250","UL17_R_gghh_M-2600","UL17_R_gghh_M-260","UL17_R_gghh_M-280","UL17_R_gghh_M-3000","UL17_R_gghh_M-300","UL17_R_gghh_M-320","UL17_R_gghh_M-350","UL17_R_gghh_M-400","UL17_R_gghh_M-450","UL17_R_gghh_M-550","UL17_R_gghh_M-600","UL17_R_gghh_M-650","UL17_R_gghh_M-700","UL17_R_gghh_M-750","UL17_R_gghh_M-800","UL17_R_gghh_M-850","UL17_R_gghh_M-900","UL17_R_gghh_M-1600","UL17_R_gghh_M-2400","UL17_R_gghh_M-1200","UL17_R_gghh_M-1900","UL17_R_gghh_M-1500","UL17_R_gghh_M-1100","UL17_R_gghh_M-1400","UL17_R_gghh_M-2800","UL17_R_gghh_M-270","UL17_R_gghh_M-2000","UL17_R_gghh_M-1300"
# "UL17_R_gghh_M-1000","UL17_R_gghh_M-1700","UL17_R_gghh_M-1800","UL17_R_gghh_M-2200","UL17_R_gghh_M-250","UL17_R_gghh_M-2600","UL17_R_gghh_M-260","UL17_R_gghh_M-280","UL17_R_gghh_M-3000","UL17_R_gghh_M-300","UL17_R_gghh_M-320","UL17_R_gghh_M-350","UL17_R_gghh_M-400","UL17_R_gghh_M-450","UL17_R_gghh_M-550","UL17_R_gghh_M-650","UL17_R_gghh_M-700","UL17_R_gghh_M-750","UL17_R_gghh_M-800","UL17_R_gghh_M-850","UL17_R_gghh_M-900","UL17_R_gghh_M-600","UL17_R_gghh_M-1600","UL17_R_gghh_M-2400","UL17_R_gghh_M-1200","UL17_R_gghh_M-1900","UL17_R_gghh_M-1500","UL17_R_gghh_M-1100","UL17_R_gghh_M-1400","UL17_R_gghh_M-2800","UL17_R_gghh_M-270","UL17_R_gghh_M-2000","UL17_R_gghh_M-1300" #  --merge_outputs
fi
if [ "$flag" = "frac_signal_FH" ]; then
   echo "running $flag"
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_FH_individual_FH_limit -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_sep_FH_SL_withPN.json" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_M-2000"  --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_signal_FH_individual_FH_limit/" --batch_system "condor"
fi
if [ "$flag" = "frac_signal_SL" ]; then
   echo "running $flag"
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_signal_SL_individual_SL_limit -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_sep_FH_SL_withPN.json" --sample_list "UL17_R_gghh_SL_M-3000","UL17_R_gghh_SL_M-2000"  --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_signal_SL_individual_SL_limit/" --batch_system "condor"
fi
if [ "$flag" = "all_signal_SL" ]; then
   echo "running $flag"
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_signal_new_SL_cat -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list   "UL17_R_gghh_SL_M-1200","UL17_R_gghh_SL_M-2600","UL17_R_gghh_SL_M-3000","UL17_R_gghh_SL_M-700","UL17_R_gghh_SL_M-1900","UL17_R_gghh_SL_M-800","UL17_R_gghh_SL_M-1600","UL17_R_gghh_SL_M-1800","UL17_R_gghh_SL_M-1100","UL17_R_gghh_SL_M-900","UL17_R_gghh_SL_M-1000","UL17_R_gghh_SL_M-350","UL17_R_gghh_SL_M-250","UL17_R_gghh_SL_M-550","UL17_R_gghh_SL_M-1300","UL17_R_gghh_SL_M-600","UL17_R_gghh_SL_M-270","UL17_R_gghh_SL_M-300","UL17_R_gghh_SL_M-650","UL17_R_gghh_SL_M-450","UL17_R_gghh_SL_M-750","UL17_R_gghh_SL_M-280","UL17_R_gghh_SL_M-2200","UL17_R_gghh_SL_M-2000","UL17_R_gghh_SL_M-2800","UL17_R_gghh_SL_M-1500","UL17_R_gghh_SL_M-850","UL17_R_gghh_SL_M-260","UL17_R_gghh_SL_M-1400","UL17_R_gghh_SL_M-320","UL17_R_gghh_SL_M-1700" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_signal_new_SL_cat" --batch_system "condor"
   # "UL17_R_gghh_SL_M-1200","UL17_R_gghh_SL_M-2400","UL17_R_gghh_SL_M-400","UL17_R_gghh_SL_M-2600","UL17_R_gghh_SL_M-3000","UL17_R_gghh_SL_M-700","UL17_R_gghh_SL_M-1900","UL17_R_gghh_SL_M-800","UL17_R_gghh_SL_M-1600","UL17_R_gghh_SL_M-1800","UL17_R_gghh_SL_M-1100","UL17_R_gghh_SL_M-900","UL17_R_gghh_SL_M-1000","UL17_R_gghh_SL_M-350","UL17_R_gghh_SL_M-250","UL17_R_gghh_SL_M-550","UL17_R_gghht_SL_M-1300","UL17_R_gghh_SL_M-600","UL17_R_gghh_SL_M-270","UL17_R_gghh_SL_M-300","UL17_R_gghh_SL_M-650","UL17_R_gghh_SL_M-450","UL17_R_gghh_SL_M-750","UL17_R_gghh_SL_M-280","UL17_R_gghh_SL_M-2200","UL17_R_gghh_SL_M-2000","UL17_R_gghh_SL_M-2800","UL17_R_gghh_SL_M-1500","UL17_R_gghh_SL_M-850","UL17_R_gghh_SL_M-260","UL17_R_gghh_SL_M-1400","UL17_R_gghh_SL_M-320","UL17_R_gghh_SL_M-1700"
fi
if [ "$flag" = "test_signal_SL" ]; then
   echo "running $flag"
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_test_signal_SL -rf
   # no merge ouput and no condor
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_R_gghh_SL_M-2000" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_test_signal_SL" --short
fi
if [ "$flag" = "test_signal_FH" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/FH_channel/hhwwgg_test_signal_FH -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_SL_M-2000" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/FH_channel/hhwwgg_test_signal_FH" --short 
fi
if [ "$flag" = "all_data" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration/*.pkl -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration/*.txt -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration/*UL17_dataB_2017* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration/*UL17_dataC_2017* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration/*UL17_dataE_2017* -rf
   

   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_dataB" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration"  --merge_outputs --yield_table 
   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_dataB" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibrationB" --batch_system "condor" --merge_outputs --yield_table 

   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_dataC" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_data_calibration" --batch_system "condor"
fi
if [ "$flag" = "all_bkg_calibration" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_TTJets_2017 -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/*.pkl
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/*.txt
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/*.json   
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_DiPhotonJetsBox_M40_80* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_DiPhotonJetsBox_MGG_80toInf* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_GluGluHToGG* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_QCD_Pt-30to40_MGG-80toInf* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_QCD_Pt-30toInf_MGG-40to80* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_QCD_Pt-40ToInf_MGG-80ToInf* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_TTGG_0Jets* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_VBFHToGG* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_VHToGG* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_W1JetsToLNu* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_W2JetsToLNu* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_W3JetsToLNu* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_W4JetsToLNu* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_WGGJets* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_WGJJToLNu_EWK_QCD* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_WGJJToLNu_EWKnotop_QCD* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_WWTo1L1Nu2Q_4f* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_ttHJetToGG* -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/UL17_ttWJets* -rf

   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list    "WJetsToQQ_HT-200to400" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/" --batch_system "condor" 

   # "UL17_W1JetsToLNu","UL17_W2JetsToLNu","UL17_W3JetsToLNu","UL17_W4JetsToLNu","UL17_WGGJets","UL17_WGJJToLNu_EWK_QCD","UL17_WWTo1L1Nu2Q_4f","UL17_ttHJetToGG","UL17_ttWJets"    
   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf","UL17_GluGluHToGG","UL17_QCD_Pt-30to40_MGG-80toInf","UL17_QCD_Pt-30toInf_MGG-40to80","UL17_QCD_Pt-40ToInf_MGG-80ToInf","UL17_TTGG_0Jets","UL17_VBFHToGG","UL17_VHToGG","UL17_W1JetsToLNu","UL17_W2JetsToLNu","UL17_W3JetsToLNu","UL17_W4JetsToLNu","UL17_WGGJets","UL17_WGJJToLNu_EWK_QCD","UL17_WGJJToLNu_EWKnotop_QCD","UL17_WWTo1L1Nu2Q_4f","UL17_ttHJetToGG","UL17_ttWJets" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/" --batch_system "condor" 
   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_TTJets" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_test/"  --merge_outputs --yield_table --batch_system "condor"

   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_DiPhotonJetsBox_M40_80" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/" --merge_outputs --yield_table --short 

   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_sep_FH_SL_withPN.json" --sample_list "UL17_dataD","UL17_dataE","UL17_dataF","UL17_dataB","UL17_dataC" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_indiviual_FHSL" --batch_system "condor" --merge_outputs
fi
if [ "$flag" = "test_bkg_calibration" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/*.pkl -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/*.json -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/*.txt -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local -rf

   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_TTJets" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/"  --yield_table 
   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_TTJets" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_test/"  --merge_outputs --yield_table --batch_system "condor"

   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_DiPhotonJetsBox_M40_80" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration/" --merge_outputs --yield_table --short 

   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_sep_FH_SL_withPN.json" --sample_list "UL17_dataD","UL17_dataE","UL17_dataF","UL17_dataB","UL17_dataC" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/Data/hhwwgg_data_indiviual_FHSL" --batch_system "condor" --merge_outputs
fi
if [ "$flag" = "test_data" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   rm -rf /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_test_data_FHSL
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_dataD" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_test_data_FHSL" --batch_system "condor"
fi
if [ "$flag" = "cal_test_data" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   rm -rf /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_test_data_cal
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json" --sample_list "UL17_dataD" --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/Data/hhwwgg_test_data_cal" --short
fi

if [ "$flag" = "test_data_FH" ]; then
   echo "running $flag"
   # frac of data
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_test_data_FH/* -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FH.json" --sample_list "test_data" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/hhwwgg_test_data_FH" --batch_system "condor"
   # test of data
   # python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FH.json" --sample_list "Data_test" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_test_data_FH"

   #  python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection.json" --sample_list "GluGluToRadionToHHTo2G2WTo2G4Q_M-250" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_testSignal" --short
fi
if [ "$flag" = "test_data_SL" ]; then
   echo "running $flag"
   # no merge ouput and no condor
    python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_SL.json" --sample_list "DoubleEG_Run2017B" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_test_data_SL" --short
fi
if [ "$flag" = "all_data_SL" ]; then
   echo "running $flag"
   # no merge ouput and no condor
    python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_SL.json" --sample_list "DoubleEG_Run2017B","DoubleEG_Run2017C","DoubleEG_Run2017D","DoubleEG_Run2017E","DoubleEG_Run2017F" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/SL_channel/hhwwgg_data_SL" --batch_system "condor" --merge_outputs
   #  python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_SL.json" --sample_list "DoubleEG_Run2017B","DoubleEG_Run2017C","DoubleEG_Run2017D","DoubleEG_Run2017E","DoubleEG_Run2017F" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_data_SL" --merge_outputs
fi

if [ "$flag" = "all_bkg" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/*.txt -rf
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/*.pkl -rf
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/*.json -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf_2017/* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80_2017/* -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/UL17_QCD_Pt-40ToInf_MGG-80ToInf_2017 -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1/UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf -rf
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_cat1" --batch_system "condor"
   # "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf","UL17_GluGluHToGG","UL17_QCD_Pt-30to40_MGG-80toInf","UL17_QCD_Pt-30toInf_MGG-40to80","UL17_QCD_Pt-40ToInf_MGG-80ToInf","UL17_TTGG_0Jets","UL17_TTGJets","UL17_TTJets","UL17_TTToHadronic","UL17_VBFHToGG","UL17_VHToGG","UL17_W1JetsToLNu","UL17_W2JetsToLNu","UL17_W3JetsToLNu","UL17_W4JetsToLNu","UL17_WGGJets","UL17_WGJJToLNu_EWK_QCD","UL17_WGJJToLNu_EWKnotop_QCD","UL17_WWG","UL17_WWTo1L1Nu2Q_4f","UL17_ttHJetToGG","UL17_ttWJets"
fi
if [ "$flag" = "all_bkg_v9" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_v9/* -rf
    python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "DiPhotonJetsBoxMGG_80toInf","DiPhotonJetsBoxM40_80" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_bkgs_v9" --short
   #   --short
   #  --batch_system "condor"
fi
if [ "$flag" = "test_bkg" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   rm /eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_test_bkgs -rf
    python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_DiPhotonJetsBox_MGG_80toInf" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_test_bkgs" --short
   #  python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection.json" --sample_list "GluGluToRadionToHHTo2G2WTo2G4Q_M-250" --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/hhwwgg_testSignal" --short
fi

if [ "$flag" = "run_bkg_diphotonjetsbox_FH" ]; then
   echo "running $flag"
   # no merge ouput and no condor
   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FH.json" --sample_list "DiPhotonJetsBox_MGG-80toInf"  --output_dir "/eos/user/z/zhenxuan/hhwwgg_parquet/FH_channel/bkg_FH/" --batch_system "condor"
fi


 # ----------------------- calibration of PN tagger!!!! ----------------------- #
if [ "$flag" = "run_bkg_cal" ]; then 
   echo "running $flag"
   # no merge ouput and no condor
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/*.json -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/*.pkl -rf
   rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/*.txt -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/UL17_QCD_Pt-30to40_MGG-80toInf -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/UL17_QCD_Pt-30toInf_MGG-40to80 -rf
   # rm /eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/UL17_W1JetsToLNu -rf

   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_top_region_calibration.json"  --sample_list "UL17_W1JetsToLNu"  --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/bkg_calibration_local/" --short 
   
   # "UL17_R_gghh_M-1000","UL17_R_gghh_M-1700","UL17_R_gghh_M-1800","UL17_R_gghh_M-2200","UL17_R_gghh_M-250","UL17_R_gghh_M-2600","UL17_R_gghh_M-260","UL17_R_gghh_M-280","UL17_R_gghh_M-3000","UL17_R_gghh_M-300","UL17_R_gghh_M-320","UL17_R_gghh_M-350","UL17_R_gghh_M-400","UL17_R_gghh_M-450","UL17_R_gghh_M-550","UL17_R_gghh_M-600","UL17_R_gghh_M-650","UL17_R_gghh_M-700","UL17_R_gghh_M-750","UL17_R_gghh_M-800","UL17_R_gghh_M-850","UL17_R_gghh_M-900","UL17_R_gghh_M-1600","UL17_R_gghh_M-2400","UL17_R_gghh_M-1200","UL17_R_gghh_M-1900","UL17_R_gghh_M-1500","UL17_R_gghh_M-1100","UL17_R_gghh_M-1400","UL17_R_gghh_M-2800","UL17_R_gghh_M-270","UL17_R_gghh_M-2000","UL17_R_gghh_M-1300","s_channel_4f_hadronic_decay","s_channel_4f_leptonic_decay","tW_antitop_5f_inclusive","tW_top_5f_inclusive","t_channel_antitop_4f_inclusive","t_channel_top_4f_inclusive"
   # "s_channel_4f_hadronic_decay","s_channel_4f_leptonic_decay", "tW_antitop_5f_inclusive", "tW_top_5f_inclusive", "t_channel_antitop_4f_inclusive","t_channel_top_4f_inclusive"
   # "TTJets_HT-1200to2500_TuneCP5_13TeV","TTJets_HT-2500toInf_TuneCP5_13TeV","TTJets_HT-600to800_TuneCP5_13TeV","TTJets_HT-800to1200_TuneCP5_13TeV","UL17_R_gghh_M-1000","UL17_R_gghh_M-1700","UL17_R_gghh_M-1800","UL17_R_gghh_M-2200","UL17_R_gghh_M-250","UL17_R_gghh_M-2600","UL17_R_gghh_M-260","UL17_R_gghh_M-280","UL17_R_gghh_M-3000","UL17_R_gghh_M-300","UL17_R_gghh_M-320","UL17_R_gghh_M-350","UL17_R_gghh_M-400","UL17_R_gghh_M-450","UL17_R_gghh_M-550","UL17_R_gghh_M-600","UL17_R_gghh_M-650","UL17_R_gghh_M-700","UL17_R_gghh_M-750","UL17_R_gghh_M-800","UL17_R_gghh_M-850","UL17_R_gghh_M-900","UL17_R_gghh_M-1600","UL17_R_gghh_M-2400","UL17_R_gghh_M-1200","UL17_R_gghh_M-1900","UL17_R_gghh_M-1500","UL17_R_gghh_M-1100","UL17_R_gghh_M-1400","UL17_R_gghh_M-2800","UL17_R_gghh_M-270","UL17_R_gghh_M-2000","UL17_R_gghh_M-1300"
fi
# ------------------------------------ -- ------------------------------------ #

#!/usr/bin

# WhichSamples=${1}

# if [ ${WhichSamples} -eq "test" ]
#   then

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list  --output_dir "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/hhwwgg_sig_test" --batch_system "condor" --yield_table 
# fi


# if [ ${WhichSamples} -eq 0 ]
#   then
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_v1" --batch_system "condor" --merge_outputs --yield_table
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_v3" --batch_system "condor" --merge_outputs --yield_table
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL_tmp.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_v2" --batch_system "condor" --merge_outputs --yield_table
  
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v3/m3000/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v3/m3000dphi/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v3/m3000dphi_v2/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v3/m3000dphi_v3/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_v5" --batch_system "condor" --merge_outputs --yield_table


#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggFH/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_M-1000","UL17_R_gghh_M-2000","UL17_R_gghh_M-1500","UL17_R_gghh_M-600" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggFH/M800" --sample_list "UL17_R_gghh_M-800" --yield_table --batch_system "condor"

#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggSL/" --sample_list "UL17_R_gghh_SL_M-500","UL17_R_gghh_SL_M-1000","UL17_R_gghh_SL_M-2000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"

#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_v1" --batch_system "condor" --merge_outputs --yield_table
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggSL3000/" --sample_list "UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggFH/X3000Y2800" --sample_list "NMSSM_XToYHTo2G2WTo2G4Q_MX-3000_MY-2800" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggFH/X3000Y28002" --sample_list "NMSSM_XToYHTo2G2WTo2G4Q_MX-3000_MY-2800" --yield_table 

#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/WWgg/SLcheck_electron/SL3000" --sample_list "UL17_R_gghh_SL_M-3000" --yield_table --short
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_latest/check_electron" --yield_table --short
# fi
# if [ ${WhichSamples} -eq 1 ]
#   then
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggFH/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_M-1000","UL17_R_gghh_M-2000","UL17_R_gghh_M-1500","UL17_R_gghh_M-800","UL17_R_gghh_M-550","UL17_R_gghh_M-600" --yield_table --batch_system "condor"

#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/WWggSL/" --sample_list "UL17_R_gghh_SL_M-500","UL17_R_gghh_SL_M-600","UL17_R_gghh_SL_M-800","UL17_R_gghh_SL_M-1000","UL17_R_gghh_SL_M-1500","UL17_R_gghh_SL_M-2000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
  
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/bkg_v1/" --sample_list "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf","UL17_GluGluHToGG","UL17_QCD_Pt-30to40_MGG-80toInf","UL17_QCD_Pt-30toInf_MGG-40to80","UL17_QCD_Pt-40ToInf_MGG-80ToInf","UL17_TTGG_0Jets","UL17_TTGJets","UL17_TTJets","UL17_TTToHadronic","UL17_W1JetsToLNu","UL17_W2JetsToLNu","UL17_W3JetsToLNu","UL17_W4JetsToLNu","UL17_WGGJets","UL17_WGJJToLNu_EWK_QCD","UL17_WWG","UL17_WWTo1L1Nu2Q_4f","UL17_ttWJets" --yield_table --batch_system "condor"

#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL_highpt.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v2/WWggFH/" --sample_list "UL17_R_gghh_M-3000","UL17_R_gghh_M-1000","UL17_R_gghh_M-2000","UL17_R_gghh_M-1500","UL17_R_gghh_M-800","UL17_R_gghh_M-550","UL17_R_gghh_M-600" --yield_table --batch_system "condor"
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL_highpt.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v2/WWggSL/" --sample_list "UL17_R_gghh_SL_M-500","UL17_R_gghh_SL_M-600","UL17_R_gghh_SL_M-800","UL17_R_gghh_SL_M-1000","UL17_R_gghh_SL_M-1500","UL17_R_gghh_SL_M-2000","UL17_R_gghh_SL_M-3000" --yield_table --batch_system "condor"
  
# fi
# if [ ${WhichSamples} -eq 2]
#   then
#   python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_M-1000","UL17_R_gghh_M-1700","UL17_R_gghh_M-1800","UL17_R_gghh_M-2200","UL17_R_gghh_M-250","UL17_R_gghh_M-2600","UL17_R_gghh_M-260","UL17_R_gghh_M-280","UL17_R_gghh_M-3000","UL17_R_gghh_M-300","UL17_R_gghh_M-320","UL17_R_gghh_M-350","UL17_R_gghh_M-400","UL17_R_gghh_M-450","UL17_R_gghh_M-550","UL17_R_gghh_M-600","UL17_R_gghh_M-650","UL17_R_gghh_M-700","UL17_R_gghh_M-750","UL17_R_gghh_M-800","UL17_R_gghh_M-850","UL17_R_gghh_M-900","UL17_R_gghh_M-1600","UL17_R_gghh_M-2400","UL17_R_gghh_M-1200","UL17_R_gghh_M-1900","UL17_R_gghh_M-1500","UL17_R_gghh_M-1100","UL17_R_gghh_M-1400","UL17_R_gghh_M-2800","UL17_R_gghh_M-270","UL17_R_gghh_M-2000","UL17_R_gghh_M-1300" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/WWgg/FH" --yield_table --batch_system "condor"
# fi
# if [ ${WhichSamples} -eq 3]
#   then

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf","UL17_GluGluHToGG","UL17_QCD_Pt-30to40_MGG-80toInf","UL17_QCD_Pt-30toInf_MGG-40to80","UL17_QCD_Pt-40ToInf_MGG-80ToInf","UL17_TTGG_0Jets","UL17_TTGJets","UL17_TTJets","UL17_TTToHadronic","UL17_VBFHToGG","UL17_VHToGG","UL17_W1JetsToLNu","UL17_W2JetsToLNu","UL17_W3JetsToLNu","UL17_W4JetsToLNu","UL17_WGGJets","UL17_WGJJToLNu_EWK_QCD","UL17_WWG","UL17_WWTo1L1Nu2Q_4f","UL17_ttHJetToGG","UL17_ttWJets","UL17_DY1JetsToLL_M-10to50","UL17_DY2JetsToLL_M-10to50","UL17_DY3JetsToLL_M-10to50","UL17_DY4JetsToLL_M-10to50","UL17_DYJetsToLL_0J","UL17_DYJetsToLL_1J","UL17_DYJetsToLL_2J","UL17_DYJetsToLL_M-50" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/bkg_latest/" --batch_system "condor" --yield_table 
# fi

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL_official.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig/cat3/m500" --sample_list "UL17_R_gghh_SL_M-500" --batch_system "condor"
# #cat2 datadriven
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/datadriven/cat2" --sample_list "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf" --batch_system "condor"
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/datadriven/cat2/data" --batch_system "condor" --merge_outputs
# ########################
# #run ZZgg sample
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig/ZZgg" --sample_list "GluGluToRadionToHHTo2G2ZTo2G4Q_M-250_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-260_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-270_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-280_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-300_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-320_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-350_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-400_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-450_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-500_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-550_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-600_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-650_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-700_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-750_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-800_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-900_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1000_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1100_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1200_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1300_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1400_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1500_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1600_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1700_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1800_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-1900_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-2000_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-2200_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-2400_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-2600_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-2800_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2ZTo2G4Q_M-3000_TuneCP5_PSWeights_narrow_13TeV" --batch_system "condor"
# #run ttgg sample
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig/ttgg" --sample_list "GluGluToRadionToHHTo2G2Tau_M-250_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-260_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-270_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-280_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-290_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-300_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-320_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-350_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-400_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-450_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-500_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-550_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-600_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-650_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-700_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-750_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-800_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-900_narrow_TuneCP5_13TeV","GluGluToRadionToHHTo2G2Tau_M-1000_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1100_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1200_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1300_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1400_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1500_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1600_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1700_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1800_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-1900_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-2000_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-2200_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-2400_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-2600_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-2800_TuneCP5_PSWeights_narrow_13TeV","GluGluToRadionToHHTo2G2Tau_M-3000_TuneCP5_PSWeights_narrow_13TeV" --batch_system "condor"

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/WWgg" --sample_list "UL17_R_gghh_SL_M-1200","UL17_R_gghh_SL_M-2400","UL17_R_gghh_SL_M-400","UL17_R_gghh_SL_M-2600","UL17_R_gghh_SL_M-3000","UL17_R_gghh_SL_M-700","UL17_R_gghh_SL_M-1900","UL17_R_gghh_SL_M-800","UL17_R_gghh_SL_M-1600","UL17_R_gghh_SL_M-1800","UL17_R_gghh_SL_M-1100","UL17_R_gghh_SL_M-900","UL17_R_gghh_SL_M-1000","UL17_R_gghh_SL_M-350","UL17_R_gghh_SL_M-250","UL17_R_gghh_SL_M-550","UL17_R_gghh_SL_M-1300","UL17_R_gghh_SL_M-600","UL17_R_gghh_SL_M-270","UL17_R_gghh_SL_M-300","UL17_R_gghh_SL_M-650","UL17_R_gghh_SL_M-450","UL17_R_gghh_SL_M-750","UL17_R_gghh_SL_M-280","UL17_R_gghh_SL_M-2200","UL17_R_gghh_SL_M-2000","UL17_R_gghh_SL_M-2800","UL17_R_gghh_SL_M-1500","UL17_R_gghh_SL_M-850","UL17_R_gghh_SL_M-260","UL17_R_gghh_SL_M-1400","UL17_R_gghh_SL_M-320","UL17_R_gghh_SL_M-1700" --batch_system "condor"

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB","UL17_dataC","UL17_dataD","UL17_dataE","UL17_dataF" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_latest/" --batch_system "condor" --merge_outputs 



# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_M-3000" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/test4/" --short --yield_table 
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_SL_M-3000" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/WWggSL3000_PFiso/" --yield_table --short

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_SL_M-3000" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/test/" --yield_table --short
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_SL_M-1000" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/test2/" --yield_table --short
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_data_FHSL.json" --sample_list "UL17_dataB" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/data_latest/test_data/" --yield_table --short
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_bkg_FHSL.json" --sample_list "UL17_DiPhotonJetsBox_M40_80","UL17_DiPhotonJetsBox_MGG_80toInf","UL17_GJet_Pt_20to40_DoubleEMEnriched_MGG_80toInf","UL17_GJet_Pt_20toInf_DoubleEMEnriched_MGG_40to80","UL17_GJet_Pt_40toInf_DoubleEMEnriched_MGG_80toInf","UL17_QCD_Pt-30to40_MGG-80toInf","UL17_QCD_Pt-30toInf_MGG-40to80","UL17_QCD_Pt-40ToInf_MGG-80ToInf","UL17_TTGG_0Jets","UL17_TTGJets","UL17_TTJets","UL17_TTToHadronic","UL17_W1JetsToLNu","UL17_W2JetsToLNu","UL17_W3JetsToLNu","UL17_W4JetsToLNu","UL17_WGGJets","UL17_WGJJToLNu_EWK_QCD","UL17_WWG","UL17_WWTo1L1Nu2Q_4f","UL17_ttWJets" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/bkg_latest/test" --yield_table --short

# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_DiPhotonJetsBox_MGG_80toInf" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/test2/" --yield_table --short
# python scripts/run_analysis.py --log-level "DEBUG" --config "metadata/tutorial/HHWW_preselection_FHSL.json" --sample_list "UL17_R_gghh_M-3000" --output_dir "/eos/user/s/shsong/combined_WWgg/parquet/sig_v1/test/" --yield_table --short