# -*- coding: utf-8 -*-
'''
Hi my friend, believe it or not, I just want to create a king script to create limit with flashggfinal fit framework from ntuples(root) with only one script
let see if the magic can work
'''
import codecs
import uproot 
import awkward as ak
import json
from collections import defaultdict
# import ROOT
import os 
import logging
import sys
import subprocess
import time
logging.basicConfig(filename='/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/logging_output.log', level=logging.INFO, filemode="w")
logging.debug('This message should go to the log file')
logging.info('So should this')
logging.warning('And this, too')
# ----------------------------  trees2WS part ---------------------------- #
def run_Tree2WS_sig(inputpath_name, inputfile_name, log_file_name, ws_path, output_sig_root_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Trees2WS")
    # no need to change the config file since all are auto set
    command = " python trees2ws.py --inputConfig config_simple.py --inputTreeFile " + inputpath_name+inputfile_name + " --inputMass 125 --productionMode gghh --year 2017  > " + log_file_name+ " 2>&1"
    logging.info("the Tree2WS command:")
    logging.info(command)
    # run trees2ws at shell
    run_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("end: {}".format(time.time()))
    # create the dir for this workspace
    command = "mkdir " +inputpath_name+ "ws_gghh_" + ws_path
    mkdir_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("mkdir for ws \n command :{0}".format(command))
    # cp the root file to the ws dir
    command = "cp " + inputpath_name +"ws_gghh/" + inputfile_name.split('.root')[0] + '_gghh' + ".root" + " " + inputpath_name +"ws_gghh_"+ws_path +"/" + output_sig_root_name
    cp_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("cp to ws dir \n command: {0}".format(command))
def run_Tree2WS_data(inputpath_name, ws_data_path, output_data_root_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Trees2WS")
    # run data tree2ws
    command = "python trees2ws_data.py --inputConfig config_simple.py --inputTreeFile " + inputpath_name + output_data_root_name
    run_data_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    # mkdir 
    command = "mkdir " +inputpath_name+  ws_data_path
    mkdir_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    # run cp data 
    command = "cp " + inputpath_name +"ws/" + output_data_root_name + " " + inputpath_name + ws_data_path +"/" + "allData.root"
    cp_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("cp to ws dir \n command: {0}".format(command))

# ------------------------------ Signal fit part ----------------------------- #
def run_ftest(ws_path, log_name, inputpath_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Signal")
    # modify the config file
    #copy config_toy.py
    command = "cp config_toy.py " + "config_" + ws_path + ".py"
    run_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("copy config_toy.py \n command :{0}".format(command))
    time.sleep(2)
    # sed config_*.py
    command1 = 'sed -i "s#ws_path#' + ws_path + '#g" ' +  "config_" + ws_path + ".py" 
    sed1_p = subprocess.call(command1, shell=True)
    logging.info("sed1 config_*.py \n command :{0}".format(command1))

    command2 = 'sed -i "s#input_path#' + inputpath_name + '#g" ' + "config_" + ws_path + ".py"

    sed2_p = subprocess.call(command2, shell=True)
    logging.info("sed2 config_*.py \n command :{0}".format(command2))
    # run ftest
    command = "python RunSignalScripts.py --inputConfig " +  "config_" + ws_path + ".py"+ " --mode 'fTest'" + " > " + log_name + " 2>&1"
    run_ftest_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run ftest \n command :{0}".format(command))
    logging.info("end: {}".format(time.time()))

def run_signalfit(ws_path, log_name, inputpath_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Signal")
    # run signalfit
    command = "python RunSignalScripts.py --inputConfig " +  "config_" + ws_path + ".py"+ " --mode 'signalFit' --modeOpts '--skipSystematics' " + " > " + log_name + " 2>&1"
    logging.info("run signalfit \n command :{0}".format(command))
    run_signalfit_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    # cp the signalfit output root file to ws
    command = "cp outdir_dcb_2017_" + ws_path + "/signalFit/output/*.root" + " " + inputpath_name +  "ws_gghh_"+ws_path
    cp_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run cp the signalfit output root file to ws \n command :{0}".format(command))

    logging.info("end: {}".format(time.time()))

def run_signal_plot(outputExt, cats, exts,log_packaged_name,ws_path,inputpath_name,log_plotter_name, cp_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Signal")
    # run packaged
    command = "python RunPackager.py --cats " + cats + " --exts " + exts + "  --batch local  --massPoints 125 --year 2017 --outputExt " + outputExt + " > " + log_packaged_name + " 2>&1"
    run_packaged_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run packaged \n command :{0}".format(command))
    # run plotter
    command = "python RunPlotter.py --cats " + cats + " --procs all --years 2017 --ext " + outputExt + " > " + log_plotter_name + " 2>&1"
    run_plotter_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run plotter \n command :{0}".format(command))
    # cp the signalfit output root file to ws
    command = "cp outdir_packaged_" + ws_path + "/*.root" + " " + inputpath_name +  "ws_gghh_"+ ws_path + "/" + cp_name
    cp_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run cp the signalfit output root file to ws \n command :{0}".format(command))
    logging.info("end: {}".format(time.time()))
# ---------------------------- Background fit part --------------------------- #
def run_backgroundfit(ws_data_path, log_name, inputpath_name, ext_name, cp_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Background")
    # modify the config file
    #copy config_toy.py
    command = "cp config_toy.py " + "config_" + ws_data_path + ".py"
    run_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("copy config_toy.py \n command :{0}".format(command))
    time.sleep(2)
    # sed config_*.py
    command1 = 'sed -i "s#ws_path#' + ws_data_path + '#g" ' +  "config_" + ws_data_path + ".py" 
    sed1_p = subprocess.call(command1, shell=True)
    logging.info("sed1 config_*.py \n command :{0}".format(command1))

    command2 = 'sed -i "s#input_path#' + inputpath_name + '#g" ' + "config_" + ws_data_path + ".py"

    sed2_p = subprocess.call(command2, shell=True)
    logging.info("sed2 config_*.py \n command :{0}".format(command2))

    command3 = 'sed -i "s#ext_name#' + ext_name + '#g" ' + "config_" + ws_data_path + ".py"

    sed2_p = subprocess.call(command3, shell=True)
    logging.info("sed2 config_*.py \n command :{0}".format(command3))
    # run ftest
    command =  "python RunBackgroundScripts.py --inputConfig " +   "config_" + ws_data_path + ".py" + " --mode fTestParallel" + " > " + log_name + " 2>&1"
    run_ftest_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run ftest \n command :{0}".format(command))
    # run cp data 
    command = "cp outdir_" + ext_name + "/*.root" +  " " + inputpath_name + ws_data_path +"/" + cp_name
    cp_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("end: {}".format(time.time()))
 
# ---------------------------- data card part --------------------------- #
def run_yields(ws_sig_path,ws_bkg_path, inputpath_name,log_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Datacard")
    # run Runyields
    command = "python RunYields.py --inputWSDirMap 2017=" +inputpath_name+ "ws_gghh_" + ws_sig_path+" --sigModelWSDir " + inputpath_name+ "ws_gghh_" + ws_sig_path+" --bkgModelWSDir " + inputpath_name+  ws_bkg_path + " --cats auto --procs auto --batch local --ext "+ ws_sig_path+ ">" + log_name+ " 2>&1"

    run_yields_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run yields \n command :{0}".format(command))
    logging.info("end: {}".format(time.time()))

def run_makeDatacard(ws_sig_path,output_card_name,log_name,channel):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Datacard")
    # run RunmakeDatacard
      
    command = "python makeDatacard.py --years 2017 --prune --ext " + ws_sig_path + " --output "+ output_card_name+  ">" + log_name+ " 2>&1"

    run_makeDatacard_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run makeDatacard \n command :{0}".format(command))
    # write the branching ratio info in it 
    if (channel=="FH") :
        add_br_note = open(output_card_name+".txt", 'a')
        add_br_note.write("CMS_wwgg_br_HH_WWgg      rateParam  *  gghh*  0.000970198 \nCMS_wwgg_br_WW_4Q     rateParam  *  gghh*  0.4489 \nnuisance  edit  freeze  CMS_wwgg_br_HH_WWgg \nnuisance  edit  freeze  CMS_wwgg_br_WW_4Q")
    elif(channel=="SL"):
        add_br_note = open(output_card_name+".txt", 'a')
        add_br_note.write("CMS_wwgg_br_HH_WWgg      rateParam  *  gghh*  0.000970198 \n CMS_wwgg_br_WW_2Qlnu     rateParam  *  gghh*  0.441 \n nuisance  edit  freeze  CMS_wwgg_br_HH_WWgg \n nuisance  edit  freeze  CMS_wwgg_br_WW_2Qlnu")
    elif(channel=="FHSL"):
        add_br_note = open(output_card_name+".txt", 'a')
        add_br_note.write("CMS_wwgg_br_HH_WWgg      rateParam  *  gghh*  0.000970198 \n CMS_wwgg_br_WW_4Q_2Qlnu     rateParam  *  gghh*  0.8899 \n nuisance  edit  freeze  CMS_wwgg_br_HH_WWgg \n nuisance  edit  freeze  CMS_wwgg_br_WW_4Q_2Qlnu")
    logging.info("end: {}".format(time.time()))

def run_combine(output_card_name,log_name,output_file_name):
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Datacard")
    # run Runcombine
    command = "combine -M AsymptoticLimits -m 125 -n " + output_file_name+" "+ output_card_name + ".txt" + " --run expected " + ">" + log_name+  " 2>&1 "
    run_combine_p = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    logging.info("run combine \n command :{0}".format(command))
    logging.info("end: {}".format(time.time()))
def run_combine_card():
    logging.info("begin: {}".format(time.time()))
    os.chdir("/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/flashggFinalFit/Datacard")
    # combine the datacard for same masspoint with categories
    command = "combineCards.py" + "FHSL_1jets_M3000_cat0=Datacard_M3000_1jets_cat0_FHSL.txt FHSL_1jets_M3000_cat1=Datacard_M3000_1jets_cat1_FHSL.txt FHSL_1jets_M3000_cat2=Datacard_M3000_1jets_cat2_FHSL.txt FHSL_1jets_M3000_cat3=Datacard_M3000_1jets_cat3_FHSL.txt   > Datacard_combined_FHSL_1jets_M3000.txt"
    run_combine_card = subprocess.call(command, shell=True, stdout=subprocess.PIPE)
if __name__ == "__main__":
    run_Tree2WS_data(inputpath_name = "/eos/user/z/zhenxuan/wwyy/scan_limit/", ws_data_path="ws_1jets_"+'_for_scan_limit', output_data_root_name="Data_FHSL_2017_cat_1jets_"+'_for_scan_limit'+ "_"+ 'tmp' + ".root")
    run_backgroundfit(ws_data_path="ws_1jets_"+'_for_scan_limit', log_name = "bkg_ws_1jets_"+'_for_scan_limit'+".log" , inputpath_name = "/eos/user/z/zhenxuan/hhwwgg_root/hhwwgg_root_FHSL_custom/" + 'tmp' + "/", ext_name="FHSL_ws_1jets_"+'_for_scan_limit', cp_name="CMS-HGG_multipdf_RECO_untagged_1jets_"+'_for_scan_limit'+"_2017.root")
    run_Tree2WS_sig(inputpath_name = "/eos/user/z/zhenxuan/wwyy/scan_limit/" + 'tmp' + "/",inputfile_name="Signal_"+'tmp'+"_FHSL_2017_1jets_"+"_for_scan_limit"+".root",log_file_name= 'tmp'+"_hhwwgg_MC_FHSL_1jets_"+"_for_scan_limit"+".log", ws_path = 'tmp'+"_1jets_"+"_for_scan_limit", output_sig_root_name="output_Signal"+'tmp'+"_1jets_"+"_for_scan_limit"+"_M125_FHSL_2017_13TeV_amcatnloFXFX_pythia8_gghh.root")
    run_ftest(ws_path = 'tmp' + "_1jets_"+"_for_scan_limit", log_name = "signal_ftest_" + 'tmp' + "_1jets_"+"_for_scan_limit"+"_FHSL.log", inputpath_name= "/eos/user/z/zhenxuan/wwyy/scan_limit/" + 'tmp' + "/")
    run_signalfit(ws_path = 'tmp' + "_1jets_"+"_for_scan_limit", log_name = "signal_signalfit_" + 'tmp' + "_1jets_"+"_for_scan_limit"+"_FHSL.log", inputpath_name="/eos/user/z/zhenxuan/wwyy/scan_limit/" + 'tmp' + "/")
    run_signal_plot(cats="RECO_untagged_1jets_"+"_for_scan_limit", exts="dcb_2017_" + 'tmp' + "_1jets_"+"_for_scan_limit", outputExt="packaged_" + 'tmp' + "_1jets_"+"_for_scan_limit", log_packaged_name = "packaged_" + 'tmp' + "_1jets_"+"_for_scan_limit"+".log", ws_path='tmp' + "_1jets_"+"_for_scan_limit", inputpath_name ="/eos/user/z/zhenxuan/wwyy/scan_limit/" + 'tmp' + "/", log_plotter_name = "plotter_" + 'tmp' + "_1jets_"+"_for_scan_limit"+".log", cp_name="CMS-HGG_sigfit_packaged_RECO_untagged_1jets_"+"_for_scan_limit"+"_2017.root")
    run_yields(ws_sig_path='tmp' + "_1jets_"+"_for_scan_limit", log_name = 'tmp' + "_FH_1jets_"+"_for_scan_limit"+"_yields.log" , inputpath_name = "/eos/user/z/zhenxuan/wwyy/scan_limit/" + 'tmp' + "/", ws_bkg_path = "ws_1jets_"+"_for_scan_limit" )
    run_makeDatacard(ws_sig_path='tmp' + "_1jets_"+"_for_scan_limit", log_name = 'tmp' + "_FH_1jets_"+"_for_scan_limit"+"_makeDatacard.log" , output_card_name = "Datacard_" + 'tmp' + "_1jets_"+"_for_scan_limit"+"_FHSL", channel="FHSL" )
    run_combine(output_file_name='tmp' + "_1jets_"+"_for_scan_limit", log_name = 'tmp' + "_FHSL_1jets_"+"_for_scan_limit"+"_combine_limit.log" , output_card_name = "Datacard_" + 'tmp' + "_1jets_"+"_for_scan_limit"+"_FHSL")

    