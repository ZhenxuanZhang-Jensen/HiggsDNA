#!/usr/bin/python
import numpy as n
from ROOT import *
import sys, getopt
from array import array
import itertools
import argparse
import operator
import os



if __name__ == '__main__':

  parser =  argparse.ArgumentParser(description='cat MVA')
  parser.add_argument('-d', '--inDir',  dest='inDir',   required=True,  type=str)
  parser.add_argument('-m', '--min',    dest='min',     required=False, type=float)
  parser.add_argument('-M', '--max',    dest='max',     required=False, type=float)
  parser.add_argument('-r', '--massMin',dest='massMin', required=False, type=float)
  parser.add_argument('-R', '--massMax',dest='massMax', required=False, type=float)
  parser.add_argument('-n', '--nStep',  dest='nStep',   required=False, type=int)
  parser.add_argument("-s", "--suffix", dest='suffix',  required=True,  type=str, help="suffix")
  # parser.add_argument('--outDir', required=False, default="/eos/user/r/rasharma/www/doubleHiggs/Categorization/HHWWyyDNN_binary_May05_NewModel5_E700_LR10em5_B100_ELU_DR0p1_Nadam_DefaultVar_CW_BalanceNonWeightedv2/", type=str, help = "Output directory for everything")
  parser.add_argument('--outDir', required=False, default="/eos/user/c/chuw/FHDNN_forSmooth/output/", type=str, help = "Output directory for everything")

  args = parser.parse_args()
  inDir = args.inDir
  #inDir='/eos/user/b/bmarzocc/HHWWgg/January_2021_Production/HHWWyyDNN_binary_noHgg_noNegWeights_BalanceYields_allBkgs_LOSignals_noPtOverM/'

  min = 0.1
  if args.min: min = args.min
  max = 1.
  if args.max: max = args.max
  massMin = 115.
  if args.massMin: massMin = args.massMin
  massMax = 135.
  if args.massMax: massMax = args.massMax
  nStep = 8
  if args.nStep: nStep = args.nStep

  print "inDir:",inDir
  print "bdtMin:",min
  print "bdtMax:",max
  print "massMin:",massMin
  print "massMin:",massMax
  print "nStep:",nStep

  local = os.getcwd()
  if not os.path.isdir('error'): os.mkdir('error')
  if not os.path.isdir('output'): os.mkdir('output')
  if not os.path.isdir('log'): os.mkdir('log')

  # Prepare condor jobs
  condor = '''executable              = run_script_%s.sh
output                  = output/strips_$(ClusterId)_$(ProcId).out
error                   = output/strips_$(ClusterId)_$(ProcId).out
log                     = log/strips_$(ClusterId)_$(ProcId).log
transfer_input_files    = run_script_%s.sh
on_exit_remove          = (ExitBySignal == False) && (ExitCode == 0)
periodic_release        = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*60))
+JobFlavour             = "espresso"
queue arguments from arguments_1_%s.txt
+JobFlavour             = "microcentury"
queue arguments from arguments_2_%s.txt
+JobFlavour             = "longlunch"
queue arguments from arguments_3_%s.txt
+JobFlavour             = "workday"
queue arguments from arguments_4_%s.txt
'''
# +AccountingGroup        = "group_u_CMS.CAF.ALCA"
# workday
# longlunch
# microcentury
# espresso

  with open("condor_job_%s.jdl"%(args.suffix), "w") as cnd_out:
     cnd_out.write(condor%(args.suffix,args.suffix,args.suffix,args.suffix,args.suffix,args.suffix))

  outputDir = os.getcwd()

  script = '''#!/bin/sh -e
echo "Starting job on " `date`
echo "Running on: `uname -a`"
echo "System software: `cat /etc/redhat-release`"
JOBID=$1;
LOCAL=$2;
INPUTDIR=$3;
NBINS=$4;
MIN=$5
MAX=$6
MASSMIN=$7
MASSMAX=$8
echo -e "evaluate"
eval `scramv1 ru -sh`
echo -e "smoothing...";
time(python ${LOCAL}/smooth_DNN.py -d ${INPUTDIR} -n ${NBINS} --min ${MIN} --max ${MAX} --massMin ${MASSMIN} --massMax ${MASSMAX} --outDir %s)
echo -e "DONE";
echo "Ending job on " `date`
'''

  # nBins = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,380,760,1520]
  # nBins1 = [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
  nBins1 = [10, 20, 30, 40]
  # nBins2 = [42,44,46,48,50,60]
  nBins2 = [50, 60]
  nBins3 = [70,80,90,100,110]
  nBins4 = [120,130,140,150,160,170,180,190,380,760]
  # nBins4 = [120,130,140,150,160,170,180,190,380,760,1520]
  #==================================================
  arguments=[]
  for iBin in nBins1:
     for i  in range(0,nStep):
        arguments.append("{jobID} {LOCAL} {indir} {nbins} {min} {max} {massMin} {massMax}".format(jobID = 1,LOCAL = local,indir = inDir+"/",nbins = iBin,min = min,max = max,massMin = massMin+i,massMax = massMax-i))
  with open("arguments_1_%s.txt"%(args.suffix), "w") as argfile:
     argfile.write("\n".join(arguments))
  #==================================================
  arguments=[]
  for iBin in nBins2:
     for i  in range(0,nStep):
        arguments.append("{jobID} {LOCAL} {indir} {nbins} {min} {max} {massMin} {massMax}".format(jobID = 2,LOCAL = local,indir = inDir+"/",nbins = iBin,min = min,max = max,massMin = massMin+i,massMax = massMax-i))
  with open("arguments_2_%s.txt"%(args.suffix), "w") as argfile:
     argfile.write("\n".join(arguments))
  #==================================================
  arguments=[]
  for iBin in nBins3:
     for i  in range(0,nStep):
        arguments.append("{jobID} {LOCAL} {indir} {nbins} {min} {max} {massMin} {massMax}".format(jobID = 3,LOCAL = local,indir = inDir+"/",nbins = iBin,min = min,max = max,massMin = massMin+i,massMax = massMax-i))
  with open("arguments_3_%s.txt"%(args.suffix), "w") as argfile:
     argfile.write("\n".join(arguments))
  #==================================================
  arguments=[]
  for iBin in nBins4:
     for i  in range(0,nStep):
        arguments.append("{jobID} {LOCAL} {indir} {nbins} {min} {max} {massMin} {massMax}".format(jobID = 4,LOCAL = local,indir = inDir+"/",nbins = iBin,min = min,max = max,massMin = massMin+i,massMax = massMax-i))
  with open("arguments_4_%s.txt"%(args.suffix), "w") as argfile:
     argfile.write("\n".join(arguments))
  #==================================================

  with open("run_script_%s.sh"%(args.suffix), "w") as scriptfile:
    scriptfile.write(script%(args.outDir))

  print('For running condor jobs do following:')
  print('1. set up proxy:')
#   print('\tvoms-proxy-init --voms cms --valid 168:00')
  print('2. Submit the condor jobs:')
#   print("\tcondor_submit condor_job_{}.jdl".format(args.suffix))