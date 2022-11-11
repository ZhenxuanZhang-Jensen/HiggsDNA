
#!/usr/bin/python
import numpy as n
import ROOT
import os
import sys, getopt
from array import array
from optparse import OptionParser
import operator
import math

def reduceTree(inTree, cut):
  small = inTree.CopyTree(str(cut))
  return small

def computeSignificance(s,b,d_noSmooth):
  significance = -999.
  if b>0 and s>0.: significance = (2*(s+b)*math.log(1+(s/b))) - 2*s
  if significance>0. and d_noSmooth>=8. and b>0.: return math.sqrt((2*(s+b)*math.log(1+(s/b))) - 2*s)
  else: return -999.

def sumSignificance(partition, h_sig_SR, h_bkg_SR,h_data_SB_noSmooth):
  sum = 0.
  for pair in partition:
    s = 2*h_sig_SR.Integral(pair[0],pair[1])
    b = h_bkg_SR.Integral(pair[0],pair[1])
    d_noSmooth = h_data_SB_noSmooth.Integral(pair[0],pair[1])
    significance = computeSignificance(s,b,d_noSmooth)
    #print( h_sig_SR.GetBinCenter(pair[0])-h_bdt_signal_SR.GetBinWidth(pair[0])/2.,significance,b
    if significance>0.: sum += significance*significance
    else: return -999.
  return math.sqrt(sum)

if __name__ == '__main__':

  ROOT.gROOT.SetBatch(ROOT.kTRUE)

  parser = OptionParser()
  parser.add_option( "-d", "--inDir",   dest="inDir",    default="",   type="string", help="inDir"  )
  parser.add_option( "-n", "--nBins",   dest="nBins",    default=190,  type="int",    help="nBins"  )
  parser.add_option( "-c", "--nCats",   dest="nCats",    default=5,    type="int",    help="nCats"  )
  parser.add_option( "-m", "--min",     dest="min",    default=5,    type="string",    help="minCut"  )
  parser.add_option( "-r", "--massMin", dest='massMin',  default=120,  type="float",  help="massMin")
  parser.add_option( "-R", "--massMax", dest='massMax',  default=130,  type="float",  help="massMax")
  parser.add_option( "-s", "--smooth",  dest='smooth',   default=1,    type="int",    help="smooth")
  parser.add_option( "-w", "--weight",  dest='weight',   default=1,    type="int",    help="weight")
  parser.add_option('-o', "--outDir", dest='outDir', default="", type="string", help="outDir")
  (options, args) = parser.parse_args()

  nBins = options.nBins
  inDir = options.inDir
  nCats = options.nCats
  massMin = options.massMin
  min=options.min
  massMax = options.massMax
  useSmoothing =  options.smooth
  useWeight = options.weight
  # options.outDir
  #inDir = '/eos/user/b/bmarzocc/HHWWgg/January_2021_Production/HHWWyyDNN_binary_withHgg_noNegWeights_BalanceYields_allBkgs_LOSignals_noPtOverM/'

  print( "inDir       :",inDir)
  print( "nBins       :",nBins)
  print( "nCats       :",nCats)
  print( "massMin     :",massMin)
  print( "massMax     :",massMax)
  print( "UseSmoothing:",useSmoothing)
  print( "useWeight   :",useWeight)

  inFile = ROOT.TFile(inDir+"/DNN_Histos_smoothing_SmoothSuper_bins"+str(nBins)+"_massMin"+str(massMin)+"_massMax"+str(massMax)+'_minCut'+min+".root","READ")
  print("/DNN_Histos_smoothing_SmoothSuper_bins"+str(nBins)+"_massMin"+str(massMin)+"_massMax"+str(massMax)+'_minCut'+min+".root")

  h_bdt_datamix_SR_weighted_smooth = 0

  h_bdt_signal_SR = inFile.Get("h_DNN_signal_SR")
  h_bdt_datamix_SR_weighted = inFile.Get("h_DNN_bkg_SR_weighted")
  h_bdt_datamix_SR_weighted_smooth = inFile.Get("h_DNN_bkg_SR_weighted_smooth")
  if useSmoothing!=1 and useWeight==1: h_bdt_datamix_SR_weighted_smooth = inFile.Get("h_DNN_bkg_SR_weighted")
  if useSmoothing==1 and useWeight!=1: h_bdt_datamix_SR_weighted_smooth = inFile.Get("h_DNN_bkg_SR_smooth")
  if useSmoothing!=1 and useWeight!=1: h_bdt_datamix_SR_weighted_smooth = inFile.Get("h_DNN_bkg_SR")
  h_bdt_datamix_SB_weighted = inFile.Get("h_DNN_bkg_SB_weighted")
  h_bdt_datamix_SB_weighted_smooth = inFile.Get("h_DNN_bkg_SB_weighted_smooth")
  h_bdt_data_SB = inFile.Get("h_DNN_data_SB")
  h_bdt_data_SB_smooth = inFile.Get("h_DNN_data_SB_smooth")

  significance_final = -999.
  partition_final = []

  #1 categories
  if nCats == 1:

   for i in range(1,nBins+1):
       partition = [[i,nBins]]
       significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
       #print( h_bdt_signal_SR.GetBinCenter(partition[0][0])-h_bdt_signal_SR.GetBinWidth(partition[0][0])/2,"1. --->",significance,h_bdt_data_SB.Integral(partition[0][0],nBins)
       #print( h_bdt_signal_SR.GetBinCenter(partition[0][0])-h_bdt_signal_SR.GetBinWidth(partition[0][0])/2,"1. --->",significance
       if significance>significance_final:
         significance_final = significance
         partition_final = partition
   print( nCats," - Best category: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2,"1. --->",significance_final)

  #2 categories
  elif nCats == 2:

   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     partition = [[1,i],[j,nBins]]
     if abs(i-j)==1:
       #print( partition
       significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
       if significance>significance_final:
         significance_final = significance
         partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2,h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2,"1. --->",significance_final)

  #3 categories
  elif nCats == 3:

   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      partition = [[1,i],[j,k-1],[k,nBins]]
      if abs(i-j)==1:
        #print( partition
        significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
        if significance>significance_final:
          significance_final = significance
          partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, "1. --->",significance_final)

  #4 categories
  elif nCats == 4:

   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      for d in range(k+1,nBins+1):
       partition = [[1,i],[j,k-1],[k,d-1],[d,nBins]]
       if abs(i-j)==1:
         #print( partition
         significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
         if significance>significance_final:
           significance_final = significance
           partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, "1. --->",significance_final)

  #5 categories
  elif nCats == 5:
   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      for d in range(k+1,nBins+1):
       for f in range(d+1,nBins+1):
        partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,nBins]]
        if abs(i-j)==1:
          #print( partition
          significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
          if significance>significance_final:
            significance_final = significance
            partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final)

  #6 categories
  elif nCats == 6:
   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      for d in range(k+1,nBins+1):
       for f in range(d+1,nBins+1):
         for g in range(f+1,nBins+1):
          partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,g-1],[g,nBins]]
          if abs(i-j)==1:
            #print( partition
            significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
            if significance>significance_final:
              significance_final = significance
              partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final)

  #7 categories
  elif nCats == 7:
   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      for d in range(k+1,nBins+1):
       for f in range(d+1,nBins+1):
         for g in range(f+1,nBins+1):
           for h in range(g+1,nBins+1):
            partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,g-1],[g,h-1],[h,nBins]]
            if abs(i-j)==1:
              #print( partition
              significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
              if significance>significance_final:
                significance_final = significance
                partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final)

  #8 categories
  elif nCats == 8:
   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      for d in range(k+1,nBins+1):
       for f in range(d+1,nBins+1):
         for g in range(f+1,nBins+1):
           for h in range(g+1,nBins+1):
             for l in range(h+1,nBins+1):
              partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,g-1],[g,h-1],[h,l-1],[l,nBins]]
              if abs(i-j)==1:
                #print( partition
                significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                if significance>significance_final:
                  significance_final = significance
                  partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final)

  #8 categories
  elif nCats == 9:
   for i in range(1,nBins+1):
    for j in range(i+1,nBins+1):
     for k in range(j+1,nBins+1):
      for d in range(k+1,nBins+1):
       for f in range(d+1,nBins+1):
         for g in range(f+1,nBins+1):
           for h in range(g+1,nBins+1):
             for l in range(h+1,nBins+1):
               for m in range(l+1,nBins+1):
                partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,g-1],[g,h-1],[h,l-1],[l,m-1],[m,nBins]]
                if abs(i-j)==1:
                  #print( partition
                  significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                  if significance>significance_final:
                    significance_final = significance
                    partition_final = partition

   print( nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final)

  else:
   print( "Number of categories not supported, choose: 1, 2, 3, 4, 5, 6, 7, 8 or 9!")
   sys.exit()

  #final details
  if options.outDir == "":
    outFileName = inDir+"/categorize_nBins_"+str(nBins)+"_nCat_"+str(nCats)+"_massMin"+str(massMin)+"_massMax"+str(massMax)+"_v2.txt"
  else:
    os.system("mkdir %s"%(options.outDir))
    outFileName = options.outDir+"/categorize_nBins_"+str(nBins)+"_nCat_"+str(nCats)+"_massMin"+str(massMin)+"_massMax"+str(massMax)+"_v2.txt"
  if useSmoothing!=1 and useWeight==1: outFileName = outFileName.replace('.txt','_noSmooth.txt')
  if useSmoothing==1 and useWeight!=1: outFileName = outFileName.replace('.txt','_noReweight.txt')
  if useSmoothing!=1 and useWeight!=1: outFileName = outFileName.replace('.txt','_noSmooth_noReweight.txt')
  outFile = open(outFileName,"w+")
  print( "Final details:")
  for pair in partition_final:
    s = 2*h_bdt_signal_SR.Integral(pair[0],pair[1])
    b = h_bdt_datamix_SR_weighted_smooth.Integral(pair[0],pair[1])
    m = h_bdt_datamix_SB_weighted_smooth.Integral(pair[0],pair[1])
    d = h_bdt_data_SB_smooth.Integral(pair[0],pair[1])
    d_noSmooth = h_bdt_data_SB.Integral(pair[0],pair[1])
    significance = computeSignificance(s,b,d_noSmooth)
    print( h_bdt_signal_SR.GetBinCenter(pair[0])-h_bdt_signal_SR.GetBinWidth(pair[0])/2., h_bdt_signal_SR.GetBinCenter(pair[1])+h_bdt_signal_SR.GetBinWidth(pair[1])/2., " --> Significance:", significance, " - N Sig:", s, "- N DataMix_SR:", b, "- N DataMix_SB:", m, "- N Data_SB:", d, "- N Data_SB (no smooth):", d_noSmooth)
    outFile.write(str(h_bdt_signal_SR.GetBinCenter(pair[0])-h_bdt_signal_SR.GetBinWidth(pair[0])/2.))
    outFile.write("  ")
    outFile.write(str(h_bdt_signal_SR.GetBinCenter(pair[1])+h_bdt_signal_SR.GetBinWidth(pair[1])/2.))
    outFile.write("  Significance: %s"%(str(significance)) )
    outFile.write(" NSig: %s"%(s))
    outFile.write(" NDataMix_SR: %s"%(b))
    outFile.write(" NDataMix_SB: %s"%(m))
    outFile.write(" NData_SB: %s"%(d))
    outFile.write(" NData_SB (no smooth): %s"%(d_noSmooth))
    outFile.write("\n")
  outFile.write("\n")
  outFile.write("Tot_Significance: %s"%(significance_final))
  outFile.close()