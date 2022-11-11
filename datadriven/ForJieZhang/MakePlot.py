import ROOT
import argparse
import os
from array import array

# python Significance_plots.py -d /eos/user/a/atishelm/www/HHWWgg/DNN/HHWWyyDNN_MultiClass_EvenSingleH_2Hgg_withKinWeight_HggClassScale_4_BkgClassScale_1_BalanceYields/Categorization/ --smooth 1 --weight 1
# python Significance_plots.py -d /eos/user/a/atishelm/www/HHWWgg/DNN/HHWWyyDNN_WithHggFactor2-200Epochs-3ClassMulticlass_EvenSingleH_2Hgg_withKinWeightCut10_BalanceYields/Categorization/ --smooth 1 --weight 1
DEBUG = False
def getExtremes(graphs):
  x_vec=[]
  y_vec=[]
  for graph in graphs:
    for point in range(0,graph.GetN()):
      x = ROOT.double(0.)
      y = ROOT.double(0.)
      graph.GetPoint(point,x,y)
      x_vec.append(x)
      y_vec.append(y)
      #print "getExtremes:",x,y
  #print "Final extremes:",[min(x_vec),min(y_vec),max(x_vec),max(y_vec)]
  return [min(x_vec),min(y_vec),max(x_vec),max(y_vec)]

def makeGraph(points):
  i=0
  graph = ROOT.TGraph()
  for key in points:
    graph.SetPoint(i,key[0],key[1])
    i+=1
  return graph

def makeGraph_vs_nCat(nBins, points_1Cat, points_2Cat, points_3Cat, points_4Cat, points_5Cat, points_6Cat, points_7Cat, points_8Cat, points_9Cat):

 graph = ROOT.TGraph()
 for val in points_1Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(0,1,val[1])
 for val in points_2Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(1,2,val[1])
 for val in points_3Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(2,3,val[1])
 for val in points_4Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(3,4,val[1])
 for val in points_5Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(4,5,val[1])
 for val in points_6Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(5,6,val[1])
 for val in points_7Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(6,7,val[1])
 for val in points_8Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(7,8,val[1])
 for val in points_9Cat:
   if int(val[0]) == nBins:
     graph.SetPoint(8,9,val[1])
 return graph

def drawGraphs(graphs,name,xTitle,yTitle,logX,massMin_pos,massMax_pos,outDir):

   legend = ROOT.TLegend(0.15,0.80,0.50,0.94)
   legend.SetFillColor(0)
   legend.SetFillStyle(1000)
   legend.SetTextFont(42)
   legend.SetTextSize(0.025)
   legend.SetNColumns(2)

   colors = (ROOT.kBlack, ROOT.kGray, ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kSpring, ROOT.kTeal, ROOT.kAzure, ROOT.kViolet, ROOT.kPink)
   for i,g in enumerate(graphs):
      g.SetLineWidth(2)
      g.SetLineColor(colors[i])
      g.SetMarkerColor(colors[i])
      g.SetMarkerStyle(20)
      legend.AddEntry(g,"["+str(int(massMin_pos[i]))+","+str(int(massMax_pos[i]))+"] GeV","PL");

   ranges = getExtremes(graphs)
   for i,g in enumerate(graphs):
     g.SetMinimum(ranges[1]*0.95)
     g.SetMaximum(ranges[3]*1.15)
   print("ranges:",ranges)

   c = ROOT.TCanvas()
   c.SetGrid()
   graphs[0].GetXaxis().SetTitle(xTitle)
   graphs[0].GetYaxis().SetTitle(yTitle)
   graphs[0].Draw("APL")
   graphs[0].GetXaxis().SetRangeUser(ranges[0]*0.95,ranges[2]*1.05)
   graphs[0].GetYaxis().SetRangeUser(ranges[1]*0.95,ranges[3]*1.15)
   if logX: c.SetLogx()
   for i,g in enumerate(graphs):
     g.GetXaxis().SetTitle(xTitle)
     g.GetYaxis().SetTitle(yTitle)
     g.Draw("PL,same")
   legend.Draw("same")
  #  outDir
  #  trainingLabel = "HHWWyyDNN_WithHggFactor2-200Epochs-3ClassMulticlass_EvenSingleH_2Hgg_withKinWeightCut10_BalanceYields"
  #  OL = "/eos/user/a/atishelm/www/HHWWgg/DNN/%s/"%(trainingLabel)
  #  OL =
   c.SaveAs(outDir + "/" + name+".png","png")
   c.SaveAs(outDir + "/" + name+".pdf","pdf")
   #c.SaveAs(outDir + "/" + name+".root","root")
   c.SaveAs(outDir + "/" + name+".C","C")
  #  c.SaveAs(name+".png","png")
  #  c.SaveAs(name+".pdf","pdf")

   c2 = ROOT.TCanvas()
   c2.SetGrid()
   graphs[0].GetXaxis().SetTitle(xTitle)
   graphs[0].GetYaxis().SetTitle(yTitle)
   graphs[0].Draw("APL")
   graphs[0].GetXaxis().SetRangeUser(ranges[0]*0.95,ranges[2]*1.05)
   graphs[0].GetYaxis().SetRangeUser(ranges[1]*0.95,ranges[3]*1.15)
   graphs[0].SetMinimum(ranges[1]*0.95)
   graphs[0].SetMaximum(ranges[3]*1.15)
   c2.SetLogy()
   if logX: c2.SetLogx()
   for i,g in enumerate(graphs):
     g.GetXaxis().SetTitle(xTitle)
     g.GetYaxis().SetTitle(yTitle)
     g.Draw("PL,same")
   legend.Draw("same")
   c2.SaveAs(outDir + "/" +name+"_log.png","png")
   c2.SaveAs(outDir + "/" +name+"_log.pdf","pdf")


if __name__ == '__main__':

 ROOT.gROOT.SetBatch(ROOT.kTRUE)

 parser =  argparse.ArgumentParser(description='plot significance')
 parser.add_argument('-d', '--inDir',   dest='inDir',   required=True, type=str)
 parser.add_argument('-s', '--smooth',  dest='smooth',  required=True, type=int)
 parser.add_argument('-w', '--weight',  dest='weight',  required=True, type=int)

 args = parser.parse_args()
 inDir = args.inDir
 useSmoothing =  args.smooth
 useWeight = args.weight


 #inDir = '/eos/user/b/bmarzocc/HHWWgg/January_2021_Production/HHWWyyDNN_binary_withHgg_noNegWeights_BalanceYields_allBkgs_LOSignals_noPtOverM/'

 os.system('ls '+inDir+'/categorize_nBins_*_nCat_*_massMin*.txt > file_dump.txt')
 with open('file_dump.txt') as f_List:
   data_List = f_List.read()
 lines_List = data_List.splitlines()

 nCats_pos = []
 nBins_pos = []
 massMin_pos = []
 massMax_pos = []

 points_1Cat_Total = []
 points_2Cat_Total = []
 points_3Cat_Total = []
 points_4Cat_Total = []
 points_5Cat_Total = []
 points_6Cat_Total = []
 points_7Cat_Total = []
 points_8Cat_Total = []
 points_9Cat_Total = []
 points_1Cat_Total_final = []
 points_2Cat_Total_final = []
 points_3Cat_Total_final = []
 points_4Cat_Total_final = []
 points_5Cat_Total_final = []
 points_6Cat_Total_final = []
 points_7Cat_Total_final = []
 points_8Cat_Total_final = []
 points_9Cat_Total_final = []

 points_1Cat_nSig_Total = []
 points_2Cat_nSig_Total = []
 points_3Cat_nSig_Total = []
 points_4Cat_nSig_Total = []
 points_5Cat_nSig_Total = []
 points_6Cat_nSig_Total = []
 points_7Cat_nSig_Total = []
 points_8Cat_nSig_Total = []
 points_9Cat_nSig_Total = []
 points_1Cat_nSig_Total_final = []
 points_2Cat_nSig_Total_final = []
 points_3Cat_nSig_Total_final = []
 points_4Cat_nSig_Total_final = []
 points_5Cat_nSig_Total_final = []
 points_6Cat_nSig_Total_final = []
 points_7Cat_nSig_Total_final = []
 points_8Cat_nSig_Total_final = []
 points_9Cat_nSig_Total_final = []

 points_1Cat_nBkg_Total = []
 points_2Cat_nBkg_Total = []
 points_3Cat_nBkg_Total = []
 points_4Cat_nBkg_Total = []
 points_5Cat_nBkg_Total = []
 points_6Cat_nBkg_Total = []
 points_7Cat_nBkg_Total = []
 points_8Cat_nBkg_Total = []
 points_9Cat_nBkg_Total = []
 points_1Cat_nBkg_Total_final = []
 points_2Cat_nBkg_Total_final = []
 points_3Cat_nBkg_Total_final = []
 points_4Cat_nBkg_Total_final = []
 points_5Cat_nBkg_Total_final = []
 points_6Cat_nBkg_Total_final = []
 points_7Cat_nBkg_Total_final = []
 points_8Cat_nBkg_Total_final = []
 points_9Cat_nBkg_Total_final = []

 points_1Cat_dataSB_Total = []
 points_2Cat_dataSB_Total = []
 points_3Cat_dataSB_Total = []
 points_4Cat_dataSB_Total = []
 points_5Cat_dataSB_Total = []
 points_6Cat_dataSB_Total = []
 points_7Cat_dataSB_Total = []
 points_8Cat_dataSB_Total = []
 points_9Cat_dataSB_Total = []
 points_1Cat_dataSB_Total_final = []
 points_2Cat_dataSB_Total_final = []
 points_3Cat_dataSB_Total_final = []
 points_4Cat_dataSB_Total_final = []
 points_5Cat_dataSB_Total_final = []
 points_6Cat_dataSB_Total_final = []
 points_7Cat_dataSB_Total_final = []
 points_8Cat_dataSB_Total_final = []
 points_9Cat_dataSB_Total_final = []

 points_1Cat_boundary_Total = []
 points_2Cat_boundary_Total = []
 points_3Cat_boundary_Total = []
 points_4Cat_boundary_Total = []
 points_5Cat_boundary_Total = []
 points_6Cat_boundary_Total = []
 points_7Cat_boundary_Total = []
 points_8Cat_boundary_Total = []
 points_9Cat_boundary_Total = []
 points_1Cat_boundary_Total_final = []
 points_2Cat_boundary_Total_final = []
 points_3Cat_boundary_Total_final = []
 points_4Cat_boundary_Total_final = []
 points_5Cat_boundary_Total_final = []
 points_6Cat_boundary_Total_final = []
 points_7Cat_boundary_Total_final = []
 points_8Cat_boundary_Total_final = []
 points_9Cat_boundary_Total_final = []

 for i,line in enumerate(lines_List):
   line_split = line.split('_')

   #print i,line_split
   nBins = -1
   nCats = -1
   massMin = -1
   massMax = -1

   if 'noSmooth_noReweight' in line:
      nBins = line_split[-8]
      nCats = line_split[-6]
      massMin = line_split[-5]
      massMin = massMin.replace('massMin','')
      massMax = line_split[-4]
      massMax = massMax.replace('massMax','')
   elif 'noSmooth' in line:
      nBins = line_split[-7]
      nCats = line_split[-5]
      massMin = line_split[-4]
      massMin = massMin.replace('massMin','')
      massMax = line_split[-3]
      massMax = massMax.replace('massMax','')
   elif 'noReweight' in line:
      nBins = line_split[-7]
      nCats = line_split[-5]
      massMin = line_split[-4]
      massMin = massMin.replace('massMin','')
      massMax = line_split[-3]
      massMax = massMax.replace('massMax','')
   else:
      nBins = line_split[-6]
      nCats = line_split[-4]
      massMin = line_split[-3]
      massMin = massMin.replace('massMin','')
      massMax = line_split[-2]
      massMax = massMax.replace('massMax','')

   if useSmoothing==1 and 'noSmooth' in line: continue
   if useWeight==1 and 'noReweight' in line: continue

   if int(nCats) not in nCats_pos: nCats_pos.append(int(nCats))
   if int(nBins) not in nBins_pos: nBins_pos.append(int(nBins))
   if float(massMin) not in massMin_pos: massMin_pos.append(float(massMin))
   if float(massMax) not in massMax_pos: massMax_pos.append(float(massMax))

 for i,mMin in enumerate(massMin_pos):

  points_1Cat = {}
  points_2Cat = {}
  points_3Cat = {}
  points_4Cat = {}
  points_5Cat = {}
  points_6Cat = {}
  points_7Cat = {}
  points_8Cat = {}
  points_9Cat = {}

  points_1Cat_nSig = {}
  points_2Cat_nSig = {}
  points_3Cat_nSig = {}
  points_4Cat_nSig = {}
  points_5Cat_nSig = {}
  points_6Cat_nSig = {}
  points_7Cat_nSig = {}
  points_8Cat_nSig = {}
  points_9Cat_nSig = {}

  points_1Cat_nBkg = {}
  points_2Cat_nBkg = {}
  points_3Cat_nBkg = {}
  points_4Cat_nBkg = {}
  points_5Cat_nBkg = {}
  points_6Cat_nBkg = {}
  points_7Cat_nBkg = {}
  points_8Cat_nBkg = {}
  points_9Cat_nBkg = {}

  points_1Cat_dataSB = {}
  points_2Cat_dataSB = {}
  points_3Cat_dataSB = {}
  points_4Cat_dataSB = {}
  points_5Cat_dataSB = {}
  points_6Cat_dataSB = {}
  points_7Cat_dataSB = {}
  points_8Cat_dataSB = {}
  points_9Cat_dataSB = {}

  points_1Cat_boundary = {}
  points_2Cat_boundary = {}
  points_3Cat_boundary = {}
  points_4Cat_boundary = {}
  points_5Cat_boundary = {}
  points_6Cat_boundary = {}
  points_7Cat_boundary = {}
  points_8Cat_boundary = {}
  points_9Cat_boundary = {}

  for nBins in nBins_pos:
   for nCats in nCats_pos:

    file = 0
    if useSmoothing==1 and useWeight==1:
     file = inDir+'/categorize_nBins_'+str(nBins)+'_nCat_'+str(nCats)+'_massMin'+str(massMin_pos[i])+'_massMax'+str(massMax_pos[i])+'_v2.txt'
    elif useSmoothing!=1 and useWeight==1:
     file = inDir+'/categorize_nBins_'+str(nBins)+'_nCat_'+str(nCats)+'_massMin'+str(massMin_pos[i])+'_massMax'+str(massMax_pos[i])+'_v2_noSmooth.txt'
    elif useSmoothing==1 and useWeight!=1:
     file = inDir+'/categorize_nBins_'+str(nBins)+'_nCat_'+str(nCats)+'_massMin'+str(massMin_pos[i])+'_massMax'+str(massMax_pos[i])+'_v2_noReweight.txt'
    elif useSmoothing!=1 and useWeight!=1:
     file = inDir+'/categorize_nBins_'+str(nBins)+'_nCat_'+str(nCats)+'_massMin'+str(massMin_pos[i])+'_massMax'+str(massMax_pos[i])+'_v2_noSmooth_noReweight.txt'
    if not os.path.exists(file): continue
    with open(inDir+'/categorize_nBins_'+str(nBins)+'_nCat_'+str(nCats)+'_massMin'+str(massMin_pos[i])+'_massMax'+str(massMax_pos[i])+'_v2.txt') as f_List:
      data_List = f_List.read()
    lines_List = data_List.splitlines()
    significance = 0.
    nSig = 0.
    nBkg = 0.
    nDataSB = 0.
    boundary = -1.
    if DEBUG: print "\n","="*51
    if DEBUG: print "nCats: ",nCats,"\n"
    for j,cat in enumerate(lines_List):
     if DEBUG: print("cat:",j,cat)
     if 'Tot_Significance:' in cat: significance = cat.split()[-1]
     if len(cat.split())>1 and float(cat.split()[1]) == 1.:
       nSig = float(cat.split()[5])
       nBkg = float(cat.split()[7])
       nDataSB = float(cat.split()[-1])
       boundary = float(cat.split()[0])
       if DEBUG: print('nSig: ',nSig)
       if DEBUG: print('nBkg: ',nBkg)
       if DEBUG: print('nDataSB: ',nDataSB)
       if DEBUG: print('boundary: ',boundary)
    if int(nCats)==1:
     points_1Cat[int(nBins)] = float(significance)
     points_1Cat_nSig[int(nBins)] = float(nSig)
     points_1Cat_nBkg[int(nBins)] = float(nBkg)
     points_1Cat_dataSB[int(nBins)] = float(nDataSB)
     points_1Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==2:
     points_2Cat[int(nBins)] = float(significance)
     points_2Cat_nSig[int(nBins)] = float(nSig)
     points_2Cat_nBkg[int(nBins)] = float(nBkg)
     points_2Cat_dataSB[int(nBins)] = float(nDataSB)
     points_2Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==3:
     points_3Cat[int(nBins)] = float(significance)
     points_3Cat_nSig[int(nBins)] = float(nSig)
     points_3Cat_nBkg[int(nBins)] = float(nBkg)
     points_3Cat_dataSB[int(nBins)] = float(nDataSB)
     points_3Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==4:
     points_4Cat[int(nBins)] = float(significance)
     points_4Cat_nSig[int(nBins)] = float(nSig)
     points_4Cat_nBkg[int(nBins)] = float(nBkg)
     points_4Cat_dataSB[int(nBins)] = float(nDataSB)
     points_4Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==5:
     points_5Cat[int(nBins)] = float(significance)
     points_5Cat_nSig[int(nBins)] = float(nSig)
     points_5Cat_nBkg[int(nBins)] = float(nBkg)
     points_5Cat_dataSB[int(nBins)] = float(nDataSB)
     points_5Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==6:
     points_6Cat[int(nBins)] = float(significance)
     points_6Cat_nSig[int(nBins)] = float(nSig)
     points_6Cat_nBkg[int(nBins)] = float(nBkg)
     points_6Cat_dataSB[int(nBins)] = float(nDataSB)
     points_6Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==7:
     points_7Cat[int(nBins)] = float(significance)
     points_7Cat_nSig[int(nBins)] = float(nSig)
     points_7Cat_nBkg[int(nBins)] = float(nBkg)
     points_7Cat_dataSB[int(nBins)] = float(nDataSB)
     points_7Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==8:
     points_8Cat[int(nBins)] = float(significance)
     points_8Cat_nSig[int(nBins)] = float(nSig)
     points_8Cat_nBkg[int(nBins)] = float(nBkg)
     points_8Cat_dataSB[int(nBins)] = float(nDataSB)
     points_8Cat_boundary[int(nBins)] = float(boundary)
    elif int(nCats)==9:
     points_9Cat[int(nBins)] = float(significance)
     points_9Cat_nSig[int(nBins)] = float(nSig)
     points_9Cat_nBkg[int(nBins)] = float(nBkg)
     points_9Cat_dataSB[int(nBins)] = float(nDataSB)
     points_9Cat_boundary[int(nBins)] = float(boundary)

  points_1Cat = sorted(points_1Cat.items())
  points_2Cat = sorted(points_2Cat.items())
  points_3Cat = sorted(points_3Cat.items())
  points_4Cat = sorted(points_4Cat.items())
  points_5Cat = sorted(points_5Cat.items())
  points_6Cat = sorted(points_6Cat.items())
  points_7Cat = sorted(points_7Cat.items())
  points_8Cat = sorted(points_8Cat.items())
  points_9Cat = sorted(points_9Cat.items())
  points_1Cat_nSig = sorted(points_1Cat_nSig.items())
  points_2Cat_nSig = sorted(points_2Cat_nSig.items())
  points_3Cat_nSig = sorted(points_3Cat_nSig.items())
  points_4Cat_nSig = sorted(points_4Cat_nSig.items())
  points_5Cat_nSig = sorted(points_5Cat_nSig.items())
  points_6Cat_nSig = sorted(points_6Cat_nSig.items())
  points_7Cat_nSig = sorted(points_7Cat_nSig.items())
  points_8Cat_nSig = sorted(points_8Cat_nSig.items())
  points_9Cat_nSig = sorted(points_9Cat_nSig.items())
  points_1Cat_nBkg = sorted(points_1Cat_nBkg.items())
  points_2Cat_nBkg = sorted(points_2Cat_nBkg.items())
  points_3Cat_nBkg = sorted(points_3Cat_nBkg.items())
  points_4Cat_nBkg = sorted(points_4Cat_nBkg.items())
  points_5Cat_nBkg = sorted(points_5Cat_nBkg.items())
  points_6Cat_nBkg = sorted(points_6Cat_nBkg.items())
  points_7Cat_nBkg = sorted(points_7Cat_nBkg.items())
  points_8Cat_nBkg = sorted(points_8Cat_nBkg.items())
  points_9Cat_nBkg = sorted(points_9Cat_nBkg.items())
  points_1Cat_dataSB = sorted(points_1Cat_dataSB.items())
  points_2Cat_dataSB = sorted(points_2Cat_dataSB.items())
  points_3Cat_dataSB = sorted(points_3Cat_dataSB.items())
  points_4Cat_dataSB = sorted(points_4Cat_dataSB.items())
  points_5Cat_dataSB = sorted(points_5Cat_dataSB.items())
  points_6Cat_dataSB = sorted(points_6Cat_dataSB.items())
  points_7Cat_dataSB = sorted(points_7Cat_dataSB.items())
  points_8Cat_dataSB = sorted(points_8Cat_dataSB.items())
  points_9Cat_dataSB = sorted(points_9Cat_dataSB.items())
  points_1Cat_boundary = sorted(points_1Cat_boundary.items())
  points_2Cat_boundary = sorted(points_2Cat_boundary.items())
  points_3Cat_boundary = sorted(points_3Cat_boundary.items())
  points_4Cat_boundary = sorted(points_4Cat_boundary.items())
  points_5Cat_boundary = sorted(points_5Cat_boundary.items())
  points_6Cat_boundary = sorted(points_6Cat_boundary.items())
  points_7Cat_boundary = sorted(points_7Cat_boundary.items())
  points_8Cat_boundary = sorted(points_8Cat_boundary.items())
  points_9Cat_boundary = sorted(points_9Cat_boundary.items())

  points_1Cat_Total.append(points_1Cat)
  points_2Cat_Total.append(points_2Cat)
  points_3Cat_Total.append(points_3Cat)
  points_4Cat_Total.append(points_4Cat)
  points_5Cat_Total.append(points_5Cat)
  points_6Cat_Total.append(points_6Cat)
  points_7Cat_Total.append(points_7Cat)
  points_8Cat_Total.append(points_8Cat)
  points_9Cat_Total.append(points_9Cat)
  points_1Cat_Total_final.append(makeGraph(points_1Cat))
  points_2Cat_Total_final.append(makeGraph(points_2Cat))
  points_3Cat_Total_final.append(makeGraph(points_3Cat))
  points_4Cat_Total_final.append(makeGraph(points_4Cat))
  points_5Cat_Total_final.append(makeGraph(points_5Cat))
  points_6Cat_Total_final.append(makeGraph(points_6Cat))
  points_7Cat_Total_final.append(makeGraph(points_7Cat))
  points_8Cat_Total_final.append(makeGraph(points_8Cat))
  points_9Cat_Total_final.append(makeGraph(points_9Cat))

  points_1Cat_nSig_Total.append(points_1Cat_nSig)
  points_2Cat_nSig_Total.append(points_2Cat_nSig)
  points_3Cat_nSig_Total.append(points_3Cat_nSig)
  points_4Cat_nSig_Total.append(points_4Cat_nSig)
  points_5Cat_nSig_Total.append(points_5Cat_nSig)
  points_6Cat_nSig_Total.append(points_6Cat_nSig)
  points_7Cat_nSig_Total.append(points_7Cat_nSig)
  points_8Cat_nSig_Total.append(points_8Cat_nSig)
  points_9Cat_nSig_Total.append(points_9Cat_nSig)
  points_1Cat_nSig_Total_final.append(makeGraph(points_1Cat_nSig))
  points_2Cat_nSig_Total_final.append(makeGraph(points_2Cat_nSig))
  points_3Cat_nSig_Total_final.append(makeGraph(points_3Cat_nSig))
  points_4Cat_nSig_Total_final.append(makeGraph(points_4Cat_nSig))
  points_5Cat_nSig_Total_final.append(makeGraph(points_5Cat_nSig))
  points_6Cat_nSig_Total_final.append(makeGraph(points_6Cat_nSig))
  points_7Cat_nSig_Total_final.append(makeGraph(points_7Cat_nSig))
  points_8Cat_nSig_Total_final.append(makeGraph(points_8Cat_nSig))
  points_9Cat_nSig_Total_final.append(makeGraph(points_9Cat_nSig))

  points_1Cat_nBkg_Total.append(points_1Cat_nBkg)
  points_2Cat_nBkg_Total.append(points_2Cat_nBkg)
  points_3Cat_nBkg_Total.append(points_3Cat_nBkg)
  points_4Cat_nBkg_Total.append(points_4Cat_nBkg)
  points_5Cat_nBkg_Total.append(points_5Cat_nBkg)
  points_6Cat_nBkg_Total.append(points_6Cat_nBkg)
  points_7Cat_nBkg_Total.append(points_7Cat_nBkg)
  points_8Cat_nBkg_Total.append(points_8Cat_nBkg)
  points_9Cat_nBkg_Total.append(points_9Cat_nBkg)
  points_1Cat_nBkg_Total_final.append(makeGraph(points_1Cat_nBkg))
  points_2Cat_nBkg_Total_final.append(makeGraph(points_2Cat_nBkg))
  points_3Cat_nBkg_Total_final.append(makeGraph(points_3Cat_nBkg))
  points_4Cat_nBkg_Total_final.append(makeGraph(points_4Cat_nBkg))
  points_5Cat_nBkg_Total_final.append(makeGraph(points_5Cat_nBkg))
  points_6Cat_nBkg_Total_final.append(makeGraph(points_6Cat_nBkg))
  points_7Cat_nBkg_Total_final.append(makeGraph(points_7Cat_nBkg))
  points_8Cat_nBkg_Total_final.append(makeGraph(points_8Cat_nBkg))
  points_9Cat_nBkg_Total_final.append(makeGraph(points_9Cat_nBkg))

  points_1Cat_dataSB_Total.append(points_1Cat_dataSB)
  points_2Cat_dataSB_Total.append(points_2Cat_dataSB)
  points_3Cat_dataSB_Total.append(points_3Cat_dataSB)
  points_4Cat_dataSB_Total.append(points_4Cat_dataSB)
  points_5Cat_dataSB_Total.append(points_5Cat_dataSB)
  points_6Cat_dataSB_Total.append(points_6Cat_dataSB)
  points_7Cat_dataSB_Total.append(points_7Cat_dataSB)
  points_8Cat_dataSB_Total.append(points_8Cat_dataSB)
  points_9Cat_dataSB_Total.append(points_9Cat_dataSB)
  points_1Cat_dataSB_Total_final.append(makeGraph(points_1Cat_dataSB))
  points_2Cat_dataSB_Total_final.append(makeGraph(points_2Cat_dataSB))
  points_3Cat_dataSB_Total_final.append(makeGraph(points_3Cat_dataSB))
  points_4Cat_dataSB_Total_final.append(makeGraph(points_4Cat_dataSB))
  points_5Cat_dataSB_Total_final.append(makeGraph(points_5Cat_dataSB))
  points_6Cat_dataSB_Total_final.append(makeGraph(points_6Cat_dataSB))
  points_7Cat_dataSB_Total_final.append(makeGraph(points_7Cat_dataSB))
  points_8Cat_dataSB_Total_final.append(makeGraph(points_8Cat_dataSB))
  points_9Cat_dataSB_Total_final.append(makeGraph(points_9Cat_dataSB))

  points_1Cat_boundary_Total.append(points_1Cat_boundary)
  points_2Cat_boundary_Total.append(points_2Cat_boundary)
  points_3Cat_boundary_Total.append(points_3Cat_boundary)
  points_4Cat_boundary_Total.append(points_4Cat_boundary)
  points_5Cat_boundary_Total.append(points_5Cat_boundary)
  points_6Cat_boundary_Total.append(points_6Cat_boundary)
  points_7Cat_boundary_Total.append(points_7Cat_boundary)
  points_8Cat_boundary_Total.append(points_8Cat_boundary)
  points_9Cat_boundary_Total.append(points_9Cat_boundary)
  points_1Cat_boundary_Total_final.append(makeGraph(points_1Cat_boundary))
  points_2Cat_boundary_Total_final.append(makeGraph(points_2Cat_boundary))
  points_3Cat_boundary_Total_final.append(makeGraph(points_3Cat_boundary))
  points_4Cat_boundary_Total_final.append(makeGraph(points_4Cat_boundary))
  points_5Cat_boundary_Total_final.append(makeGraph(points_5Cat_boundary))
  points_6Cat_boundary_Total_final.append(makeGraph(points_6Cat_boundary))
  points_7Cat_boundary_Total_final.append(makeGraph(points_7Cat_boundary))
  points_8Cat_boundary_Total_final.append(makeGraph(points_8Cat_boundary))
  points_9Cat_boundary_Total_final.append(makeGraph(points_9Cat_boundary))

 vec_bins = []
 for val in points_1Cat_Total[0]:
   vec_bins.append(int(val[0]))

 for nBins in vec_bins:

   points_Significance = []
   points_nSig = []
   points_nBkg = []
   points_dataSB = []
   points_boundary = []

   for i in range(0,len(points_1Cat_Total)):
     points_Significance.append(makeGraph_vs_nCat(nBins, points_1Cat_Total[i], points_2Cat_Total[i], points_3Cat_Total[i], points_4Cat_Total[i], points_5Cat_Total[i], points_6Cat_Total[i], points_7Cat_Total[i], points_8Cat_Total[i], points_9Cat_Total[i]))
     points_nSig.append(makeGraph_vs_nCat(nBins, points_1Cat_nSig_Total[i], points_2Cat_nSig_Total[i], points_3Cat_nSig_Total[i], points_4Cat_nSig_Total[i], points_5Cat_nSig_Total[i], points_6Cat_nSig_Total[i], points_7Cat_nSig_Total[i], points_8Cat_nSig_Total[i], points_9Cat_nSig_Total[i]))
     points_nBkg.append(makeGraph_vs_nCat(nBins, points_1Cat_nBkg_Total[i], points_2Cat_nBkg_Total[i], points_3Cat_nBkg_Total[i], points_4Cat_nBkg_Total[i], points_5Cat_nBkg_Total[i], points_6Cat_nBkg_Total[i], points_7Cat_nBkg_Total[i], points_8Cat_nBkg_Total[i], points_9Cat_nBkg_Total[i]))
     points_dataSB.append(makeGraph_vs_nCat(nBins, points_1Cat_dataSB_Total[i], points_2Cat_dataSB_Total[i], points_3Cat_dataSB_Total[i], points_4Cat_dataSB_Total[i], points_5Cat_dataSB_Total[i], points_6Cat_dataSB_Total[i], points_7Cat_dataSB_Total[i], points_8Cat_dataSB_Total[i], points_9Cat_dataSB_Total[i]))
     points_boundary.append(makeGraph_vs_nCat(nBins, points_1Cat_boundary_Total[i], points_2Cat_boundary_Total[i], points_3Cat_boundary_Total[i], points_4Cat_boundary_Total[i], points_5Cat_boundary_Total[i], points_6Cat_boundary_Total[i], points_7Cat_boundary_Total[i], points_8Cat_boundary_Total[i], points_9Cat_boundary_Total[i]))

   # print "length:",len(points_Significance)
   # print("minMass: %i, maxMass: %i, bin: %i"%(massMin_pos,massMax_pos,nBins))
   if useSmoothing==1 and useWeight==1:
     drawGraphs(points_Significance,"Significance_vs_nCats_nBins_"+str(nBins),"nCats","Significance",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nSig,"nSig_vs_nCats_nBins_"+str(nBins),"nCats","nSig",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nBkg,"nBkg_vs_nCats_nBins_"+str(nBins),"nCats","nBkg",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_dataSB,"dataSB_vs_nCats_nBins_"+str(nBins),"nCats","dataSB",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_boundary,"boundary_vs_nCats_nBins_"+str(nBins),"nCats","boundary",False,massMin_pos,massMax_pos, inDir)
   elif useSmoothing!=1 and useWeight==1:
     drawGraphs(points_Significance,"Significance_vs_nCats_nBins_"+str(nBins)+"_noSmooth","nCats","Significance",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nSig,"nSig_vs_nCats_nBins_"+str(nBins)+"_noSmooth","nCats","nSig",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nBkg,"nBkg_vs_nCats_nBins_"+str(nBins)+"_noSmooth","nCats","nBkg",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_dataSB,"dataSB_vs_nCats_nBins_"+str(nBins)+"_noSmooth","nCats","dataSB",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_boundary,"boundary_vs_nCats_nBins_"+str(nBins)+"_noSmooth","nCats","boundary",False,massMin_pos,massMax_pos, inDir)
   elif useSmoothing==1 and useWeight!=1:
     drawGraphs(points_Significance,"Significance_vs_nCats_nBins_"+str(nBins)+"_noReweight","nCats","Significance",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nSig,"nSig_vs_nCats_nBins_"+str(nBins)+"_noReweight","nCats","nSig",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nBkg,"nBkg_vs_nCats_nBins_"+str(nBins)+"_noReweight","nCats","nBkg",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_dataSB,"dataSB_vs_nCats_nBins_"+str(nBins)+"_noReweight","nCats","dataSB",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_boundary,"boundary_vs_nCats_nBins_"+str(nBins)+"_noReweight","nCats","boundary",False,massMin_pos,massMax_pos, inDir)
   elif useSmoothing!=1 and useWeight!=1:
     drawGraphs(points_Significance,"Significance_vs_nCats_nBins_"+str(nBins)+"_noSmooth_noReweight","nCats","Significance",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nSig,"nSig_vs_nCats_nBins_"+str(nBins)+"_noSmooth_noReweight","nCats","nSig",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_nBkg,"nBkg_vs_nCats_nBins_"+str(nBins)+"_noSmooth_noReweight","nCats","nBkg",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_dataSB,"dataSB_vs_nCats_nBins_"+str(nBins)+"_noSmooth_noReweight","nCats","dataSB",False,massMin_pos,massMax_pos, inDir)
     drawGraphs(points_boundary,"boundary_vs_nCats_nBins_"+str(nBins)+"_noSmooth_noReweight","nCats","boundary",False,massMin_pos,massMax_pos, inDir)

 if useSmoothing==1 and useWeight==1:
   drawGraphs(points_1Cat_Total_final,"Significance_vs_nBins_nCats_1","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nSig_Total_final,"nSig_vs_nBins_nCats_1","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_1","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_1","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_boundary_Total_final,"boundary_vs_nBins_nCats_1","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_2Cat_Total_final,"Significance_vs_nBins_nCats_2","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nSig_Total_final,"nSig_vs_nBins_nCats_2","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_2","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_2","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_boundary_Total_final,"boundary_vs_nBins_nCats_2","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_3Cat_Total_final,"Significance_vs_nBins_nCats_3","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nSig_Total_final,"nSig_vs_nBins_nCats_3","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_3","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_3","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_boundary_Total_final,"boundary_vs_nBins_nCats_3","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_4Cat_Total_final,"Significance_vs_nBins_nCats_4","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nSig_Total_final,"nSig_vs_nBins_nCats_4","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_4","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_4","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_boundary_Total_final,"boundary_vs_nBins_nCats_4","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_5Cat_Total_final,"Significance_vs_nBins_nCats_5","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nSig_Total_final,"nSig_vs_nBins_nCats_5","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_5","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_5","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_boundary_Total_final,"boundary_vs_nBins_nCats_5","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   # drawGraphs(points_6Cat_Total_final,"Significance_vs_nBins_nCats_6","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_6Cat_nSig_Total_final,"nSig_vs_nBins_nCats_6","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_6Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_6","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_6Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_6","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_6Cat_boundary_Total_final,"boundary_vs_nBins_nCats_6","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   # drawGraphs(points_7Cat_Total_final,"Significance_vs_nBins_nCats_7","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_7Cat_nSig_Total_final,"nSig_vs_nBins_nCats_7","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_7Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_7","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_7Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_7","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_7Cat_boundary_Total_final,"boundary_vs_nBins_nCats_7","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   # drawGraphs(points_8Cat_Total_final,"Significance_vs_nBins_nCats_8","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_8Cat_nSig_Total_final,"nSig_vs_nBins_nCats_8","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_8Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_8","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_8Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_8","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_8Cat_boundary_Total_final,"boundary_vs_nBins_nCats_8","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   # drawGraphs(points_9Cat_Total_final,"Significance_vs_nBins_nCats_9","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_9Cat_nSig_Total_final,"nSig_vs_nBins_nCats_9","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_9Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_9","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_9Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_9","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   # drawGraphs(points_9Cat_boundary_Total_final,"boundary_vs_nBins_nCats_9","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

 elif useSmoothing!=1 and useWeight==1:
   drawGraphs(points_1Cat_Total_final,"Significance_vs_nBins_nCats_1_noSmooth","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nSig_Total_final,"nSig_vs_nBins_nCats_1_noSmooth","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_1_noSmooth","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_1_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_boundary_Total_final,"boundary_vs_nBins_nCats_1_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_2Cat_Total_final,"Significance_vs_nBins_nCats_2_noSmooth","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nSig_Total_final,"nSig_vs_nBins_nCats_2_noSmooth","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_2_noSmooth","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_2_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_boundary_Total_final,"boundary_vs_nBins_nCats_2_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_3Cat_Total_final,"Significance_vs_nBins_nCats_3_noSmooth","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nSig_Total_final,"nSig_vs_nBins_nCats_3_noSmooth","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_3_noSmooth","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_3_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_boundary_Total_final,"boundary_vs_nBins_nCats_3_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_4Cat_Total_final,"Significance_vs_nBins_nCats_4_noSmooth","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nSig_Total_final,"nSig_vs_nBins_nCats_4_noSmooth","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_4_noSmooth","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_4_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_boundary_Total_final,"boundary_vs_nBins_nCats_4_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_5Cat_Total_final,"Significance_vs_nBins_nCats_5_noSmooth","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nSig_Total_final,"nSig_vs_nBins_nCats_5_noSmooth","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_5_noSmooth","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_5_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_boundary_Total_final,"boundary_vs_nBins_nCats_5_noSmooth","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

 elif useSmoothing==1 and useWeight!=1:
   drawGraphs(points_1Cat_Total_final,"Significance_vs_nBins_nCats_1_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nSig_Total_final,"nSig_vs_nBins_nCats_1_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_1_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_1_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_boundary_Total_final,"boundary_vs_nBins_nCats_1_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_2Cat_Total_final,"Significance_vs_nBins_nCats_2_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nSig_Total_final,"nSig_vs_nBins_nCats_2_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_2_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_2_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_boundary_Total_final,"boundary_vs_nBins_nCats_2_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_3Cat_Total_final,"Significance_vs_nBins_nCats_3_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nSig_Total_final,"nSig_vs_nBins_nCats_3_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_3_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_3_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_boundary_Total_final,"boundary_vs_nBins_nCats_3_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_4Cat_Total_final,"Significance_vs_nBins_nCats_4","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nSig_Total_final,"nSig_vs_nBins_nCats_4","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_4","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_4","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_boundary_Total_final,"boundary_vs_nBins_nCats_4_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_5Cat_Total_final,"Significance_vs_nBins_nCats_5_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nSig_Total_final,"nSig_vs_nBins_nCats_5_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_5_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_5_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_boundary_Total_final,"boundary_vs_nBins_nCats_5_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

 elif useSmoothing!=1 and useWeight!=1:
   drawGraphs(points_1Cat_Total_final,"Significance_vs_nBins_nCats_1_noSmooth_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nSig_Total_final,"nSig_vs_nBins_nCats_1_noSmooth_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_1_noSmooth_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_1_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_1Cat_boundary_Total_final,"boundary_vs_nBins_nCats_1_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_2Cat_Total_final,"Significance_vs_nBins_nCats_2_noSmooth_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nSig_Total_final,"nSig_vs_nBins_nCats_2_noSmooth_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_2_noSmooth_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_2_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_2Cat_boundary_Total_final,"boundary_vs_nBins_nCats_2_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_3Cat_Total_final,"Significance_vs_nBins_nCats_3_noSmooth_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nSig_Total_final,"nSig_vs_nBins_nCats_3_noSmooth_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_3_noSmooth_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_3_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_3Cat_boundary_Total_final,"boundary_vs_nBins_nCats_3_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_4Cat_Total_final,"Significance_vs_nBins_nCats_4_noSmooth_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nSig_Total_final,"nSig_vs_nBins_nCats_4_noSmooth_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_4_noSmooth_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_4_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_4Cat_boundary_Total_final,"boundary_vs_nBins_nCats_4_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

   drawGraphs(points_5Cat_Total_final,"Significance_vs_nBins_nCats_5_noSmooth_noReweight","nBins","Significance",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nSig_Total_final,"nSig_vs_nBins_nCats_5_noSmooth_noReweight","nBins","nSig",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_nBkg_Total_final,"nBkg_vs_nBins_nCats_5_noSmooth_noReweight","nBins","nBkg",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_dataSB_Total_final,"dataSB_vs_nBins_nCats_5_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)
   drawGraphs(points_5Cat_boundary_Total_final,"boundary_vs_nBins_nCats_5_noSmooth_noReweight","nBins","dataSB",True,massMin_pos,massMax_pos, inDir)

 os.system('rm file_dump.txt')