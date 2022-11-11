import ROOT
import sys
fname=sys.argv[1]
isSignal=sys.argv[2] #1 signal 2 singleH 3 Data
ext=sys.argv[3]
# print(isSignal)
# cuts=[0.1,0.79,0.91,0.94,1]
cuts=[0.6,0.85,0.91,0.94,0.96,1]
file=ROOT.TFile.Open(fname)
tree=file.Get("parquettree")
# if(isSignal=="1"):
	# morecuts="&& (CMS_hgg_mass>=0)"
# else:
    # morecuts="&& (CMS_hgg_mass>=0)"
outfile=ROOT.TFile.Open(fname.split(".")[0]+"_"+ext+"_categorised.root","RECREATE")
if(isSignal=="WWgg"):
    newtree1=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fath4qjet_deepTagMD_H4qvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s) "%(0.7,0.6))
if(isSignal=="1"):
    newtree1=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s) "%(cuts[0],cuts[1]))
    newtree2=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[1],cuts[2]))
    newtree3=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[2],cuts[3]))
    newtree4=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[3],cuts[4]))
    newtree5=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[4],cuts[5]))
else:
    newtree1=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s) "%(cuts[0],cuts[1]))
    newtree2=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[1],cuts[2]))
    newtree3=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[2],cuts[3]))
    newtree4=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[3],cuts[4]))
    newtree5=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fathbbjet_deepTagMD_HbbvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s ) "%(cuts[4],cuts[5]))
    newtree6=tree.CopyTree("(Leading_Photon_MVA>-0.7 && Subleading_Photon_MVA>-0.7) && (fath4qjet_deepTagMD_H4qvsQCD>%s && fathbbjet_deepTagMD_HbbvsQCD<=%s) "%(0.7,0.6))

if(isSignal=="1"):
	newtree1.SetName("GluGluToHHTo2B2G_"+ext+"_13TeV_HHTag_BBgg_0")
	newtree2.SetName("GluGluToHHTo2B2G_"+ext+"_13TeV_HHTag_BBgg_1")
	newtree3.SetName("GluGluToHHTo2B2G_"+ext+"_13TeV_HHTag_BBgg_2")
	newtree4.SetName("GluGluToHHTo2B2G_"+ext+"_13TeV_HHTag_BBgg_3")
	newtree5.SetName("GluGluToHHTo2B2G_"+ext+"_13TeV_HHTag_BBgg_4")
if(isSignal=="WWgg"):
	newtree1.SetName("GluGluToHHTo2G4Q_"+ext+"_13TeV_HHTag_WWgg_0")
elif(isSignal=="0"):
    newtree1.SetName("Data_13TeV_HHTag_BBgg_0")
    newtree2.SetName("Data_13TeV_HHTag_BBgg_1")
    newtree3.SetName("Data_13TeV_HHTag_BBgg_2")
    newtree4.SetName("Data_13TeV_HHTag_BBgg_3")
    newtree5.SetName("Data_13TeV_HHTag_BBgg_4")
    newtree6.SetName("Data_13TeV_HHTag_WWgg_0")
elif(isSignal=="2"):
    newtree1.SetName(ext+"_125_13TeV_HHTag_BBgg_0")
    newtree2.SetName(ext+"_125_13TeV_HHTag_BBgg_1")
    newtree3.SetName(ext+"_125_13TeV_HHTag_BBgg_2")
    newtree4.SetName(ext+"_125_13TeV_HHTag_BBgg_3")
    newtree5.SetName(ext+"_125_13TeV_HHTag_BBgg_4")
    newtree6.SetName(ext+"_125_13TeV_HHTag_WWgg_0")
if(isSignal=="WWgg"):
    print(newtree1.GetEntries())
    newtree1.Write("",ROOT.TObject.kOverwrite)
elif(isSignal=="1"):
    print(newtree1.GetEntries())
    print(newtree2.GetEntries())
    print(newtree3.GetEntries())
    print(newtree4.GetEntries())
    print(newtree5.GetEntries())
    newtree1.Write("",ROOT.TObject.kOverwrite)
    newtree2.Write("",ROOT.TObject.kOverwrite)
    newtree3.Write("",ROOT.TObject.kOverwrite)
    newtree4.Write("",ROOT.TObject.kOverwrite)
    newtree5.Write("",ROOT.TObject.kOverwrite)
else:
    print(newtree1.GetEntries())
    print(newtree2.GetEntries())
    print(newtree3.GetEntries())
    print(newtree4.GetEntries())
    print(newtree5.GetEntries())
    print(newtree6.GetEntries())
    newtree1.Write("",ROOT.TObject.kOverwrite)
    newtree2.Write("",ROOT.TObject.kOverwrite)
    newtree3.Write("",ROOT.TObject.kOverwrite)
    newtree4.Write("",ROOT.TObject.kOverwrite)
    newtree5.Write("",ROOT.TObject.kOverwrite)
    newtree6.Write("",ROOT.TObject.kOverwrite)
outfile.Close()
# output_M125_GluGluHHTo2B2G_M2000_HHTag_BBgg_0.root


