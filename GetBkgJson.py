import sys
import os
import json
input=sys.argv[1]

xsMap={"DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa":84.4,
       'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8':872.1,
       'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8':232.8,
       'GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8':3164.0,
       'QCD_Pt-40ToInf_DoubleEMEnriched_MGG-80ToInf_TuneCP5_13TeV-pythia8':118100.0,
       'GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8':48.58,
       'ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8':0.5071,
       'VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8':2.257,
       'VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8':3.782,
    #    'QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8':24810.0,
       'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8':831.76,
       'TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8':3.055,
       'TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8':4.078,
    #    'TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin-pythia8':0.01687
       }
Data={}
fd = open(input)
outname=input.split(".")[0]
for line in fd:
    sampleName=line.split("/")[1]
    campaign=line.split("/")[2]
    if("UL16" in campaign):
        if("APV" in campaign):
            y="2016preVFP"
        else:
            y="2016postVFP"
    elif("UL17" in campaign):
        y="2017"
    elif("UL18" in campaign):
        y="2018"
    if not Data.has_key(sampleName) :
        # print(Data.has_key(sampleName))
        Data[sampleName]={}
        Data[sampleName]["xs"]=xsMap[sampleName]
        Data[sampleName]["br"]=1.0
        Data[sampleName]["files"]={}
    print(line)
    print(y)
    text=os.popen('dasgoclient -query="file dataset=/%s/%s/NANOAODSIM"'%(sampleName,campaign)).read()
    files=text.replace('/store/mc','root://cms-xrd-global.cern.ch///store/mc')
    OneList=files.split("\n") 
    OneList = filter(None, OneList)               
    # print(OneList)
    Data[sampleName]["files"][y]=OneList
with open(outname+'.json', 'w') as fp:
    json.dump(Data, fp)




    # print('dasgoclient -query="file dataset=/%s/%s/NANOAODSIM"'%(sampleName,campaign))
    # text=os.popen('dasgoclient -query="file dataset=/%s/%s/NANOAODSIM"'%(sampleName,campaign)).read()
    # print(type(text))
    # text=text.replace('.root','.root",')
    # text=text.replace('.root",\n  ','.root"')
    # files=text.replace('/store/mc','            "root://cms-xrd-global.cern.ch///store/mc')
    # print(files)
    # printTail()
    # print(",")
