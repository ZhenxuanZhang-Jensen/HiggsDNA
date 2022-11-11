import sys
import os
import json
input=sys.argv[1]


Data={}
fd = open(input)
out=input.split(".")[0]
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
        print(Data.has_key(sampleName))
        Data[sampleName]={}
        Data[sampleName]["xs"]=0.001
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
with open(out+'.json', 'w') as fp:
    json.dump(Data, fp)
# for line in fd:
    # sampleName=line.split("/")[1]
    # campaign=line.split("/")[2]
    # print sampleName
    # xs=xsMap[sampleName]
    # year=""
    # if ("20UL16" in campaign and "APV" not in campaign):
    #     year="2016"
    #     # printHead(xs,sampleName,year)
    # elif("20UL16" in campaign and "APV" in campaign):
    #     year="2016APV"
    #     # printHead(xs,sampleName,year)
    # elif("20UL17" in campaign):
    #     year="2017"
    #     # printHead(xs,sampleName,year)
    # elif("20UL18" in campaign):
    #     year="2018"



    # print('dasgoclient -query="file dataset=/%s/%s/NANOAODSIM"'%(sampleName,campaign))
    # text=os.popen('dasgoclient -query="file dataset=/%s/%s/NANOAODSIM"'%(sampleName,campaign)).read()
    # print(type(text))
    # text=text.replace('.root','.root",')
    # text=text.replace('.root",\n  ','.root"')
    # files=text.replace('/store/mc','            "root://cms-xrd-global.cern.ch///store/mc')
    # print(files)
    # printTail()
    # print(",")
