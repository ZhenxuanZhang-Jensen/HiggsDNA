import sys
import os
import json
input=sys.argv[1]


xsMap={"Data":None,
       }
# years=["20UL16NanoAODAPVv9","20UL16NanoAODv9","20UL17NanoAODv9","20UL18NanoAODv9"]
years=["2016","2017","2018"]
Data={}
Data["Data"]={}
Data["Data"]["files"]={}
for year in years:
    fd = open(input)
    YearList=[]
    for line in fd:
        if( "Run"+year in line):
            sampleName=line.split("/")[1]
            campaign=line.split("/")[2]
            # print(sampleName,campaign)
            text=os.popen('dasgoclient -query="file dataset=/%s/%s/NANOAOD"'%(sampleName,campaign)).read()
            files=text.replace('/store/data','root://cms-xrd-global.cern.ch///store/data')
            OneList=files.split("\n") 
            OneList = filter(None, OneList)               
            print(OneList)
            YearList.extend(OneList)
    Data["Data"]["files"][year]=YearList
Data["Data"]["fpo"]=10
with open('Data_UL.json', 'w') as fp:
    json.dump(Data, fp)

        # OneList.append(files)
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
