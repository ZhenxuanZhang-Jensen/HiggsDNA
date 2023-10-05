#!/bin/python

import ROOT as r
from array import array
import numpy # alternative to array

path="/eos/user/z/zhenxuan/wwyy/root_files/" 
procs = ["Data_cat1_forDataDriven"]
component = ["cat1"]

select = ["min(LeadPhoton_mvaID,SubleadPhoton_mvaID)>=-0.95 && min(LeadPhoton_mvaID,SubleadPhoton_mvaID)<= -0.9 && max(LeadPhoton_mvaID,SubleadPhoton_mvaID)>=-0.9"]

i=0
for proc in procs:
   file=proc+".root"
   f=r.TFile(path+file)
   inputtree="cat1"
   t=f.Get(inputtree) 
   if not t: 
       print( 'Could not find tree '+inputtree+'. Skipping...' )
       i=i+1 
       continue 
   print( '   JTao : input Nevt = ',t.GetEntries()  )
   newfile = "New_"+component[i]+".root"
   nf = r.TFile(newfile,"recreate");
   nt = t.CopyTree(select[i]);
   nt.SetName(component[i])
   print( '   JTao : selected Nevt = ',nt.GetEntries()  )

   
#    r.gDirectory.Delete(ipt+';*') 
   nt.Write();
   
   nf.Write() 
   nf.Close() 
   del t 
   del f
   del nt
   del nf
   print( 'Renamed trees for '+file+" with new file "+newfile)
   i=i+1
