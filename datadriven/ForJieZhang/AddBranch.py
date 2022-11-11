#!/usr/bin/python
import numpy as n
from ROOT import TH1F, TFile, TCanvas, THStack
from ROOT import gStyle, TTree
import ROOT
from ROOT import gROOT, TChain, TH1F, kBlue, kRed, TLegend, TCanvas, TLorentzVector, TList
from array import array
import operator
import math
import sys
import os
import argparse
import random
from math import *







fname=sys.argv[1]
scale=sys.argv[2]
name=fname.split(".")[0]
out = TFile.Open(name+"_reweight.root", "RECREATE")
out.Close()
tfile = ROOT.TFile(fname)
# if ("Data" in filen):
    # tfile.cd("tagsDumper/trees")
# tfile.cd("tagsDumper/trees")
for key in ROOT.gDirectory.GetListOfKeys():
    tname = key.GetName()
    print ("tree name is :", tname)
    tree = tfile.Get(tname)
    weight = array('f', [0])
    weight[0] = -999
    dZ = array('f', [0])
    dZ[0] = -999
    out = TFile.Open(name+"_reweight.root", "UPDATE")

    # newtree = ch_0.CopyTree(
        #  "((N_goodElectrons>1 && N_goodMuons==0) || (N_goodElectrons==0 && N_goodMuons>1) || (N_goodElectrons>0 && N_goodMuons>0))")
    newtree = tree.CopyTree("(category==1)")
    tree.SetBranchStatus("*", 1)
    _weight = newtree.Branch(
        "weight", weight, "weight/F")
    _dZ = newtree.Branch(
        "dZ", dZ, "dZ/F")
    print (newtree.GetEntries())
    #  for i in range(0,ch_0.GetEntries()):
    #for i in range(0, newtree.GetEntries()):
    for i in range(0, newtree.GetEntries()):
        newtree.GetEntry(i)
        # print(float(newtree.weight_central))
        if("Data" in fname):
            weight[0]=(newtree.weight_central)*float(scale)
        else:
            weight[0]=(newtree.weight_central_no_lumi)*float(scale)
        dZ[0]=float(2.0)
        if(i==20):
            print(weight[0])
        #print FL_Leading_lepton_pt[0],FL_Leading_lepton_eta[0],FL_Leading_lepton_phi[0],FL_Leading_lepton_E[0],FL_Subleading_lepton_pt,FL_Subleading_lepton_eta,FL_Subleading_lepton_phi,FL_Subleading_lepton_E
        #print newtree.goodElectrons_0_pt,newtree.goodElectrons_1_pt,newtree.goodElectrons_3_pt,newtree.goodElectrons_4_pt
        #print newtree.goodMuons_0_pt,newtree.goodMuons_1_pt,newtree.goodMuons_3_pt,newtree.goodMuons_4_pt
        _weight.Fill()
        _dZ.Fill()
        #print "Before:",newtree.MET_phi,"After:",METCor_phi[0]
        #print Flavor,"DiLepMass",DiLeptonMass[0],"DR:",FLDR[0]
    newtree.Write()
    print (newtree.GetEntries())
    out.Close()


