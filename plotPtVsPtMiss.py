#!/usr/bin/env python

# Draws hist for met pt vs. pt for some input root

from ROOT import TFile, TTree, TH2F, TCanvas, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

testMode = False 
c = TCanvas("c")

bkgdFile = TFile.Open("/eos/user/a/alkaloge/HLLHC/Skims/v3/DESY_pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/TTJets_DiLept_ntuple_20.root", "READ")
bkgdTree = bkgdFile.Get("AC1B")
nentries = bkgdTree.GetEntries()
print("nentries={0:d}".format(nentries))
histBkgd = TH2F("histBkgd", "Pt miss vs. Pt", 100, 0, 400, 100, 0, 400)
nMax = nentries
if testMode: nMax = 5000
for count, entry in enumerate(bkgdTree):
    if count > nMax : break
    if count % 500000 == 0: print("count={0:d}".format(count))

    mIndex = -1
    maxPt = 0
    for im in range(entry.muon_count):
        pt = list(entry.muon_pt)[im]
        if pt > 25 and pt > maxPt and abs(list(entry.muon_eta)[im]) < 2.4: 
            maxPt = pt
            mIndex = im
    if mIndex == -1: continue
    histBkgd.Fill(list(entry.muon_pt)[mIndex], entry.pfmet_pt) 
histBkgd.SetDirectory(0)
histBkgd.Draw()

sigFile = TFile.Open("/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/single-stop14TeV_R_220_A.root", "READ")
sigTree = sigFile.Get("AC1B")
nentries = sigTree.GetEntries()
print("nentries={0:d}".format(nentries))
histSig = TH2F("histSig", "Pt miss vs. Pt", 100, 0, 400, 100, 0, 400)
nMax = nentries
if testMode: nMax = 5000
for count, entry in enumerate(sigTree):
    if count > nMax : break
    if count % 500000 == 0: print("count={0:d}".format(count))

    mIndex = -1
    maxPt = 0
    for im in range(entry.muon_count):
        pt = list(entry.muon_pt)[im]
        if pt > 25 and pt > maxPt and abs(list(entry.muon_eta)[im]) < 2.4: 
            maxPt = pt
            mIndex = im
    if mIndex == -1: continue
    histSig.Fill(list(entry.muon_pt)[mIndex], entry.pfmet_pt) 
histSig.SetDirectory(0)
histSig.SetLineColor(2)
histSig.SetMarkerColor(2)
histSig.Draw("same")

gStyle.SetOptStat(0) # don't show any stats
legend = TLegend(.75,.80,.95,.95)
legend.AddEntry(histBkgd, histBkgd.GetTitle())
legend.AddEntry(histSig, histSig.GetTitle())
legend.Draw("same")

histBkgd.GetXaxis().SetTitle("Pt")
histBkgd.GetYaxis().SetTitle("Pt miss")

c.Update()
print "Done"

raw_input()
