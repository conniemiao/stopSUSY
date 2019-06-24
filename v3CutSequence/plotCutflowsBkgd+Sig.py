#!/usr/bin/env python

# Overlays the cutflow hists for bkgd and each of the signal files.
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from sig_SingleStop_files

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TImage, TLegend, TText
from ROOT import gSystem, gStyle, gROOT, kTRUE
import numpy as np
from math import sqrt

assert len(sys.argv) == 5, "need 4 command line args: testMode{0,1}, cutMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}"

# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
# applying cuts
cutMode = bool(int(sys.argv[2]))
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
findingSameFlavor = bool(int(sys.argv[3]))
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
muPreference = bool(int(sys.argv[4]))
# copy in the bkgd and sigs filenames from makeNtupleBkgd.py and makeNtupleSigs.py
if findingSameFlavor:
    if muPreference: 
        l1Flav = "muon"
        l2Flav = "muon"
    else: 
        l1Flav = "electron"
        l2Flav = "electron"
else: 
    l1Flav = "muon"
    l2Flav = "electron"

# assemble the sigsNtupleAdr and bkgdNtupleAdr
# number of files to process
numBkgdFiles = 27  # need to loop over all the files in order to have correct xsec
if testMode: 
    numBkgdFiles = 2 
numSigFiles = 2 # max 25
baseDir = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/"
bkgdNtupleAdr = baseDir+"stopCut_"
sigsNtupleAdr = baseDir+"stopCut_"
if numSigFiles < 10: sigsNtupleAdr += "0"+str(numSigFiles)
else: sigsNtupleAdr += str(numSigFiles)
if numBkgdFiles < 10: bkgdNtupleAdr += "0"+str(numBkgdFiles)
else: bkgdNtupleAdr += str(numBkgdFiles)
bkgdNtupleAdr += "Bkgd_TTDiLept_"+l1Flav[:2]+l2Flav[:2]
sigsNtupleAdr += "Sig_"+l1Flav[:2]+l2Flav[:2]
if not cutMode: 
    sigsNtupleAdr += "_baseline.root"
    bkgdNtupleAdr += "_baseline.root"
else:
    sigsNtupleAdr += "_withcuts.root"
    bkgdNtupleAdr += "_withcuts.root"
if not testMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases

# assert bkgdNtupleAdr[50:54] == "Bkgd", "bkgdNtupleAdr not bkgd"
# assert sigsNtupleAdr[50:53] == "Sig", "sigsNtupleAdr not sigs"
# assert bkgdNtupleAdr[-18:] == sigsNtupleAdr[-18:], "sigs/bkgd settings don't match"
print "Plotting from",bkgdNtupleAdr,"and",sigsNtupleAdr

lumi = 3000000 # luminosity = 3000 /fb = 3,000,000 /fb
c1 = TCanvas("c1","Plot",10,20,1000,700)
gStyle.SetOptStat(0) # don't show any stats

xsec = 67.75 # production cross section
nEvtsLabels = []

bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
hBkgd = bkgdFile.Get("bkgd_cutflow")
nCuts = hBkgd.GetNbinsX()
hBkgd.Sumw2()
hBkgd.SetTitle("Cutflow ("+bkgdNtupleAdr[-18:-14]+", norm to 3000 /fb)")
hBkgd.GetYaxis().SetTitle("Number of Events")
hBkgd.Scale(xsec*lumi/hBkgd.GetSumOfWeights())
hBkgd.SetLineColor(1) # black
hBkgd.SetMinimum(1)
hBkgd.SetMaximum(10**12)
hBkgd.Draw("hist")
c1.Update()


sigDataListFile = open("sig_SingleStop_files")
coloropts = [2,4,3,6,7,9,28,46] # some good colors for lines
markeropts = [1,20,21,22,23] # some good marker styles for lines
linestyleopts = [1,2,3,7,9] # some good styles for lines

hSigArr = []
print
for fileNum, line in enumerate(sigDataListFile):
    if fileNum + 1 > numSigFiles: break
    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)

    sigFile = TFile.Open(sigsNtupleAdr, "READ")
    hSigArr.append(sigFile.Get("sig_"+filename[21:31]+"_cutflow"))
    hSigArr[fileNum].SetDirectory(0)

    hcolor = coloropts[fileNum % len(coloropts)]
    hSigArr[fileNum].SetLineColor(hcolor) 
    hmarkerstyle = markeropts[(fileNum/len(coloropts)) % len(markeropts)]
    hSigArr[fileNum].SetMarkerStyle(hmarkerstyle)
    hSigArr[fileNum].SetMarkerColor(hcolor)
    hlinestyle = linestyleopts[(fileNum/len(coloropts)/len(markeropts)) % \
            len(linestyleopts)]
    hSigArr[fileNum].SetLineStyle(hlinestyle)

    hSigArr[fileNum].Sumw2()
    hSigArr[fileNum].Scale(xsec * lumi /
            hSigArr[fileNum].GetSumOfWeights())
    hSigArr[fileNum].SetMinimum(1)
    hSigArr[fileNum].SetMaximum(10**12)
    hSigArr[fileNum].Draw("hist same") # same pad, draw marker
    c1.Update()

    print "Num surviving events after each cut from sig %s:" % filename 
    for i in range(0,nCuts):
        # print hBkgd.GetXaxis().GetBinLabel(i+1),"S/sqrt(B):",\
        #         hSigArr[fileNum].GetBinContent(i+1)/sqrt(hBkgd.GetBinContent(i+1))
        hSig = hSigArr[fileNum]
        print hBkgd.GetXaxis().GetBinLabel(i+1),hSig.GetBinContent(i+1)
        nEvtsLabel = TText()
        nEvtsLabel.SetNDC()
        nEvtsLabel.SetTextSize(0.02)
        nEvtsLabel.SetTextAlign(22)
        nEvtsLabel.SetTextAngle(0)
        nEvtsLabel.SetTextColor(hcolor)
        nEvtsLabel.DrawText(0.1+0.4/nCuts+0.8*float(i)/nCuts, \
                0.1+(1+fileNum)*0.02, str(int(hSig.GetBinContent(i+1))))
        nEvtsLabels.append(nEvtsLabel)
    print

print "Num surviving events after each cut from bkgd:" 
for i in range(0,nCuts):
    print hBkgd.GetXaxis().GetBinLabel(i+1),hBkgd.GetBinContent(i+1)
    nEvtsLabel = TText()
    nEvtsLabel.SetNDC()
    nEvtsLabel.SetTextSize(0.02)
    nEvtsLabel.SetTextAlign(22)
    nEvtsLabel.SetTextAngle(0)
    nEvtsLabel.SetTextColor(1)
    nEvtsLabel.DrawText(0.1+0.4/nCuts+0.8*float(i)/nCuts, \
            0.1+(numSigFiles+1)*0.02, str(int(hBkgd.GetBinContent(i+1))))
    nEvtsLabels.append(nEvtsLabel)
print

legend = TLegend(.70,.75,.90,.90)
legend.AddEntry(hBkgd, "bkgd_cutflow")
legend.SetTextSize(0.02)
for fileNum, h in enumerate(hSigArr):
    legend.AddEntry(hSigArr[fileNum], hSigArr[fileNum].GetTitle())
legend.Draw("same")

c1.SetLogy()
c1.Update()

if testMode:
    print "Done. Press enter to finish."
    raw_input()
else:
    imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
            "plots/v3CutSequence/cutflow_"+ l1Flav[:2]+l2Flav[:2]+"_"+\
            sigsNtupleAdr[-13:-5]+".png"
    print "Saving image", imgName
    gSystem.ProcessEvents()
    img = TImage.Create()
    img.FromPad(c1)
    img.WriteImage(imgName)
    print "Done."
