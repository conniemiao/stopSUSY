#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, cutMode, findingSameFlavor, muPreference, plotVar
# Draws 1D hist for data for some variable, for the summed bkgd data and for each 
# of the signal files.
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from sig_SingleStop_files

import sys
from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

assert len(sys.argv) == 6, "need 5 command line args: testMode{0,1}, cutMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}, plotVar"
plotVar = sys.argv[5]

plotSettings = { # [nBins,xMin,xMax,listForm]
        "lep1_pt":[100,0,400,False],
        "lep1_eta":[100,-3,3,False],
        "lep1_phi":[100,-4,4,False],
        "lep1_relIso":[10,0,0.2,False],
        "lep2_pt":[100,0,400,False],
        "lep2_eta":[100,-4,4,False],
        "lep2_phi":[100,-4,4,False],
        "lep2_relIso":[10,0,0.2,False],
        "njets":[10,0,10,False],
        "jet_pt":[100,0,400,True], 
        "jet_eta":[100,-3,3,True],
        "jet_phi":[100,-4,4,True],
        "nbtag":[10,0,10,False],
        "deltaR_lep1_jet":[100,0,7,False],
        "deltaR_lep2_jet":[100,0,7,False],
        "mtlep1":[100,0,500,False],
        "mtlep2":[100,0,500,False],
        "met_pt":[100,0,500,False],
        "met_phi":[100,-4,-4,False],
        }

assert (plotVar in plotSettings), "invalid plotVar"

# Determining adr of bkgd and sig ntuples.
# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
# applying cuts
cutMode = bool(int(sys.argv[2]))
print "Cut mode:", cutMode
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
findingSameFlavor = bool(int(sys.argv[3]))
print "Finding same flavor:", findingSameFlavor
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
muPreference = bool(int(sys.argv[4]))
print "Mu preference:", muPreference
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

print "Plotting",plotVar,"from",bkgdNtupleAdr,"and",sigsNtupleAdr

numSigFiles = int(sigsNtupleAdr[48:50])
nBins = plotSettings[plotVar][0]
# if not testMode and nBins > 20: nBins = nBins * 5
xMin = plotSettings[plotVar][1]
xMax = plotSettings[plotVar][2]
listForm = plotSettings[plotVar][3] # only for some of the jet variables 

binwidth = (xMax - xMin)/nBins # include overflow bin
hBkgd = TH1D(plotVar + "_bkgd", plotVar + "_bkgd", nBins + 1, \
        xMin, xMax + binwidth)
lumi = 3000000 # luminosity = 3000 /fb = 3,000,000 /pb

c1 = TCanvas("c1","Plot",10,20,1000,700)
gStyle.SetOptStat(0) # don't show any stats

#--------------------------------------------------------------------------------#
# *************** Filling bkgd data summed together  ************
print "Plotting", plotVar, "from background."
xsec = 67.75 # production cross section
bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
tBkgd = bkgdFile.Get("tBkgd")

nentries = tBkgd.GetEntries()
print("nentries={0:d}".format(nentries))
assert nentries > 0, "You have no events in your tree..."

for count, event in enumerate(tBkgd):
    if count % 500000 == 0: print("count={0:d}".format(count))
    val = getattr(event, plotVar)
    if listForm:
        numjets = event.njets
        val = np.reshape(val, 20)
        for j in range(numjets):
            jetval = val[j]
            if jetval <= xMax:
                hBkgd.Fill(jetval, 1)
            else: # overflow
                hBkgd.Fill(xMax + binwidth/2, 1)
    else:
        if val <= xMax:
            hBkgd.Fill(val, 1)
        else: # overflow
            hBkgd.Fill(xMax + binwidth/2, 1)
hBkgd.Sumw2()

# rebinning
# print hBkgd.GetBinError((xMax+xMin)/2)/hBkgd.GetSumOfWeights()
# while hBkgd.GetBinError((xMax+xMin)/2) > 0.1*hBkgd.GetSumOfWeights():
#     hBkgd.Rebin(2)
#     print rebinned

title = plotVar + " ("+bkgdNtupleAdr[-18:-14]+", normalized to 3000 /fb)"
title += ", "+bkgdNtupleAdr[-13:-5]
hBkgd.SetTitle(title)
hBkgd.GetXaxis().SetTitle(plotVar+" [GeV]")
hBkgd.GetYaxis().SetTitle("Number of Events")
hBkgd.Scale(xsec*lumi/hBkgd.GetSumOfWeights())
hBkgd.SetMinimum(1)
hBkgd.SetMaximum(10**12)
hBkgd.SetLineColor(1) # black
# hBkgd.SetStats(0)
hBkgd.Draw("hist")
c1.Update()
bkgdFile.Close()

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVar + " from signal."
sigDataListFile = open("sig_SingleStop_files")

coloropts = [2,4,3,6,7,9,28,46] # some good colors for lines
markeropts = [1,20,21,22,23] # some good marker styles for lines
linestyleopts = [1,2,3,7,9] # some good styles for lines

hSigArr = [] # one hist for each signal file 
for fileNum, line in enumerate(sigDataListFile):
    if fileNum + 1 > numSigFiles: break

    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename

    sigFile = TFile.Open(sigsNtupleAdr, "READ")
    tSig = sigFile.Get("tSig"+str(fileNum))
    nentries = tSig.GetEntries()
    print("nentries={0:d}".format(nentries))
    assert nentries > 0, "You have no events in your tree..."

    hSigArr.append(TH1D(plotVar + "_sig_" + filename, plotVar + "_sig_" + \
            filename[19:24], nBins + 1, xMin, xMax + binwidth))
    
    for count, event in enumerate(tSig):
        if count % 500000 == 0: print("count={0:d}".format(count))
        val = getattr(event, plotVar)
        if listForm:
            numjets = event.njets
            val = np.reshape(val, 20)
            for j in range(numjets):
                jetval = val[j]
                if jetval <= xMax:
                    hSigArr[fileNum].Fill(jetval, 1)
                else: # overflow
                    hSigArr[fileNum].Fill(xMax + binwidth/2, 1)
        else:
            if val <= xMax:
                hSigArr[fileNum].Fill(val, 1)
            else: # overflow
                hSigArr[fileNum].Fill(xMax + binwidth/2, 1)

    hSigArr[fileNum].SetDirectory(0) # necessary to keep hist from closing

    hcolor = coloropts[fileNum % len(coloropts)]
    hSigArr[fileNum].SetLineColor(hcolor) 
    hmarkerstyle = markeropts[(fileNum/len(coloropts)) % len(markeropts)]
    hSigArr[fileNum].SetMarkerStyle(hmarkerstyle)
    hSigArr[fileNum].SetMarkerColor(hcolor)
    hlinestyle = linestyleopts[(fileNum/len(coloropts)/len(markeropts)) % \
            len(linestyleopts)]
    hSigArr[fileNum].SetLineStyle(hlinestyle)

    hSigArr[fileNum].Sumw2()
    hSigArr[fileNum].Scale(xsec * lumi / hSigArr[fileNum].GetSumOfWeights())
    hSigArr[fileNum].SetMinimum(1)
    hSigArr[fileNum].SetMaximum(10**12)
    hSigArr[fileNum].Draw("hist same") # same pad
    c1.Update()
    sigFile.Close()


#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************

legend = TLegend(.70,.75,.90,.90)
legend.AddEntry(hBkgd, plotVar + "_bkgd")
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
            "plots/v2CutSequence/"+plotVar+"_"+l1Flav[:2]+l2Flav[:2]+"_"+\
            sigsNtupleAdr[-13:-5]+".png"
    print "Saving image", imgName
    gSystem.ProcessEvents()
    img = TImage.Create()
    img.FromPad(c1)
    img.WriteImage(imgName)
    print "Done."

