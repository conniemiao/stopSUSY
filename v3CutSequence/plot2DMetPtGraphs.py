#!/usr/bin/env python

# NOTE: NEEDS 7 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, displayMode, findingSameFlavor, muPreference, lastcut, plotVarX, 
# plotVarY
# True testMode plots only a few events; True displayMode displays rather than 
# saves w/o displaying the hists.
# Implements additional cuts and then draws a 2D hist for the summed bkgd data 
# and for each of the signal files (e.g. met vs. pt(l), met vs. mt(l)).
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses nentries info from bkgd_TTDiLept_files_nentries
# Uses xsec and nentries info from sig_SingleStop_files_nentries
# Possible lastcuts: listed in cuts below.

import sys
from ROOT import TFile, TTree, TH2F, TCanvas, TImage, TLegend, TPaletteAxis
from ROOT import gSystem, gStyle, gROOT, kTRUE, gPad
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
import numpy as np

assert len(sys.argv) == 8, "need 7 command line args: testMode{0,1}, displayMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}, lastcut, plotVarX, plotVarY"

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("nbtag<2",2), ("MET>80",3),\
        ("no3rdlept",4), ("njets<4",5)])
lastcut = sys.argv[5]
assert lastcut in cuts, "invalid last cut %s" % lastcut
nCuts = cuts[lastcut]+1

plotSettings = { # [nBins,xMin,xMax]
        "lep1_pt":[20,0,400,"[Gev]"],
        "lep2_pt":[20,0,400,"[GeV]"],
        "lep1_mt":[20,0,500,"[GeV]"],
        "lep2_mt":[20,0,500,"[GeV]"],
        "met_pt":[20,0,500,"[GeV]"],
        "lep1_eta":[20,-4,4,""],
        "lep2_eta":[20,-4,4,""],
        }

plotVarsXY = sys.argv[6:8] # x, y
for plotVar in plotVarsXY:
    assert (plotVar in plotSettings), "invalid plotVar %s" % plotVar

testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
displayMode = bool(int(sys.argv[2]))
print "Display mode:", displayMode
if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases
# applying cuts
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
channelName = l1Flav[:2] + l2Flav[:2]

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
bkgdNtupleAdr += "Bkgd_TTDiLept_"+channelName+".root"
sigsNtupleAdr += "Sig_"+channelName+".root"
print "Plotting",str(plotVarsXY),"from",bkgdNtupleAdr,"and",sigsNtupleAdr
print "Cutting events up to and including", lastcut

numSigFiles = int(sigsNtupleAdr[48:50])

testMode = True 
if numSigFiles > 10: testMode = False 
nBinsX = plotSettings[plotVarsXY[0]][0]
if not testMode: nBinsX = nBinsX * 5
xMin = plotSettings[plotVarsXY[0]][1]
xMax = plotSettings[plotVarsXY[0]][2]
binwidthX = (xMax - xMin)/nBinsX # include overflow bin
nBinsY = plotSettings[plotVarsXY[1]][0]
if not testMode: nBinsY = nBinsY * 5
yMin = plotSettings[plotVarsXY[1]][1]
yMax = plotSettings[plotVarsXY[1]][2]
binwidthY = (yMax - yMin)/nBinsY # include overflow bin

hBkgd = TH2F(plotVarsXY[1]+"_"+plotVarsXY[0]+"_bkgd", \
        plotVarsXY[1]+"_"+plotVarsXY[0]+"_bkgd", \
        nBinsX, xMin, xMax, nBinsY, yMin, yMax)
lumi = 3000000 # luminosity = 3000 /fb = 3,000,000 /fb

c = TCanvas("c","Plot",10,20,1000,700)
gStyle.SetPalette(1)
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.9)

#--------------------------------------------------------------------------------#
# *************** Filling bkgd data summed together  ************
print
print "Plotting " + plotVarsXY[1] + " vs. " + plotVarsXY[0] + " from background."
xsec = 67.75 # production cross section
bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
inTree = bkgdFile.Get("tBkgd")

totOrigNentries = 0
bkgdDataListFile = open("bkgd_TTDiLept_files_nentries")
for line in bkgdDataListFile:
    line = line.rstrip('\n')
    filename, origNentries = line.split(" ")
    totOrigNentries += int(origNentries)

nentries = inTree.GetEntries()
print("nentries={0:d}".format(nentries))
assert nentries > 0, "You have no events in your tree..."

# ********** Looping over events. ***********
for count, event in enumerate(inTree):
    if count % 500000 == 0: print("count={0:d}".format(count))
    # ********** Additional cuts. ***********
    if findingSameFlavor:
        if muPreference:
            lepIndices = selectMuMu(event)
            # l1Flav, l2Flav set at runtime
        else: 
            lepIndices = selectElEl(event)
            # l1Flav, l2Flav set at runtime
        if lepIndices is None: continue
    else:
        lepIndices = selectMuEl(event)
        l1Flav = "muon"
        l2Flav = "electron"
        if lepIndices is None:
            lepIndices = selectElMu(event)
            if lepIndices is None: continue
            l1Flav = "electron"
            l2Flav = "muon"
    l1Index = lepIndices[0]
    l2Index = lepIndices[1]

    # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

    if nCuts > cuts["nbtag<2"]:
        if event.nbtag > 1: continue

    if nCuts > cuts["MET>80"]:
        if event.met_pt < 80: continue

    if nCuts > cuts["no3rdlept"]:
        # veto (3rd lepton) checks:
        if findingSameFlavor:
            # event should not give valid muel or elmu pair
            if selectMuEl(event) is not None: continue
            if selectElMu(event) is not None: continue
        else:
            # event should not give valid mumu or elel pair
            if selectMuMu(event) is not None: continue
            if selectElEl(event) is not None: continue
        
    if nCuts > cuts["njets<4"]:
        if event.njets >= 4: continue

    # ********** Plotting. ***********
    valXY = []
    for plotVar in plotVarsXY:
        if plotVar[:4] == "lep1": 
            valXY.append(np.reshape(getattr(event, l1Flav+plotVar[4:]),\
                    20)[l1Index])
        elif plotVar[:4] == "lep2": 
            valXY.append(np.reshape(getattr(event, l2Flav+plotVar[4:]), \
                    20)[l2Index])
        else: valXY.append(getattr(event, plotVar))

    if valXY[0] > xMax:
        if valXY[1] > yMax: # x and y overflow
            hBkgd.Fill(xMax - binwidthX/2, yMax - binwidthY/2, 1)
        else: # x overflow, y in range
            hBkgd.Fill(xMax - binwidthX/2, valXY[1], 1)
    else:
        if valXY[1] > yMax: # x in range, y overflow
            hBkgd.Fill(valXY[0], yMax - binwidthY/2, 1)
        else: # x in range, y in range
            hBkgd.Fill(valXY[0], valXY[1], 1)
hBkgd.Sumw2()
c.cd()
title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" ("+channelName+\
        ", cuts to "+lastcut+")"
hBkgd.SetTitle(title)
unitsLabelX = plotSettings[plotVarsXY[0]][3]
unitsLabelY = plotSettings[plotVarsXY[1]][3]
hBkgd.GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
hBkgd.GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
hBkgd.Scale(xsec*lumi/totOrigNentries)
hBkgd.SetLineColor(1) # black
hBkgd.Draw("colz")
gPad.Update()
hBkgd.GetZaxis().SetLabelSize(0.02)
c.Update()
if displayMode:
    print "Done plotting bkgd 2d hist. Press enter to continue."
    raw_input()
else:
    gSystem.ProcessEvents()
    imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
            "plots/v3CutSequence/bkgd_"+plotVarsXY[1]+"_"+plotVarsXY[0]+"_"+\
            channelName+"_"+lastcut+".png"
    print "Saving image", imgName
    img = TImage.Create()
    img.FromPad(c)
    img.WriteImage(imgName)
    print "Done plotting bkgd 2d hist."
bkgdFile.Close()

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVarsXY[1] + " vs. " + plotVarsXY[0] + " from signal."
sigDataListFile = open("sig_SingleStop_files_nentries")

hSigArr = []
for fileNum, line in enumerate(sigDataListFile):
    print
    if fileNum + 1 > numSigFiles: break

    filename, xsec, origNentries = line.split(" ")
    xsec = float(xsec)
    origNentries = int(origNentries)
    print filename

    sigFile = TFile.Open(sigsNtupleAdr, "READ")
    inTree = sigFile.Get("tSig"+str(fileNum))
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    assert nentries > 0, "You have no events in your tree..."

    hSigArr.append(TH2F(plotVarsXY[1]+"_"+plotVarsXY[0]+"_sig_"+filename[21:31], \
            plotVarsXY[1]+"_"+plotVarsXY[0] +"_sig_"+filename[21:31], \
            nBinsX, xMin, xMax, nBinsY, yMin, yMax))
    
    for count, event in enumerate(inTree):
    # ********** Additional cuts. ***********
        if findingSameFlavor:
            if muPreference:
                lepIndices = selectMuMu(event)
                # l1Flav, l2Flav set at runtime
            else: 
                lepIndices = selectElEl(event)
                # l1Flav, l2Flav set at runtime
            if lepIndices is None: continue
        else:
            lepIndices = selectMuEl(event)
            l1Flav = "muon"
            l2Flav = "electron"
            if lepIndices is None:
                lepIndices = selectElMu(event)
                if lepIndices is None: continue
                l1Flav = "electron"
                l2Flav = "muon"
        l1Index = lepIndices[0]
        l2Index = lepIndices[1]

        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue

        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue

        if nCuts > cuts["no3rdlept"]:
            # veto (3rd lepton) checks:
            if findingSameFlavor:
                # event should not give valid muel or elmu pair
                if selectMuEl(event) is not None: continue
                if selectElMu(event) is not None: continue
            else:
                # event should not give valid mumu or elel pair
                if selectMuMu(event) is not None: continue
                if selectElEl(event) is not None: continue
            
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue

        # ********** Plotting. ***********
        valXY = []
        for plotVar in plotVarsXY:
            if plotVar[:4] == "lep1": 
                valXY.append(np.reshape(getattr(event, l1Flav+plotVar[4:]),\
                        20)[l1Index])
            elif plotVar[:4] == "lep2": 
                valXY.append(np.reshape(getattr(event, l2Flav+plotVar[4:]), \
                        20)[l2Index])
            else: valXY.append(getattr(event, plotVar))

        if valXY[0] > xMax:
            if valXY[1] > yMax: # x and y overflow
                hSigArr[fileNum].Fill(xMax - binwidthX/2, yMax - binwidthY/2, 1)
            else: # x overflow, y in range
                hSigArr[fileNum].Fill(xMax - binwidthX/2, valXY[1], 1)
        else:
            if valXY[1] > yMax: # x in range, y overflow
                hSigArr[fileNum].Fill(valXY[0], yMax - binwidthY/2, 1)
            else: # x in range, y in range
                hSigArr[fileNum].Fill(valXY[0], valXY[1], 1)

    hSigArr[fileNum].Sumw2()
    title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" ("+channelName+\
            ", cuts to "+lastcut+")"
    hSigArr[fileNum].SetTitle(title)
    unitsLabelX = plotSettings[plotVarsXY[0]][3]
    unitsLabelY = plotSettings[plotVarsXY[1]][3]
    hSigArr[fileNum].GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
    hSigArr[fileNum].GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
    hSigArr[fileNum].Scale(xsec*lumi/origNentries)
    hSigArr[fileNum].Draw("colz")
    hSigArr[fileNum].GetZaxis().SetLabelSize(0.02)
    c.Update()
    if displayMode:
        print "Done plotting sig%i 2d hist. Press enter to continue." %fileNum
        raw_input()
    else:
        gSystem.ProcessEvents()
        imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
                "plots/v3CutSequence/sig"+str(fileNum)+"_"+plotVarsXY[1]+"_"+\
                plotVarsXY[0]+"_"+channelName+"_"+lastcut+".png"
        print "Saving image", imgName
        img = TImage.Create()
        img.FromPad(c)
        img.WriteImage(imgName)
        print "Done plotting sig%i 2d hist." %fileNum

    sigFile.Close()


#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************

# legend = TLegend(.65,.75,.90,.90)
# legend.AddEntry(hBkgd, plotVarsXY[1]+"_"+plotVarsXY[0]+"_bkgd")
# legend.SetTextSize(0.02)
# for fileNum, h in enumerate(hSigArr):
#     legend.AddEntry(hSigArr[fileNum], hSigArr[fileNum].GetTitle())
# legend.Draw("same")

print "Done."

