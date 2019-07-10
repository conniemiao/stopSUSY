#!/usr/bin/env python

# NOTE: NEEDS 8 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, displayMode, findingSameFlavor, muPreference, lastcut, process, 
# plotVarX, plotVarY
# True testMode plots only a few events; True displayMode displays rather than 
# saving w/o displaying the hists.
# Implements additional cuts and then draws a 2D hist for the chosen process of
# bkgd data and for each of the signal files (e.g. met vs. pt(l), met vs. mt(l)).
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from bkgd_files
# Uses xsec info from sig_SingleStop_files
# Possible lastcuts: listed in cuts below.

import sys
from ROOT import TFile, TTree, TH2F, TCanvas, TImage, TLegend, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE, gPad
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
from math import sqrt, cos
import numpy as np
import time

assert len(sys.argv) == 9, "need 8 command line args: testMode{0,1}, displayMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}, lastcut, process, plotVarX, plotVarY"

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("nbtag<2",2), ("MET>80",3),\
        ("no3rdlept",4), ("njets<4",5)])
lastcut = sys.argv[5]
assert lastcut in cuts, "invalid last cut %s" % lastcut
nCuts = cuts[lastcut]+1

plotSettings = { # [nBins,xMin,xMax]
        "lep1_pt":[50,0,400,"[Gev]"],
        "lep2_pt":[50,0,400,"[GeV]"],
        "lep1_mt":[50,0,500,"[GeV]"],
        "lep2_mt":[50,0,500,"[GeV]"],
        "met_pt":[50,0,500,"[GeV]"],
        "lep1_eta":[50,-4,4,""],
        "lep2_eta":[50,-4,4,""],
        "jet_ht":[50,0,800,"[GeV]"],
        "mt_tot":[50,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        "mt_sum":[50,0,1000,"[GeV]"], # mt1 + mt2
        "m_eff":[50,0,1000,"[GeV]"], # ht + met + pt1 + pt2
        }

plotVarsXY = sys.argv[7:9] # x, y
for plotVar in plotVarsXY:
    assert (plotVar in plotSettings), "invalid plotVar %s" % plotVar

testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
displayMode = bool(int(sys.argv[2]))
print "Display mode:", displayMode
if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases
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

processes = {"W-Jets", "Drell-Yan", "Diboson", "Single-Top", "TT+X"}

thisProcess = sys.argv[6]
assert thisProcess in processes, "invalid process %s" % thisProcess
print "Bkgd process:", thisProcess

# assemble the sigsNtupleAdr and bkgdNtupleAdr
# number of files to process
numBkgdFiles = 27  # need to loop over all the files in order to have correct xsec
if testMode: 
    numBkgdFiles = 2 
numSigFiles = 2 # just use the first signal one
baseDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/"
print "Plotting",str(plotVarsXY)
print "Cutting events up to and including", lastcut

nBinsX = plotSettings[plotVarsXY[0]][0]
xMin = plotSettings[plotVarsXY[0]][1]
xMax = plotSettings[plotVarsXY[0]][2]
binwidthX = (xMax - xMin)/nBinsX
nBinsY = plotSettings[plotVarsXY[1]][0]
yMin = plotSettings[plotVarsXY[1]][1]
yMax = plotSettings[plotVarsXY[1]][2]
binwidthY = (yMax - yMin)/nBinsY

hBkgdSubprocessesDict = {} 
# 1 hBkgd for each subprocess for this process which will be stacked into 1 hstack
with open("bkgd_files") as bkgdSubprocessesListFile:
    for subprocessLine in bkgdSubprocessesListFile:
        subprocessLine = subprocessLine.rstrip('\n')
        subprocess, process, xsec = subprocessLine.split(" ")
        if subprocess[0] == "#": continue # problematic input files
        if process != thisProcess: continue
        hBkgd = TH2F("bkgd_"+subprocess, "bkgd_"+subprocess, nBinsX, xMin, xMax, \
                nBinsY, yMin, yMax)
        hBkgd.SetDirectory(0) # necessary to keep hist from closing
        hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
        hBkgdSubprocessesDict.update({subprocess:hBkgd})

lumi = 3000000 # luminosity = 3 /ab = 3000 /fb = 3,000,000 /fb

c = TCanvas("c","Plot",10,20,1000,700) # same canvas used for signal and bkgd
gStyle.SetPalette(1)
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.9)

start_time = time.time()

#--------------------------------------------------------------------------------#
# *************** Filling bkgd data summed together  ************
print
print "Plotting", plotVarsXY[1], "vs.", plotVarsXY[0], "from bkgd", thisProcess

title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" ("+thisProcess+" bkgd, "+channelName\
        +", cuts to "+lastcut+")"
hBkgdStack = THStack("hBkgdStack", title)
bkgdSubprocessesListFile = open("bkgd_files")
firstFile = True # indicates if this is the first subprocess file for this process
for subprocessLine in bkgdSubprocessesListFile:
    subprocessLine = subprocessLine.rstrip('\n')
    subprocess, process, xsec = subprocessLine.split(" ")
    xsec = float(xsec)

    if subprocess[0] == "#": continue # problematic input files
    if process != thisProcess: continue # only take subprocess from this process
    if not process in processes: continue

    # assemble the bkgdNtupleAdr
    bkgdNtupleAdr = baseDir+"stopCut_"
    if testMode: bkgdNtupleAdr += "test_"
    else: bkgdNtupleAdr += "all_"
    bkgdNtupleAdr += "Bkgd_"+subprocess+"_"+channelName+".root"
    print "Plotting from", bkgdNtupleAdr

    bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
    tBkgd = bkgdFile.Get("tBkgd")

    nentries = tBkgd.GetEntries()
    print("nentries={0:d}".format(nentries))
    assert nentries > 0, "You have no events in your tree..."
    
    hBkgdGenweights = bkgdFile.Get("genweights")
    # tot for this subprocess:
    bkgdTotGenweight = hBkgdGenweights.GetSumOfWeights()

    hBkgd = hBkgdSubprocessesDict[subprocess]

    nMax = 10000
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        if count > nMax : break
        if count % 100000 == 0: print("count={0:d}".format(count))
        genwt = event.genweight
    
        if not findingSameFlavor: # if findingSameFlavor, l1/l2Flav set at runtime
            if event.lep1_isMu: l1Flav = "muon"
            else: l1Flav = "electron"
            if event.lep2_isMu: l1Flav = "muon"
            else: l2Flav = "electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
        # ********** Additional cuts. ***********
        if nCuts > cuts["dilepton"]:
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
            if event.found3rdLept: continue
            
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue
    
        # ********** Plotting. ***********
        valXY = []
        if plotVar == "jet_ht" and event.njets == 0: continue # to next event
        for plotVar in plotVarsXY:
            if plotVar[:4] == "lep1": 
                valXY.append(np.reshape(getattr(event, l1Flav+plotVar[4:]),\
                        20)[l1Index])
            elif plotVar[:4] == "lep2": 
                valXY.append(np.reshape(getattr(event, l2Flav+plotVar[4:]), \
                        20)[l2Index])
            elif plotVar == "jet_ht":
                valXY.append(event.jet_ht)
            elif plotVar[:6] == "mt_tot":
                valXY.append(sqrt((np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index])**2 + (np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index])**2))
            elif plotVar == "mt_sum":
                valXY.append(np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index])
            elif plotVar[:5] == "m_eff":
                val = event.met_pt + np.reshape(getattr(event, \
                        l1Flav+"_pt"),20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_pt"),20)[l2Index]
                if event.njets > 0: val += event.jet_ht
                valXY.append(val)
            else: valXY.append(getattr(event, plotVar))
    
        if valXY[0] > xMax:
            if valXY[1] > yMax: # x and y overflow
                hBkgd.Fill(xMax - binwidthX/2, yMax - binwidthY/2, genwt)
            else: # x overflow, y in range
                hBkgd.Fill(xMax - binwidthX/2, valXY[1], genwt)
        else:
            if valXY[1] > yMax: # x in range, y overflow
                hBkgd.Fill(valXY[0], yMax - binwidthY/2, genwt)
            else: # x in range, y in range
                hBkgd.Fill(valXY[0], valXY[1], genwt)
    # hBkgd.Sumw2()
    c.cd()
    hBkgd.Scale(xsec*lumi/bkgdTotGenweight)
    hBkgdStack.Add(hBkgd)

    print
    bkgdFile.Close()

# hBkgdStack.Draw("colz")
# hBkgdStack.Draw("lego2")
hBkgdStack.Draw("col")
unitsLabelX = plotSettings[plotVarsXY[0]][3]
unitsLabelY = plotSettings[plotVarsXY[1]][3]
hBkgdStack.GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
hBkgdStack.GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
gPad.Update()
c.Update()

if displayMode:
    print "Done plotting bkgd 2d hist. Press enter to continue."
    raw_input()
else:
    gSystem.ProcessEvents()
    imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
            "plots/v3CutSequence/bkgd_"+thisProcess+"_"+plotVarsXY[1]+"_v_"+\
            plotVarsXY[0]+"_"+channelName+"_"+lastcut+".png"
    print "Saving image", imgName
    img = TImage.Create()
    img.FromPad(c)
    img.WriteImage(imgName)
    print "Done plotting bkgd 2d hist."

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVarsXY[1] + " vs. " + plotVarsXY[0] + " from signal."
# assemble the sigsNtupleAdr
sigsNtupleAdr = baseDir+"stopCut_"
if numSigFiles < 10: sigsNtupleAdr += "0"+str(numSigFiles)
else: sigsNtupleAdr += str(numSigFiles)
sigsNtupleAdr += "Sig_"+channelName+".root"
sigDataListFile = open("sig_SingleStop_files")

hSigArr = []
for fileNum, line in enumerate(sigDataListFile):
    print
    if fileNum + 1 > numSigFiles: break
    if fileNum + 1 > 1: break # just want to plot from 1 sig for now

    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename

    sigFile = TFile.Open(sigsNtupleAdr, "READ")
    inTree = sigFile.Get("tSig"+str(fileNum))
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    assert nentries > 0, "You have no events in your tree..."

    hSigArr.append(TH2F("sig_"+filename[18:31], "sig_"+filename[18:31], \
            nBinsX, xMin, xMax, nBinsY, yMin, yMax))
    hSigArr[fileNum].SetDefaultSumw2() # automatically sum w^2 while filling

    hSigGenweights = sigFile.Get("genweights")
    sigTotGenweight = hSigGenweights.GetSumOfWeights()
    
    # ********** Looping over events. ***********
    for count, event in enumerate(inTree):
        if count % 100000 == 0: print("count={0:d}".format(count))
        genwt = event.genweight

        if not findingSameFlavor: # if findingSameFlavor, l1/l2Flav set at runtime
            if event.lep1_isMu: l1Flav = "muon"
            else: l1Flav = "electron"
            if event.lep2_isMu: l1Flav = "muon"
            else: l2Flav = "electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index

        # ********** Additional cuts. ***********
        if nCuts > cuts["dilepton"]:
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
            if event.found3rdLept: continue
            
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue

        # ********** Plotting. ***********
        valXY = []
        if plotVar == "jet_ht" and event.njets == 0: continue # to next event
        for plotVar in plotVarsXY:
            if plotVar[:4] == "lep1": 
                valXY.append(np.reshape(getattr(event, l1Flav+plotVar[4:]),\
                        20)[l1Index])
            elif plotVar[:4] == "lep2": 
                valXY.append(np.reshape(getattr(event, l2Flav+plotVar[4:]), \
                        20)[l2Index])
            elif plotVar[:6] == "jet_ht":
                valXY.append(event.jet_ht)
            elif plotVar[:6] == "mt_tot":
                valXY.append(sqrt((np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index])**2 + (np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index])**2))
            elif plotVar == "mt_sum":
                valXY.append(np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index])
            elif plotVar[:5] == "m_eff":
                val = event.met_pt + np.reshape(getattr(event, \
                        l1Flav+"_pt"),20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_pt"),20)[l2Index]
                if event.njets > 0: val += event.jet_ht
                valXY.append(val)

        if valXY[0] > xMax:
            if valXY[1] > yMax: # x and y overflow
                hSigArr[fileNum].Fill(xMax-binwidthX/2, yMax-binwidthY/2, genwt)
            else: # x overflow, y in range
                hSigArr[fileNum].Fill(xMax-binwidthX/2, valXY[1], genwt)
        else:
            if valXY[1] > yMax: # x in range, y overflow
                hSigArr[fileNum].Fill(valXY[0], yMax-binwidthY/2, genwt)
            else: # x in range, y in range
                hSigArr[fileNum].Fill(valXY[0], valXY[1], genwt)

    # hSigArr[fileNum].Sumw2()
    title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" (sig_"+filename[18:31]+", "+\
            channelName+", cuts to "+lastcut+")"
    hSigArr[fileNum].SetTitle(title)
    unitsLabelX = plotSettings[plotVarsXY[0]][3]
    unitsLabelY = plotSettings[plotVarsXY[1]][3]
    hSigArr[fileNum].GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
    hSigArr[fileNum].GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
    hSigArr[fileNum].Scale(xsec*lumi/sigTotGenweight)
    hSigArr[fileNum].Draw("colz")
    hSigArr[fileNum].GetZaxis().SetLabelSize(0.02)
    c.Update()
    if displayMode:
        print "Done plotting sig%i 2d hist. Press enter to continue." %fileNum
        raw_input()
    else:
        gSystem.ProcessEvents()
        imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
                "plots/v3CutSequence/sig"+str(fileNum)+"_"+plotVarsXY[1]+"_v_"+\
                plotVarsXY[0]+"_"+channelName+"_"+lastcut+".png"
        print "Saving image", imgName
        img = TImage.Create()
        img.FromPad(c)
        img.WriteImage(imgName)
        print "Done plotting sig%i 2d hist." %fileNum

    sigFile.Close()


#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print int(time.time()-start_time), "secs of processing."

# legend = TLegend(.65,.75,.90,.90)
# legend.AddEntry(hBkgd, plotVarsXY[1]+"_"+plotVarsXY[0]+"_bkgd")
# legend.SetTextSize(0.02)
# for fileNum, h in enumerate(hSigArr):
#     legend.AddEntry(hSigArr[fileNum], hSigArr[fileNum].GetTitle())
# legend.Draw("same")

print "Done."

