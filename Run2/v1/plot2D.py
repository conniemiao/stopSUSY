#!/usr/bin/env python

# NOTE: NEEDS 7 CMD LINE ARGS with values:
# testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut,
# process, plotVarX, plotVarY
#
# Implements additional cuts and then draws a 2D hist for the chosen process of
# bkgd data, for each of the signal files, and for data (e.g. met vs. pt(l), met vs. 
# mt(l)).
# ** WILL NOT SUPPORT PLOTTING WJETS OR DYJETS **
#
# Uses the root files outputted by hadding the output from makeNtuple.py
# Uses xsec info from bkgd_fileRedirector
# Uses xsec info from sig_fileRedirector
# Possible lastcuts: listed in cuts below.

print "Importing modules."
import sys, os
from ROOT import TFile, TTree, TH2F, TCanvas, TImage, TLegend, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE, gPad
from stopSelection import deltaR
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
from math import sqrt, cos
import numpy as np
import time
print "Beginning execution of", sys.argv

assert len(sys.argv) == 8, "need 7 command line args: testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut, process, plotVarX, plotVarY"

if sys.argv[1] == "test": testMode = True
elif sys.argv[1] == "all": testMode = False
else: assert False, "invalid test mode, need {test, all}"

if sys.argv[2] == "show": displayMode = True
elif sys.argv[2] == "save": displayMode = False
else: assert False, "invalid display mode, need {show, save}"
if not displayMode: gROOT.SetBatch(kTRUE) # prevent displaying canvases

# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
if sys.argv[3] == "mumu":
    findingSameFlavor = True
    muPreference = True
    l1Flav = "Muon"
    l2Flav = "Muon"
    dataProcess = "DoubleMuon"
elif sys.argv[3] == "elel":
    findingSameFlavor = True
    muPreference = False
    l1Flav = "Electron"
    l2Flav = "Electron"
    dataProcess = "DoubleEG"
elif sys.argv[3] == "muel":
    findingSameFlavor = False
    muPreference = False
    l1Flav = "Muon"
    l2Flav = "Electron"
    dataProcess = "MuonEG"
else: assert False, "invalid channel, need {mumu, elel, muel}"
channelName = l1Flav[:2] + l2Flav[:2]

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
        ("MET>80",4),("nJet<4",5)])
lastcut = sys.argv[4]
assert lastcut in cuts, "invalid last cut %s" % lastcut
nCuts = cuts[lastcut]+1

plotSettings = { # [nBins,xMin,xMax]
        "lep1_pt":[50,0,400,"[Gev]"],
        "lep2_pt":[50,0,400,"[GeV]"],
        "lep1_mt":[50,0,500,"[GeV]"],
        "lep2_mt":[50,0,500,"[GeV]"],
        "MET_pt":[50,0,500,"[GeV]"],
        "lep1_eta":[50,-4,4,""],
        "lep2_eta":[50,-4,4,""],
        "Jet_ht":[50,0,800,"[GeV]"],
        "mt_tot":[50,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        "mt_sum":[50,0,1000,"[GeV]"], # mt1 + mt2
        "m_eff":[50,0,1000,"[GeV]"], # ht + met + pt1 + pt2
        }

processes = {"TTBar":30, "Diboson":41, "Single-Top":40, "TT+X":7}
thisProcess = sys.argv[5]
assert thisProcess in processes, "invalid process %s" % thisProcess
print "Bkgd process:", thisProcess

plotVarsXY = sys.argv[6:8] # x, y
for plotVar in plotVarsXY:
    assert (plotVar in plotSettings), "invalid plotVar %s" % plotVar

myDataDir = "/eos/user/c/cmiao/private/myDataSusy/Run2/"
# limit the number of files to process (other than what is commented out in the file
# redirector)
numBkgdFiles = 1
numSigFiles = 1
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
        "plots/Run2/v1/plot1D/"
if not os.path.exists(imgDir): os.makedirs(imgDir) 

nBinsX = plotSettings[plotVarsXY[0]][0]
xMin = plotSettings[plotVarsXY[0]][1]
xMax = plotSettings[plotVarsXY[0]][2]
binwidthX = (xMax - xMin)/nBinsX
nBinsY = plotSettings[plotVarsXY[1]][0]
yMin = plotSettings[plotVarsXY[1]][1]
yMax = plotSettings[plotVarsXY[1]][2]
binwidthY = (yMax - yMin)/nBinsY
unitsLabelX = plotSettings[plotVarsXY[0]][3]
unitsLabelY = plotSettings[plotVarsXY[1]][3]

c = TCanvas("c","Plot",10,20,1000,700) # same canvas used for signal and bkgd
gStyle.SetPalette(1)
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.9)

lumi = 35921 # 2016 lumi in /pb

#--------------------------------------------------------------------------------#
# *************** Filling bkgd data summed together  ************
print
print "----------- Plotting from background. -----------"

hBkgdSubprocessesDict = {} 
# 1 hBkgd for each subprocess for this process which will be stacked into 1 hstack
with open("bkgd_fileRedirector") as bkgdSubprocessesListFile:
    for subprocessLine in bkgdSubprocessesListFile:
        subprocessLine = subprocessLine.rstrip('\n').split(" ")
        subprocess = subprocessLine[0]
        if subprocess[0] == "#": continue
        process = subprocessLine[1]
        if process != thisProcess: continue
        hBkgd = TH2F("bkgd_"+subprocess, "bkgd_"+subprocess, nBinsX, xMin, xMax, \
                nBinsY, yMin, yMax)
        hBkgd.SetDirectory(0) # necessary to keep hist from closing
        hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
        hBkgdSubprocessesDict.update({subprocess:hBkgd})

title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" ("+thisProcess+" bkgd, "+channelName\
        +", cuts to "+lastcut+")"
hBkgdStack = THStack("hBkgdStack", title)

bkgdSubprocessesListFile = open("bkgd_fileRedirector")
for subprocessLine in bkgdSubprocessesListFile:
    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    if subprocess[0] == "#": continue
    process = subprocessLine[1]
    if process != thisProcess: continue
    xsec = float(subprocessLine[2])

    # assemble the bkgdNtupleAdr
    bkgdNtupleAdr = myDataDir+"bkgd/"+process+"/"+subprocess+"/"+subprocess+"_"
    if testMode: bkgdNtupleAdr += "test_"
    else: bkgdNtupleAdr += "all_"
    bkgdNtupleAdr += channelName+".root"
    print "Plotting from", bkgdNtupleAdr

    try:
        bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
        tBkgd = bkgdFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+bkgdNtupleAdr+\
                ", skipping\n")
        continue

    try: nentries = tBkgd.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+bkgdNtupleAdr+\
                ", skipping\n")
        continue
    print("nentries={0:d}".format(nentries))
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+bkgdNtupleAdr+" has no entries!"+\
                " Skipping\n")
        continue
    
    hBkgdGenweights = bkgdFile.Get("genWeights")
    # tot for this subprocess:
    bkgdTotGenweight = hBkgdGenweights.GetSumOfWeights()

    hBkgd = hBkgdSubprocessesDict[subprocess]

    nMax = nentries
    if testMode: nMax = 10000
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        if count > nMax : break
        if count == 0: start_time = time.time()
        if count % 100000 == 0: print "count =", count
        genwt = event.genWeight
        puwt = event.puWeight
        evtwt = genwt*puwt
    
        # ********** Additional cuts. ***********

        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index

        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if list(getattr(event, l1Flav+"_charge"))[l1Index] * \
                    list(getattr(event, l2Flav+"_charge"))[l2Index] >= 0: continue
            if findingSameFlavor:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.1: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.1: continue
            else:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.2: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.2: continue
    
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
    
        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
    
        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
    
        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue
            
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
    
        # ********** Plotting. ***********
        valXY = []
        if "Jet_ht" in plotVarsXY and event.nJet == 0: continue # to next event
        for plotVar in plotVarsXY:
            if plotVar[:4] == "lep1":
                valXY.append(list(getattr(event, l1Flav+plotVar[4:]))[l1Index])
            elif plotVar[:4] == "lep2":
                valXY.append(list(getattr(event, l2Flav+plotVar[4:]))[l2Index])
            else:
                valXY.append(getattr(event, plotVar))
    
        if valXY[0] > xMax:
            if valXY[1] > yMax: # x and y overflow
                hBkgd.Fill(xMax - binwidthX/2, yMax - binwidthY/2, evtwt)
            else: # x overflow, y in range
                hBkgd.Fill(xMax - binwidthX/2, valXY[1], evtwt)
        else:
            if valXY[1] > yMax: # x in range, y overflow
                hBkgd.Fill(valXY[0], yMax - binwidthY/2, evtwt)
            else: # x in range, y in range
                hBkgd.Fill(valXY[0], valXY[1], evtwt)
    # hBkgd.Sumw2()
    c.cd()
    hBkgd.Scale(xsec*lumi/bkgdTotGenweight)
    hBkgdStack.Add(hBkgd)

    bkgdFile.Close()

# hBkgdStack.Draw("colz")
# hBkgdStack.Draw("lego2")
hBkgdStack.Draw("col")
hBkgdStack.GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
hBkgdStack.GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
gPad.Update()
c.Update()

if displayMode:
    print "Done plotting bkgd 2d hist. Press enter to continue."
    raw_input()
else:
    gSystem.ProcessEvents()
    imgName = imgDir+"bkgd_"+thisProcess+"_"+plotVarsXY[1]+\
            "_v_"+plotVarsXY[0]+"_"+channelName+"_"+lastcut+".png"
    print "Saving image", imgName
    img = TImage.Create()
    img.FromPad(c)
    img.WriteImage(imgName)
    print "Done plotting bkgd 2d hist."
    print

#--------------------------------------------------------------------------------#
# *************** Filling each signal in a separate hist  ************
print
print "----------- Plotting from signal. -----------"

sig_redirector = open("sig_fileRedirector")

hSigDict = {}

for fileNum, subprocessLine in enumerate(sig_redirector):
    if fileNum + 1 > numSigFiles: break

    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    process = subprocessLine[1]
    xsec = float(subprocessLine[2])
    print subprocess 
    
    # assemble the sigNtupleAdr
    sigNtupleAdr = myDataDir+"sig/"+process+"/"+subprocess+"/"+subprocess+"_"
    if testMode: sigNtupleAdr += "test_"
    else: sigNtupleAdr += "all_"
    sigNtupleAdr += channelName+".root"

    try:
        sigFile = TFile.Open(sigNtupleAdr, "READ")
        tSig = sigFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+sigNtupleAdr+\
                ", skipping\n")
        continue
    try:
        nentries = tSig.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+sigNtupleAdr+\
                ", skipping\n")
        continue
    print("nentries={0:d}".format(nentries))
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+sigNtupleAdr+" has no entries!"+\
                " Skipping\n")
        continue

    hSig = TH2F("sig_"+subprocess[10:27], "sig_"+subprocess[10:27], nBinsX, \
            xMin, xMax, nBinsY, yMin, yMax)
    hSig.SetDefaultSumw2() # automatically sum w^2 while filling
    hSig.SetDirectory(0)
    hSigDict.update({subprocess:hSig})

    hSigGenweights = sigFile.Get("genWeights")
    sigTotGenweight = hSigGenweights.GetSumOfWeights()
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tSig):
        if count % 100000 == 0: print "count =", count
        genwt = event.genWeight
        puwt = event.puWeight
        evtwt = genwt*puwt

        # ********** Additional cuts. ***********

        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
    
        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if list(getattr(event, l1Flav+"_charge"))[l1Index] * \
                    list(getattr(event, l2Flav+"_charge"))[l2Index] >= 0: continue
            if findingSameFlavor:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.1: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.1: continue
            else:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.2: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.2: continue

        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue

        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
        
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue

        # ********** Plotting. ***********
        valXY = []
        if "Jet_ht" in plotVarsXY and event.nJet == 0: continue # to next event
        for plotVar in plotVarsXY:
            if plotVar[:4] == "lep1":
                valXY.append(list(getattr(event, l1Flav+plotVar[4:]))[l1Index])
            elif plotVar[:4] == "lep2":
                valXY.append(list(getattr(event, l2Flav+plotVar[4:]))[l2Index])
            else:
                valXY.append(getattr(event, plotVar))

        if valXY[0] > xMax:
            if valXY[1] > yMax: # x and y overflow
                hSig.Fill(xMax-binwidthX/2, yMax-binwidthY/2, evtwt)
            else: # x overflow, y in range
                hSig.Fill(xMax-binwidthX/2, valXY[1], evtwt)
        else:
            if valXY[1] > yMax: # x in range, y overflow
                hSig.Fill(valXY[0], yMax-binwidthY/2, evtwt)
            else: # x in range, y in range
                hSig.Fill(valXY[0], valXY[1], evtwt)

    sigFile.Close()

    # hSig.Sumw2()
    title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" (sig_"+subprocess[10:27]+", "+\
            channelName+", cuts to "+lastcut+")"
    hSig.SetTitle(title)
    hSig.GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
    hSig.GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
    hSig.Scale(xsec*lumi/sigTotGenweight)
    hSig.Draw("colz")
    hSig.GetZaxis().SetLabelSize(0.02)
    c.Update()

    if displayMode:
        print "Done plotting sig", subprocess,"2d hist. Press enter to continue."
        raw_input()
    else:
        gSystem.ProcessEvents()
        imgName = imgDir+"sig_"+subprocess+"_"+plotVarsXY[1]+\
                "_v_"+plotVarsXY[0]+"_"+channelName+"_"+lastcut+".png"
        print "Saving image", imgName
        img = TImage.Create()
        img.FromPad(c)
        img.WriteImage(imgName)
        print "Done plotting sig", subprocess, "2d hist."

#--------------------------------------------------------------------------------#
# *************** Filling data in a separate hist  ************
print
print "----------- Plotting from data. -----------"
data_redirector = open("data_fileRedirector")

hData = TH2F("data", "data", nBinsX, xMin, xMax, nBinsY, yMin, yMax)
hData.SetDefaultSumw2() # automatically sum w^2 while filling
hData.SetDirectory(0)

for fileNum, subprocessLine in enumerate(data_redirector):
    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    if subprocess[0] == "#": continue
    process = subprocessLine[1]
    if process != dataProcess: continue

    # assemble the dataNtupleAdr
    dataNtupleAdr = myDataDir+"data/"+process+"/"+subprocess+"/"+subprocess+"_"
    if testMode: dataNtupleAdr += "test_"
    else: dataNtupleAdr += "all_"
    dataNtupleAdr += channelName+".root"
    print dataNtupleAdr
    
    try:
        dataFile = TFile.Open(dataNtupleAdr, "READ")
        tData = dataFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+dataNtupleAdr+\
                ", skipping\n")
        continue
    try:
        nentries = tData.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+dataNtupleAdr+\
                ", skipping\n")
        continue
    print("nentries={0:d}".format(nentries))
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+dataNtupleAdr+" has no entries!"+\
                " Skipping\n")
        continue
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tData):
        if count % 100000 == 0: print "count =", count
    
        # ********** Additional cuts. ***********
    
        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
    
        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if list(getattr(event, l1Flav+"_charge"))[l1Index] * \
                    list(getattr(event, l2Flav+"_charge"))[l2Index] >= 0: continue
            if findingSameFlavor:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.1: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.1: continue
            else:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.2: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.2: continue
    
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
    
        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
    
        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
    
        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
        
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
    
        # ********** Plotting. ***********
        valXY = []
        if "Jet_ht" in plotVarsXY and event.nJet == 0: continue # to next event
        for plotVar in plotVarsXY:
            if plotVar[:4] == "lep1":
                valXY.append(list(getattr(event, l1Flav+plotVar[4:]))[l1Index])
            elif plotVar[:4] == "lep2":
                valXY.append(list(getattr(event, l2Flav+plotVar[4:]))[l2Index])
            else:
                valXY.append(getattr(event, plotVar))

        if valXY[0] > xMax:
            if valXY[1] > yMax: # x and y overflow
                hData.Fill(xMax-binwidthX/2, yMax-binwidthY/2, 1)
            else: # x overflow, y in range
                hData.Fill(xMax-binwidthX/2, valXY[1], 1)
        else:
            if valXY[1] > yMax: # x in range, y overflow
                hData.Fill(valXY[0], yMax-binwidthY/2, 1)
            else: # x in range, y in range
                hData.Fill(valXY[0], valXY[1], 1)
    dataFile.Close()
    
title = plotVarsXY[1]+" v. "+plotVarsXY[0]+" (data, "+channelName+", cuts to "+\
        lastcut+")"
hData.SetTitle(title)
hData.GetXaxis().SetTitle(plotVarsXY[0]+" "+unitsLabelX)
hData.GetYaxis().SetTitle(plotVarsXY[1]+" "+unitsLabelY)
hData.Draw("colz")
hData.GetZaxis().SetLabelSize(0.02)
c.Update()

if displayMode:
    print "Done plotting data 2d hist. Press enter to continue."
    raw_input()
else:
    gSystem.ProcessEvents()
    imgName = imgDir+"data_"+plotVarsXY[1]+"_v_"+plotVarsXY[0]+"_"+channelName+\
            "_"+lastcut+".png"
    print "Saving image", imgName
    img = TImage.Create()
    img.FromPad(c)
    img.WriteImage(imgName)
    print "Done plotting data 2d hist."

# *************** Wrap up. *******************
print int(time.time()-start_time), "secs of processing."

print "Done."

