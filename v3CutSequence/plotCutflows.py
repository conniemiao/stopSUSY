#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, displayMode, findingSameFlavor, muPreference
# True testMode plots only a few events; True displayMode displays rather than 
# saves w/o displaying the hists.
# Implements additional cuts and then draws the cutflow. 
# Uses xsec info from bkgd_files
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from sig_SingleStop_files

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TImage, TLegend, TText, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
import numpy as np
import time

assert len(sys.argv) == 5, "needs 4 command line args: testMode{0,1}, displayMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1} ..."

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("nbtag<2",2), ("MET>80",3),\
        ("no3rdlept",4), ("njets<4",5)])
nCuts = len(cuts)

# Determining adr of bkgd and sig ntuples.
# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
displayMode = bool(int(sys.argv[2]))
print "Display mode:", displayMode
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
findingSameFlavor = bool(int(sys.argv[3]))
print "Finding same flavor:", findingSameFlavor
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
muPreference = bool(int(sys.argv[4]))
print "Mu preference:", muPreference
channelName = ""
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

# bkgd process name : color for plotting
processes = OrderedDict([("TT+X",30), ("Diboson",38), ("W-Jets",41), \
        ("Drell-Yan",46), ("Single-Top",40)])
if testMode:
    processes = OrderedDict([("TT+X",30), ("Diboson",38)])

nEvtsLabels = []

if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases

baseDir = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/"
# number of files to process
numBkgdFiles = float("inf")  # note: must loop over all files to have correct xsec
if testMode: 
    numBkgdFiles = 2 
numSigFiles = 3 # max 25

#--------------------------------------------------------------------------------#
start_time = time.time()
# *************** Filling bkgd data summed together  ************
print
print "Plotting from background."
print

# hBkgdDict maps every subprocess to an hBkgd which contains data from all the 
# ntuples for that subprocess.
hBkgdDict = {}
# hBkgdCutsCountDict maps every process to an array of size nCuts that keeps track
# of the num evts remaining after each cut for that process
hBkgdCutsCountDict = {}
with open("bkgd_files") as bkgdSubprocessesListFile:
    for subprocessLine in bkgdSubprocessesListFile:
        subprocessLine = subprocessLine.rstrip('\n')
        subprocess, process, xsec = subprocessLine.split(" ")
        if subprocess[0] == "#": continue # problematic input files
        hBkgd = TH1F("cutflow_"+subprocess+"_bkgd", \
                "cutflow_"+subprocess+"_bkgd", nCuts, 0, nCuts)
        hBkgd.SetDirectory(0) # necessary to keep hist from closing
        hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
        hBkgdDict.update({subprocess:hBkgd})
for process in processes:
    hBkgdCutsCountDict.update({process:[0]*nCuts})
c = TCanvas("c","c",10,20,1000,700)
legend = TLegend(.65,.75,.90,.90)
title = "cutflow ("+channelName+")"
hBkgdStack = THStack("cutflow_bkgdStack", title)

lumi = 3000000 # luminosity = 3000 /fb = 3,000,000 /pb

gStyle.SetOptStat(0) # don't show any stats

# ********** Looping over each subprocess. ***********
prevProcess = "" # to determine when you got to the next process
processNum = 0
bkgdSubprocessesListFile = open("bkgd_files")
for subprocessLine in bkgdSubprocessesListFile:
    subprocessLine = subprocessLine.rstrip('\n')
    subprocess, process, xsec = subprocessLine.split(" ")
    xsec = float(xsec)

    if subprocess[0] == "#": continue # problematic input files
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

    hBkgd = hBkgdDict[subprocess] 
    for i, cut in enumerate(cuts, start=1):
        if i>nCuts: break
        hBkgd.GetXaxis().SetBinLabel(i, cut)
    
    hBkgdGenweights = bkgdFile.Get("genweights")
    # tot for this subprocess:
    bkgdSubprocessGenweight = hBkgdGenweights.GetSumOfWeights()
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        if count % 100000 == 0: print("count={0:d}".format(count))
        genwt = event.genweight
    
        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "muon"
            else: l1Flav = "electron"
            if event.lep2_isMu: l1Flav = "muon"
            else: l2Flav = "electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
        hBkgd.Fill(cuts["nocut"], genwt)
    
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
            hBkgd.Fill(cuts["dilepton"], genwt)
    
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hBkgd.Fill(cuts["deltaR(ll)>0.3"], genwt)
    
        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hBkgd.Fill(cuts["nbtag<2"], genwt)
    
        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue
            hBkgd.Fill(cuts["MET>80"], genwt)
    
        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hBkgd.Fill(cuts["no3rdlept"], genwt)
            
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue
            hBkgd.Fill(cuts["njets<4"], genwt)
    
    newProcess = False
    if not prevProcess == process:
        prevProcess = process
        processNum += 1
        newProcess = True

    c.cd()
    # hBkgd.Sumw2() # already summed while filling
    hBkgd.GetXaxis().SetTitle("cutflow")
    hBkgd.GetYaxis().SetTitle("Number of Events, norm to 3000 /fb")
    hBkgd.Scale(xsec*lumi/bkgdSubprocessGenweight)
    hBkgd.SetFillColor(processes[process])
    hBkgd.SetLineColor(processes[process])
    hBkgdStack.Add(hBkgd)
    if newProcess: legend.AddEntry(hBkgd, process+"_bkgd")

    for i, cut in enumerate(cuts):
        hBkgdCutsCountDict[process][i] += hBkgd.GetBinContent(i+1)

    bkgdFile.Close()

hBkgdStack.Draw("hist")
hBkgdStack.SetMinimum(1)
hBkgdStack.SetMaximum(10**12)

# show the number of events left over after each cut
processNum = 0
for process, color in processes.items():
    print
    print "Num surviving events after each cut from bkgd %s:" % process 
    for i, cut in enumerate(cuts):
        nEvtsLabel = TText()
        nEvtsLabel.SetNDC()
        nEvtsLabel.SetTextSize(0.02)
        nEvtsLabel.SetTextAlign(22)
        nEvtsLabel.SetTextAngle(0)
        nEvtsLabel.SetTextColor(color)
        nEvtsLabel.DrawText(0.1+0.4/nCuts+0.8*float(i)/nCuts, \
                0.7-(processNum)*0.02, \
                str(int(hBkgdCutsCountDict[process][i])))
        print cut, hBkgdCutsCountDict[process][i]
        nEvtsLabels.append(nEvtsLabel)
    processNum += 1
print

#--------------------------------------------------------------------------------#
# *************** Filling each signal in a separate hist  ************
print "Plotting from signal."

# assemble the sigsNtupleAdr
sigsNtupleAdr = baseDir+"stopCut_"
if numSigFiles < 10: sigsNtupleAdr += "0"+str(numSigFiles)
else: sigsNtupleAdr += str(numSigFiles)
sigsNtupleAdr += "Sig_"+channelName+".root"

sigDataListFile = open("sig_SingleStop_files")

coloropts = [2,4,3,6,7,9,28,46] # some good colors for lines
markeropts = [1,20,21,22,23] # some good marker styles for lines
linestyleopts = [1,2,3,7,9] # some good styles for lines

hSigArr = []
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

    hSig = TH1F("cutflow_sig_" + filename, "cutflow_sig_" + \
            filename[19:31], nCuts, 0, nCuts)
    hSig.SetDirectory(0)
    hSig.SetDefaultSumw2() # automatically sum w^2 while filling
    hSigArr.append(hSig)

    for i, cut in enumerate(cuts, start=1):
        if i>nCuts: break
        hSig.GetXaxis().SetBinLabel(i, cut)

    hSigGenweights = sigFile.Get("genweights")
    sigTotGenweight = hSigGenweights.GetSumOfWeights()

    # ********** Looping over events. ***********
    for count, event in enumerate(tSig):
        if count % 100000 == 0: print("count={0:d}".format(count))
        genwt = event.genweight

        if not findingSameFlavor: # if findingSameFlavor, l1/l2Flav set at runtime
            if event.lep1_isMu: l1Flav = "muon"
            else: l1Flav = "electron"
            if event.lep2_isMu: l1Flav = "muon"
            else: l2Flav = "electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
        hSig.Fill(cuts["nocut"], genwt)

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
            hSig.Fill(cuts["dilepton"], genwt)


        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hSig.Fill(cuts["deltaR(ll)>0.3"], genwt)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hSig.Fill(cuts["nbtag<2"], genwt)

        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue
            hSig.Fill(cuts["MET>80"], genwt)

        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hSig.Fill(cuts["no3rdlept"], genwt)
        
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue
            hSig.Fill(cuts["njets<4"], genwt)

    hcolor = coloropts[fileNum % len(coloropts)]
    hmarkerstyle = markeropts[(fileNum/len(coloropts)) % len(markeropts)]

    hSig.SetLineColor(hcolor) 
    hSig.SetMarkerStyle(hmarkerstyle)
    hSig.SetMarkerColor(hcolor)
    hlinestyle = linestyleopts[(fileNum/len(coloropts)/len(markeropts)) % \
            len(linestyleopts)]
    hSig.SetLineStyle(hlinestyle)

    hSig.Scale(xsec*lumi/sigTotGenweight)
    legend.AddEntry(hSig, hSig.GetTitle())
    hSig.Draw("hist same") # same pad

    # show the number of events left over after each cut
    print "Num surviving events after each cut from sig %s:" % filename 
    for i, cut in enumerate(cuts):
        print cut, hSig.GetBinContent(i+1)
        nEvtsLabel = TText()
        nEvtsLabel.SetNDC()
        nEvtsLabel.SetTextSize(0.02)
        nEvtsLabel.SetTextAlign(22)
        nEvtsLabel.SetTextAngle(0)
        nEvtsLabel.SetTextColor(hcolor)
        nEvtsLabel.DrawText(0.1+0.4/nCuts+0.8*float(i)/nCuts, \
                0.7-(processNum)*0.02-(1+fileNum)*0.02, \
                str(int(hSig.GetBinContent(i+1))))
        nEvtsLabels.append(nEvtsLabel)
    print

    sigFile.Close()

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print int(time.time()-start_time), "secs of processing."
print "Drawing."

legend.SetTextSize(0.02)
legend.Draw("same")
c.SetLogy()
c.Update()

if displayMode:
    print "Done. Press enter to finish."
    raw_input()
else:
    gSystem.ProcessEvents()
    imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
            "plots/v3CutSequence/cutflow_"+channelName+"_"+lastcut+".png"
    print "Saving image", imgName
    img = TImage.Create()
    img.FromPad(c)
    img.WriteImage(imgName)
    print "Done."

