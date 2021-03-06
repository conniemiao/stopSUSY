#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, experimental, findingSameFlavor, muPreference
# True testMode plots only a few events. True experimental mode makes extra cuts
# for cut optimization.
# Implements additional cuts and then writes the counts of each sig/bkgd type
# after each cut to an output file of the form cutflow_stats_[channel].txt in
# the plots folder.
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from bkgd_files
# Uses xsec info from sig_SingleStop_files

print "Importing modules."
import sys
from ROOT import TFile, TTree, TH1D, TCanvas, TImage, TLegend, TText, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
import numpy as np
import time
print "Beginning execution of", sys.argv

assert len(sys.argv) == 5, "needs 4 command line args: testMode{0,1}, experimental{0,1}, findingSameFlavor{0,1}, muPreference{0,1}"

# Determining adr of bkgd and sig ntuples.
# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
experimental = bool(int(sys.argv[2]))
print "Experimental mode:", experimental
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

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
        ("MET>80",4),("njets<4",5)])
if experimental:
    cuts = OrderedDict([("nocut",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3),\
            ("MET>20",4), ("MET>50",5), ("MET>80",6), ("MET>110",7), ("njets<4",8)])
nCuts = len(cuts)


# bkgd process name : color for plotting
processes = OrderedDict([("W-Jets",38), ("Drell-Yan",46), ("Diboson",41), \
        ("Single-Top",30), ("TT+X",7)])

baseDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/"
# number of files to process
numBkgdFiles = float("inf")  # note: must loop over all files to have correct xsec
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
        hBkgd = TH1D("cutflow_"+subprocess+"_bkgd", \
                "cutflow_"+subprocess+"_bkgd", nCuts, 0, nCuts)
        hBkgd.SetDirectory(0) # necessary to keep hist from closing
        hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
        hBkgdDict.update({subprocess:hBkgd})
for process in processes:
    hBkgdCutsCountDict.update({process:[0]*nCuts})
title = "cutflow ("+channelName+")"
hBkgdStack = THStack("cutflow_bkgdStack", title)

lumi = 3000000 # luminosity = 3 /ab = 3000 /fb = 3,000,000 /fb

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
    bkgdNtupleAdr += "all_"
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
    
    nMax = nentries
    if testMode: nMax = 10000

    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        if count > nMax : break
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
        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if findingSameFlavor:
                if event.lep1_relIso >= 0.1: continue
                if event.lep2_relIso >= 0.1: continue
            else:
                if event.lep1_relIso >= 0.2: continue
                if event.lep2_relIso >= 0.2: continue
            hBkgd.Fill(cuts["dilepton"], genwt)
    
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hBkgd.Fill(cuts["deltaR(ll)>0.3"], genwt)
    
        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hBkgd.Fill(cuts["no3rdlept"], genwt)
    
        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hBkgd.Fill(cuts["nbtag<2"], genwt)

        if experimental:
            if nCuts > cuts["MET>20"]:
                if event.met_pt < 20: continue
                hBkgd.Fill(cuts["MET>20"], genwt)
            if nCuts > cuts["MET>50"]:
                if event.met_pt < 50: continue
                hBkgd.Fill(cuts["MET>50"], genwt)

        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue
            hBkgd.Fill(cuts["MET>80"], genwt)

        if experimental:
            if nCuts > cuts["MET>110"]:
                if event.met_pt < 110: continue
                hBkgd.Fill(cuts["MET>110"], genwt)
            
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue
            hBkgd.Fill(cuts["njets<4"], genwt)
    
    # hBkgd.Sumw2() # already summed while filling
    hBkgd.Scale(xsec*lumi/bkgdSubprocessGenweight)

    for i, cut in enumerate(cuts):
        hBkgdCutsCountDict[process][i] += int(hBkgd.GetBinContent(i+1))

    bkgdFile.Close()

#--------------------------------------------------------------------------------#
# *************** Filling each signal in a separate hist  ************
print
print "Plotting from signal."
print

# assemble the sigsNtupleAdr
sigsNtupleAdr = baseDir+"stopCut_"
if numSigFiles < 10: sigsNtupleAdr += "0"+str(numSigFiles)
else: sigsNtupleAdr += str(numSigFiles)
sigsNtupleAdr += "Sig_"+channelName+".root"

sigDataListFile = open("sig_SingleStop_files")

hSigArr = []

# technically an array of arrays, but can use as a dictionary mapping fileNum 
# for a sig to an array containing num evts remaining after cut i:
hSigCutsCountDict = [] 

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

    hSig = TH1D("sig_" + filename, "sig_" + filename[18:31], nCuts, 0, nCuts)
    hSig.SetDirectory(0)
    hSig.SetDefaultSumw2() # automatically sum w^2 while filling
    hSigArr.append(hSig)
    hSigCutsCountDict.append([0]*nCuts)

    for i, cut in enumerate(cuts, start=1):
        if i>nCuts: break
        hSig.GetXaxis().SetBinLabel(i, cut)

    hSigGenweights = sigFile.Get("genweights"+str(fileNum))
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
        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if findingSameFlavor:
                if event.lep1_relIso >= 0.1: continue
                if event.lep2_relIso >= 0.1: continue
            else:
                if event.lep1_relIso >= 0.2: continue
                if event.lep2_relIso >= 0.2: continue
            hSig.Fill(cuts["dilepton"], genwt)


        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hSig.Fill(cuts["deltaR(ll)>0.3"], genwt)

        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hSig.Fill(cuts["no3rdlept"], genwt)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hSig.Fill(cuts["nbtag<2"], genwt)

        if experimental:
            if nCuts > cuts["MET>20"]:
                if event.met_pt < 20: continue
                hSig.Fill(cuts["MET>20"], genwt)
            if nCuts > cuts["MET>50"]:
                if event.met_pt < 50: continue
                hSig.Fill(cuts["MET>50"], genwt)

        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue
            hSig.Fill(cuts["MET>80"], genwt)

        if experimental:
            if nCuts > cuts["MET>110"]:
                if event.met_pt < 110: continue
                hSig.Fill(cuts["MET>110"], genwt)
        
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue
            hSig.Fill(cuts["njets<4"], genwt)
    hSig.Scale(xsec*lumi/sigTotGenweight)

    for i, cut in enumerate(cuts):
        hSigCutsCountDict[fileNum][i] = int(hSig.GetBinContent(i+1))

    sigFile.Close()

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print int(time.time()-start_time), "secs of processing."

statsHeader = "Cut_name      "
statsStack = np.array(cuts.keys()).reshape(nCuts, 1)
for process in processes:
    statsHeader +=  process + ''.join([' ']*(16-len(process)))
    statsStack = np.append(statsStack, \
            np.array(hBkgdCutsCountDict[process]).reshape(nCuts, 1), axis=1)
for fileNum in range(numSigFiles):
    sigName = hSigArr[fileNum].GetTitle()[4:]
    statsHeader += sigName + ''.join([' ']*(16-len(sigName)))
    statsStack = np.append(statsStack, \
            np.array(hSigCutsCountDict[fileNum]).reshape(nCuts, 1), axis=1)
statsFileName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
        "plots/HLLHC/v3CutSequence/cutflow_stats/cutflow_stats_"+channelName
if experimental: statsFileName += "_experimental"
statsFileName += ".txt"
np.savetxt(statsFileName, statsStack, delimiter='   ', header=statsHeader, \
        fmt='%-13s')
print "Saved file", statsFileName

print "Done."

