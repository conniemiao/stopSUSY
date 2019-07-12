#!/usr/bin/env python

# NOTE: NEEDS >= 6 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, displayMode, findingSameFlavor, muPreference, lastcut, plotVar1, 
# plotVar2, ...
# True testMode plots only a few events; True displayMode displays rather than 
# saves w/o displaying the hists.
# Implements additional cuts and then draws 1D hist for data for some 
# variable(s), for the summed bkgd data and for each of the signal files.
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from bkgd_files
# Uses xsec info from sig_SingleStop_files
# Possible plotVars: listed in plotSettings
# Possible lastcuts: listed in cuts below.

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TImage, TLegend, TText, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
from math import sqrt, cos
import numpy as np
import time

assert len(sys.argv) >= 7, "need at least 6 command line args: testMode{0,1}, displayMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}, lastcut, plotVar1, plotVar2, ..."

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("nbtag<2",2), ("MET>80",3),\
        ("no3rdlept",4), ("njets<4",5)])
lastcut = sys.argv[5]
assert lastcut in cuts, "invalid last cut %s" % lastcut
nCuts = cuts[lastcut]+1

plotVarArr = sys.argv[6:]

plotSettings = { # [nBins,xMin,xMax,units]
        "lep1_pt":[100,0,400,"[Gev]"],
        "lep1_eta":[100,-4,4,""],
        "lep1_phi":[100,-4,4,""],
        "lep1_relIso":[100,0,0.2,""],
        "lep1_mt":[100,0,500,"[GeV]"],
        "lep2_pt":[100,0,400,"[GeV]"],
        "lep2_eta":[100,-4,4,""],
        "lep2_phi":[100,-4,4,""],
        "lep2_relIso":[100,0,0.2,""],
        "lep2_mt":[100,0,500,"[GeV]"],
        "njets":[10,0,10,""],
        "jet_pt":[100,0,400,"[GeV]"], 
        "jet_eta":[100,-3,3,""],
        "jet_phi":[100,-4,4,""],
        "jet_ht":[100,0,800,"[GeV]"],
        "nbtag":[5,0,5,""],
        "nbtagLoose":[10,0,10,""],
        "nbtagTight":[5,0,5,""],
        "dR_lep1_jet":[100,0,7,""],
        "dR_lep2_jet":[100,0,7,""],
        "met_pt":[100,0,500,"[GeV]"],
        "mt_tot":[100,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        "mt_sum":[100,0,1000,"[GeV]"], # mt1 + mt2
        "m_eff":[100,0,1000,"[GeV]"], # ht + met + pt1 + pt2
        "jet_ht_div_sqrt_met":[100,0,100,""],
        "mt_tot_div_sqrt_met":[100,0,100,""],
        "m_eff_div_sqrt_met":[100,0,100,""]
        }

for plotVar in plotVarArr:
    assert (plotVar in plotSettings), "invalid plotVar %s" % plotVar
print "Plotting", str(plotVarArr)

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

print "Cutting events up to and including", lastcut

# bkgd process name : color for plotting
processes = {"W-Jets":38, "Drell-Yan":46, "Diboson":41, "Single-Top":30, \
        "TT+X":7}
if testMode:
    processes = {"Diboson":41, "TT+X":7}


canvasDict = {}
legendDict = {}
hBkgdStacksDict = {} # maps plotVar to the stack of background
nEvtsLabels = []

if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases

baseDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/"
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

# For each plotVar in hBkgdPlotVarSubprocessesDict, there is a dictionary 
# which maps every subprocess to an hBkgd which contains data from all the 
# ntuples for that subprocess.
hBkgdPlotVarSubprocessesDict = {} 
for plotVar in plotVarArr: # add an entry to the plotVar:hist dictionary
    nBins = plotSettings[plotVar][0]
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    binwidth = (xMax - xMin)/nBins
    hBkgdPlotVarSubprocessesDict.update({plotVar:{}})
    with open("bkgd_files") as bkgdSubprocessesListFile:
        for subprocessLine in bkgdSubprocessesListFile:
            subprocessLine = subprocessLine.rstrip('\n')
            subprocess, process, xsec = subprocessLine.split(" ")
            if subprocess[0] == "#": continue # problematic input files
            hBkgd = TH1F(subprocess+"_bkgd", \
                    subprocess+"_bkgd", nBins, xMin, xMax)
            hBkgd.SetDirectory(0) # necessary to keep hist from closing
            hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
            hBkgdPlotVarSubprocessesDict[plotVar].update({subprocess:hBkgd})
    c = TCanvas("c_"+plotVar,"Plot",10,20,1000,700)
    canvasDict.update({plotVar:c})
    legendDict.update({plotVar:TLegend(.70,.70,.90,.90)})
    title = plotVar+" ("+channelName+", cuts to "+lastcut+")"
    hBkgdStacksDict.update({plotVar:THStack(plotVar+"_bkgdStack", title)})
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

    if subprocess[:1] == "#": continue # problematic input files
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

    # nMax = 10000
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        # if count > nMax : break
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
        if event.njets > 0:
            jMaxPt = 0
            for j in range(event.njets):
                if np.reshape(event.jet_pt,20)[j] > \
                        np.reshape(event.jet_pt,20)[jMaxPt]: jMaxPt = j
            dR_lep1_jet = deltaR(event, l1Flav, l1Index, "jet", jMaxPt)
            dR_lep2_jet = deltaR(event, l2Flav, l2Index, "jet", jMaxPt)

        for plotVar in plotVarArr:
            hBkgd = hBkgdPlotVarSubprocessesDict[plotVar][subprocess]
            nBins = plotSettings[plotVar][0]
            xMin = plotSettings[plotVar][1]
            xMax = plotSettings[plotVar][2]
            binwidth = (xMax - xMin)/nBins
            if plotVar[:4] == "lep1": 
                val = np.reshape(getattr(event, l1Flav+plotVar[4:]),\
                        20)[l1Index]
            elif plotVar[:4] == "lep2": 
                val = np.reshape(getattr(event, l2Flav+plotVar[4:]),\
                        20)[l2Index]
            elif plotVar[:6] == "jet_ht":
                if event.njets == 0: continue # to the next plotVar
                val = event.jet_ht
                if plotVar == "jet_ht_div_sqrt_met":
                    val /= sqrt(event.met_pt)
            elif plotVar[:3] == "jet":
                if event.njets == 0: continue # to the next plotVar
                val = np.reshape(getattr(event, "jet_"+plotVar[4:]),\
                        20)[jMaxPt]
            elif plotVar == "dR_lep1_jet": 
                if event.njets == 0: continue # to the next plotVar
                val = dR_lep1_jet
            elif plotVar == "dR_lep2_jet": 
                if event.njets == 0: continue # to the next plotVar
                val = dR_lep2_jet
            elif plotVar[:6] == "mt_tot":
                val = sqrt((np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index])**2 + (np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index])**2)
                if plotVar == "mt_tot_div_sqrt_met":
                    val /= sqrt(event.met_pt)
            elif plotVar == "mt_sum":
                val = np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index]
            elif plotVar[:5] == "m_eff":
                val = event.met_pt + np.reshape(getattr(event, \
                        l1Flav+"_pt"),20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_pt"),20)[l2Index]
                if event.njets > 0: val += event.jet_ht
                if plotVar == "m_eff_div_sqrt_met":
                    val /= sqrt(event.met_pt)
            else: val = getattr(event, plotVar)
    
            if val <= xMax:
                hBkgd.Fill(val, genwt)
            else: # overflow
                hBkgd.Fill(xMax - binwidth/2, genwt)
    
    newProcess = False
    if not prevProcess == process:
        prevProcess = process
        processNum += 1
        newProcess = True

    for plotVar in plotVarArr:
        c = canvasDict[plotVar]
        c.cd()
        hBkgd = hBkgdPlotVarSubprocessesDict[plotVar][subprocess]
        # hBkgd.Sumw2() # already summed while filling
        hBkgd.Scale(xsec*lumi/bkgdTotGenweight)
        hBkgd.SetFillColor(processes[process])
        hBkgd.SetLineColor(processes[process])
        hBkgdStacksDict[plotVar].Add(hBkgd)
        if newProcess:
            legend = legendDict[plotVar]
            legend.AddEntry(hBkgd, process+"_bkgd")
    print
    bkgdFile.Close()

for plotVar in plotVarArr:
    c = canvasDict[plotVar]
    c.cd()
    hBkgdStack = hBkgdStacksDict[plotVar]
    hBkgdStack.Draw("hist")
    unitsLabel = plotSettings[plotVar][3]
    hBkgdStack.GetXaxis().SetTitle(plotVar+" "+unitsLabel)
    hBkgdStack.GetYaxis().SetTitle("Number of Events, norm to 3000 /fb")
    hBkgdStack.SetMinimum(1)
    hBkgdStack.SetMaximum(10**12)

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

hSigArrDict = {}
for plotVar in plotVarArr: # add an entry to the plotVar:hist dictionary
    hSigArrDict.update({plotVar:[]})

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

    for plotVar in plotVarArr:
        nBins = plotSettings[plotVar][0]
        xMin = plotSettings[plotVar][1]
        xMax = plotSettings[plotVar][2]
        hSigArr = hSigArrDict[plotVar]  # one hist for each signal file
        hSig = TH1F("sig_" + filename, "sig_" + \
                filename[18:31], nBins, xMin, xMax)
        hSig.SetDirectory(0)
        hSig.SetDefaultSumw2() # automatically sum w^2 while filling
        hSigArr.append(hSig)

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
        if event.njets > 0:
            jMaxPt = 0
            for j in range(event.njets):
                if np.reshape(event.jet_pt,20)[j] > \
                        np.reshape(event.jet_pt,20)[jMaxPt]: jMaxPt = j
            dR_lep1_jet = deltaR(event, l1Flav, l1Index, "jet", jMaxPt)
            dR_lep2_jet = deltaR(event, l2Flav, l2Index, "jet", jMaxPt)

        for plotVar in plotVarArr:
            hSig = hSigArrDict[plotVar][fileNum]
            nBins = plotSettings[plotVar][0]
            xMin = plotSettings[plotVar][1]
            xMax = plotSettings[plotVar][2]
            binwidth = (xMax - xMin)/nBins
            if plotVar[:4] == "lep1": 
                val = np.reshape(getattr(event, l1Flav+plotVar[4:]),\
                        20)[l1Index]
            elif plotVar[:4] == "lep2": 
                val = np.reshape(getattr(event, l2Flav+plotVar[4:]),\
                        20)[l2Index]
            elif plotVar[:6] == "jet_ht":
                if event.njets == 0: continue # to the next plotVar
                val = event.jet_ht
                if plotVar == "jet_ht_div_sqrt_met":
                    val /= sqrt(event.met_pt)
            elif plotVar[:3] == "jet":
                if event.njets == 0: continue # to the next plotVar
                val = np.reshape(getattr(event, "jet_"+plotVar[4:]),\
                        20)[jMaxPt]
            elif plotVar == "dR_lep1_jet": 
                if event.njets == 0: continue # to the next plotVar
                val = dR_lep1_jet
            elif plotVar == "dR_lep2_jet": 
                if event.njets == 0: continue # to the next plotVar
                val = dR_lep2_jet
            elif plotVar[:6] == "mt_tot":
                val = sqrt((np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index])**2 + (np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index])**2)
                if plotVar == "mt_tot_div_sqrt_met":
                    val /= sqrt(event.met_pt)
            elif plotVar == "mt_sum":
                val = np.reshape(getattr(event, l1Flav+"_mt"),\
                        20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_mt"),20)[l2Index]
            elif plotVar[:5] == "m_eff":
                val = event.met_pt + np.reshape(getattr(event, \
                        l1Flav+"_pt"),20)[l1Index] + np.reshape(getattr(event, \
                        l2Flav+"_pt"),20)[l2Index]
                if event.njets > 0: val += event.jet_ht
                if plotVar == "m_eff_div_sqrt_met":
                    val /= sqrt(event.met_pt)
            else: val = getattr(event, plotVar)

            if val <= xMax:
                hSig.Fill(val, genwt)
            else: # overflow
                hSig.Fill(xMax - binwidth/2, genwt)


    hcolor = coloropts[fileNum % len(coloropts)]
    hmarkerstyle = markeropts[(fileNum/len(coloropts)) % len(markeropts)]

    for plotVar in plotVarArr:
        c = canvasDict[plotVar]
        c.cd()
        hSig = hSigArrDict[plotVar][fileNum]
        hSig.SetLineColor(hcolor) 
        hSig.SetMarkerStyle(hmarkerstyle)
        hSig.SetMarkerColor(hcolor)
        hlinestyle = linestyleopts[(fileNum/len(coloropts)/len(markeropts)) % \
                len(linestyleopts)]
        hSig.SetLineStyle(hlinestyle)

        hSig.Scale(xsec*lumi/sigTotGenweight)
        legend = legendDict[plotVar]
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad
    print
    sigFile.Close()

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print int(time.time()-start_time), "secs of processing."
print "Drawing."

for plotVar in plotVarArr:
    c = canvasDict[plotVar]
    c.cd()
    legend = legendDict[plotVar]
    legend.SetTextSize(0.02)
    legend.Draw("same")
    c.SetLogy()
    c.Update()

if displayMode:
    print "Done. Press enter to finish."
    raw_input()
else:
    gSystem.ProcessEvents()
    for plotVar in plotVarArr:
        imgName = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
                "plots/v3CutSequence/plot1D/"+plotVar+"_"+channelName+"_"+\
                lastcut+".png"
        print "Saving image", imgName
        img = TImage.Create()
        img.FromPad(canvasDict[plotVar])
        img.WriteImage(imgName)
    print "Done."

