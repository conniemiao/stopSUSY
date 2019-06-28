#!/usr/bin/env python

# NOTE: NEEDS >= 6 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, displayMode, findingSameFlavor, muPreference, lastcut, plotVar1, 
# plotVar2, ...
# True testMode plots only a few events; True displayMode displays rather than 
# saves w/o displaying the hists.
# Implements additional cuts and then draws 1D hist for data for some 
# variable(s), for the summed bkgd data and for each of the signal files. Also,
# draws the cutflow. 
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from sig_SingleStop_files
# Possible plotVars: listed in plotSettings, other than cutflow (which is always
# plotted)
# Possible lastcuts: listed in cuts below.

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TImage, TLegend, TText
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
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
        "nbtag":[5,0,5,""],
        "nbtagLoose":[10,0,10,""],
        "nbtagTight":[5,0,5,""],
        "dR_lep1_jet":[100,0,7,""],
        "dR_lep2_jet":[100,0,7,""],
        "met_pt":[100,0,500,"[GeV]"],
        "cutflow":[nCuts, 0, nCuts,""]
        }

for plotVar in plotVarArr:
    assert (plotVar in plotSettings), "invalid plotVar %s" % plotVar
plotVarArr.append("cutflow")

# Determining adr of bkgd and sig ntuples.
# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
displayMode = bool(int(sys.argv[2]))
print "Display mode:", displayMode
# applying cuts
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

print "Plotting",str(plotVarArr),"from",bkgdNtupleAdr,"and",sigsNtupleAdr
print "Cutting events up to and including", lastcut

numSigFiles = int(sigsNtupleAdr[48:50])

canvasDict = {}
legendDict = {}
nEvtsLabels = []

hBkgdDict = {} 
if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases
for plotVar in plotVarArr: # add an entry to the plotVar:hist dictionary
    nBins = plotSettings[plotVar][0]
    # if not testMode and nBins > 20: nBins = nBins * 5
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    binwidth = (xMax - xMin)/nBins
    hBkgd = TH1F(plotVar + "_bkgd", plotVar + "_bkgd", nBins, xMin, xMax)
    hBkgd.SetDirectory(0) # necessary to keep hist from closing
    hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
    hBkgdDict.update({plotVar:hBkgd})
    c = TCanvas("c_"+plotVar,"Plot",10,20,1000,700)
    canvasDict.update({plotVar:c})
    legendDict.update({plotVar:TLegend(.65,.75,.90,.90)})
lumi = 3000000 # luminosity = 3000 /fb = 3,000,000 /pb

gStyle.SetOptStat(0) # don't show any stats

start_time = time.time()
#--------------------------------------------------------------------------------#

# *************** Filling bkgd data summed together  ************
print
print "Plotting from background."
xsec = 67.75 # ttbar production cross section
bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
tBkgd = bkgdFile.Get("tBkgd")

nentries = tBkgd.GetEntries()
print("nentries={0:d}".format(nentries))
assert nentries > 0, "You have no events in your tree..."

hBkgdCutflow = hBkgdDict["cutflow"] 
for i, cut in enumerate(cuts, start=1):
    if i>nCuts: break
    hBkgdCutflow.GetXaxis().SetBinLabel(i, cut)

hBkgdGenweights = bkgdFile.Get("genweights")
bkgdTotGenweight = hBkgdGenweights.GetSumOfWeights()

# ********** Looping over events. ***********
for count, event in enumerate(tBkgd):
    if count % 100000 == 0: print("count={0:d}".format(count))
    genwt = event.genweight

    if not findingSameFlavor: # if findingSameFlavor, l1/l2Flav set at runtime
        if event.lep1_isMu: l1Flav = "muon"
        else: l1Flav = "electron"
        if event.lep2_isMu: l1Flav = "muon"
        else: l2Flav = "electron"
    l1Index = event.lep1_index
    l2Index = event.lep2_index
    hBkgdCutflow.Fill(cuts["nocut"], genwt)

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
        hBkgdCutflow.Fill(cuts["dilepton"], genwt)

    # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
    # hBkgdCutflow.Fill(cuts["deltaR(ll)>0.3"], genwt)

    if nCuts > cuts["nbtag<2"]:
        if event.nbtag > 1: continue
        hBkgdCutflow.Fill(cuts["nbtag<2"], genwt)

    if nCuts > cuts["MET>80"]:
        if event.met_pt < 80: continue
        hBkgdCutflow.Fill(cuts["MET>80"], genwt)

    if nCuts > cuts["no3rdlept"]:
        if event.found3rdLept: continue
        hBkgdCutflow.Fill(cuts["no3rdlept"], genwt)
        
    if nCuts > cuts["njets<4"]:
        if event.njets >= 4: continue
        hBkgdCutflow.Fill(cuts["njets<4"], genwt)

    # ********** Plotting. ***********
    if event.njets > 0:
        jMaxPt = 0
        for j in range(event.njets):
            if np.reshape(event.jet_pt,20)[j] > \
                    np.reshape(event.jet_pt,20)[jMaxPt]: jMaxPt = j
        dR_lep1_jet = deltaR(event, l1Flav, l1Index, "jet", jMaxPt)
        dR_lep2_jet = deltaR(event, l2Flav, l2Index, "jet", jMaxPt)
    for plotVar in plotVarArr:
        if plotVar == "cutflow": break

        hBkgd = hBkgdDict[plotVar]
        xMin = plotSettings[plotVar][1]
        xMax = plotSettings[plotVar][2]
        if plotVar[:4] == "lep1": 
            val = np.reshape(getattr(event, l1Flav+plotVar[4:]),20)[l1Index]
        elif plotVar[:4] == "lep2": 
            val = np.reshape(getattr(event, l2Flav+plotVar[4:]),20)[l2Index]
        elif plotVar[:3] == "jet":
            if event.njets == 0: continue # to the next plotVar
            val = np.reshape(getattr(event, "jet_"+plotVar[4:]),20)[jMaxPt]
        elif plotVar == "dR_lep1_jet": 
            if event.njets == 0: continue # to the next plotVar
            val = dR_lep1_jet
        elif plotVar == "dR_lep2_jet": 
            if event.njets == 0: continue # to the next plotVar
            val = dR_lep2_jet
        else: val = getattr(event, plotVar)

        if val <= xMax:
            hBkgd.Fill(val, genwt)
        else: # overflow
            hBkgd.Fill(xMax - binwidth/2, genwt)

for plotVar in plotVarArr:
    c = canvasDict[plotVar]
    c.cd()
    hBkgd = hBkgdDict[plotVar]
    # hBkgd.Sumw2() # already summed while filling
    title = plotVar+" ("+channelName+", cuts to "+lastcut+")"
    hBkgd.SetTitle(title)
    unitsLabel = plotSettings[plotVar][3]
    hBkgd.GetXaxis().SetTitle(plotVar+" "+unitsLabel)
    hBkgd.GetYaxis().SetTitle("Number of Events, norm to 3000 /fb")
    hBkgd.Scale(xsec*lumi/bkgdTotGenweight)
    hBkgd.SetMinimum(1)
    hBkgd.SetMaximum(10**12)
    hBkgd.SetLineColor(1) # black
    legend = legendDict[plotVar]
    legend.AddEntry(hBkgd, plotVar+"_bkgd")
    hBkgd.Draw("hist")
print "Num surviving events after each cut from bkgd:" 
for i, cut in enumerate(cuts):
    print cut, hBkgdCutflow.GetBinContent(i+1)
    nEvtsLabel = TText()
    nEvtsLabel.SetNDC()
    nEvtsLabel.SetTextSize(0.02)
    nEvtsLabel.SetTextAlign(22)
    nEvtsLabel.SetTextAngle(0)
    nEvtsLabel.SetTextColor(1)
    nEvtsLabel.DrawText(0.1+0.4/nCuts+0.8*float(i)/nCuts, 0.1+(numSigFiles+1)*\
            0.02, str(int(hBkgdCutflow.GetBinContent(i+1))))
    nEvtsLabels.append(nEvtsLabel)
print
bkgdFile.Close()

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting from signal."
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
        # if not testMode and nBins > 20: nBins = nBins * 5
        xMin = plotSettings[plotVar][1]
        xMax = plotSettings[plotVar][2]
        binwidth = (xMax - xMin)/nBins
        hSigArr = hSigArrDict[plotVar]  # one hist for each signal file
        hSig = TH1F(plotVar + "_sig_" + filename, plotVar + "_sig_" + \
                filename[21:31], nBins, xMin, xMax)
        hSig.SetDirectory(0)
        hSig.SetDefaultSumw2() # automatically sum w^2 while filling
        hSigArr.append(hSig)

    hSigCutflow = hSigArrDict["cutflow"][fileNum]
    for i, cut in enumerate(cuts, start=1):
        if i>nCuts: break
        hSigCutflow.GetXaxis().SetBinLabel(i, cut)

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
        hSigCutflow.Fill(cuts["nocut"], genwt)

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
            hSigCutflow.Fill(cuts["dilepton"], genwt)


        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hSigCutflow.Fill(cuts["deltaR(ll)>0.3"], genwt)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hSigCutflow.Fill(cuts["nbtag<2"], genwt)

        if nCuts > cuts["MET>80"]:
            if event.met_pt < 80: continue
            hSigCutflow.Fill(cuts["MET>80"], genwt)

        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hSigCutflow.Fill(cuts["no3rdlept"], genwt)
        
        if nCuts > cuts["njets<4"]:
            if event.njets >= 4: continue
            hSigCutflow.Fill(cuts["njets<4"], genwt)

        # ********** Plotting. ***********
        if event.njets > 0:
            jMaxPt = 0
            for j in range(event.njets):
                if np.reshape(event.jet_pt,20)[j] > \
                        np.reshape(event.jet_pt,20)[jMaxPt]: jMaxPt = j
            dR_lep1_jet = deltaR(event, l1Flav, l1Index, "jet", jMaxPt)
            dR_lep2_jet = deltaR(event, l2Flav, l2Index, "jet", jMaxPt)

        for plotVar in plotVarArr:

            if plotVar == "cutflow": break

            hSig = hSigArrDict[plotVar][fileNum]
            xMin = plotSettings[plotVar][1]
            xMax = plotSettings[plotVar][2]
            if plotVar[:4] == "lep1": 
                val = np.reshape(getattr(event, l1Flav+plotVar[4:]),20)[l1Index]
            elif plotVar[:4] == "lep2": 
                val = np.reshape(getattr(event, l2Flav+plotVar[4:]),20)[l2Index]
            elif plotVar[:3] == "jet":
                if event.njets == 0: continue # to the next plotVar
                val = np.reshape(getattr(event, "jet_"+plotVar[4:]),20)[jMaxPt]
            elif plotVar == "dR_lep1_jet": 
                if event.njets == 0: continue # to the next plotVar
                val = dR_lep1_jet
            elif plotVar == "dR_lep2_jet": 
                if event.njets == 0: continue # to the next plotVar
                val = dR_lep2_jet
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

        # print "here3"
        # sys.stdout.flush()
        # hSig.Sumw2()
        # print "here4"
        # sys.stdout.flush()
        hSig.Scale(xsec*lumi/sigTotGenweight)
        hSig.SetMinimum(1)
        hSig.SetMaximum(10**12)
        legend = legendDict[plotVar]
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad

    print "Num surviving events after each cut from sig %s:" % filename 
    for i, cut in enumerate(cuts):
        print cut, hSigCutflow.GetBinContent(i+1)
        nEvtsLabel = TText()
        nEvtsLabel.SetNDC()
        nEvtsLabel.SetTextSize(0.02)
        nEvtsLabel.SetTextAlign(22)
        nEvtsLabel.SetTextAngle(0)
        nEvtsLabel.SetTextColor(hcolor)
        nEvtsLabel.DrawText(0.1+0.4/nCuts+0.8*float(i)/nCuts, \
                0.1+(1+fileNum)*0.02, str(int(hSigCutflow.GetBinContent(i+1))))
        nEvtsLabels.append(nEvtsLabel)
    print

    sigFile.Close()


#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print int(time.time()-start_time), "secs of processing."

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
                "plots/v3CutSequence/"+plotVar+"_"+channelName+"_"+lastcut+".png"
        print "Saving image", imgName
        img = TImage.Create()
        img.FromPad(canvasDict[plotVar])
        img.WriteImage(imgName)
    print "Done."

