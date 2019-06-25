#!/usr/bin/env python

# NOTE: NEEDS >= 5 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, displayMode, findingSameFlavor, muPreference, plotVar1, plotVar2, ...
# True testMode plots only a few events; True displayMode displays rather than 
# saves w/o displaying the hists.
# Implements additional cuts and then draws 1D hist for data for some 
# variable(s), for the summed bkgd data and for each of the signal files. Also,
# draws the cutflow. 
# Uses the root files outputted by makeNtupleBkgd.py and makeNtupleSigs.py
# Uses xsec info from sig_SingleStop_files
# Possible plotVars: listed in plotSettings, other than cutflow (which is always
# plotted)

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TImage, TLegend
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
import numpy as np


assert len(sys.argv) >= 6, "need at least 5 command line args: testMode{0,1}, displayMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}, plotVar1, plotVar2, ..."
plotVarArr = sys.argv[5:]

cuts = OrderedDict([("no cut",0), ("dilepton",1), ("nbtag<2",2), ("MET>80",3),\
        ("no 3rd lepton",4), ("njets<4",5)])
nCuts = len(cuts)

plotSettings = { # [nBins,xMin,xMax]
        "lep1_pt":[100,0,400],
        "lep1_eta":[100,-4,4],
        "lep1_phi":[100,-4,4],
        "lep1_relIso":[100,0,0.2],
        "lep1_mt":[100,0,500],
        "lep2_pt":[100,0,400],
        "lep2_eta":[100,-4,4],
        "lep2_phi":[100,-4,4],
        "lep2_relIso":[100,0,0.2],
        "lep2_mt":[100,0,500],
        "njets":[10,0,10],
        "jet_pt":[100,0,400], 
        "jet_eta":[100,-3,3],
        "jet_phi":[100,-4,4],
        "nbtag":[5,0,5],
        "nbtagLoose":[10,0,10],
        "nbtagTight":[5,0,5],
        "dR_lep1_jet":[100,0,7],
        "dR_lep2_jet":[100,0,7],
        "met_pt":[100,0,500],
        "cutflow":[nCuts, 0, nCuts]
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

print "Plotting",str(plotVarArr),"from",bkgdNtupleAdr,"and",sigsNtupleAdr

numSigFiles = int(sigsNtupleAdr[48:50])

canvasDict = {}
canvasDict.update({"cutflow":TCanvas("c_cutflow","Plot",10,20,1000,700)})
legendDict = {}
legendDict.update({"cutflow":TLegend(.70,.75,.90,.90)})

hBkgdDict = {} 
if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases
for plotVar in plotVarArr: # add an entry to the plotVar:hist dictionary
    nBins = plotSettings[plotVar][0]
    # if not testMode and nBins > 20: nBins = nBins * 5
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    binwidth = (xMax - xMin)/nBins
    hBkgd = TH1F(plotVar + "_bkgd", plotVar + "_bkgd", \
            nBins + 1, xMin, xMax)
    hBkgd.SetDirectory(0) # necessary to keep hist from closing
    hBkgdDict.update({plotVar:hBkgd})
    c = TCanvas("c_"+plotVar,"Plot",10,20,1000,700)
    canvasDict.update({plotVar:c})
    legendDict.update({plotVar:TLegend(.70,.75,.90,.90)})
lumi = 3000000 # luminosity = 3000 /fb = 3,000,000 /pb

gStyle.SetOptStat(0) # don't show any stats

#--------------------------------------------------------------------------------#

# *************** Filling bkgd data summed together  ************
print "Plotting from background."
xsec = 67.75 # ttbar production cross section
bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
tBkgd = bkgdFile.Get("tBkgd")

nentries = tBkgd.GetEntries()
print("nentries={0:d}".format(nentries))
assert nentries > 0, "You have no events in your tree..."

hBkgdCutflow = hBkgdDict["cutflow"] 
for i, cut in enumerate(cuts, start=1):
    hBkgdCutflow.GetXaxis().SetBinLabel(i, cut)

for count, event in enumerate(tBkgd):
    if count % 500000 == 0: print("count={0:d}".format(count))
    hBkgdCutflow.Fill(cuts["no cut"])

    # ********** Additional cuts. ***********
    if findingSameFlavor:
        if muPreference:
            lepIndices = selectMuMu(event)
        else: lepIndices = selectElEl(event)
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
    hBkgdCutflow.Fill(cuts["dilepton"])

    # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
    # hBkgdCutflow.Fill(cuts["deltaR(ll)>0.3"])

    if event.nbtag > 1: continue
    hBkgdCutflow.Fill(cuts["nbtag<2"])

    if event.met_pt < 80: continue
    hBkgdCutflow.Fill(cuts["MET>80"])

    # veto (3rd lepton) checks:
    if findingSameFlavor:
        # event should not give valid muel or elmu pair
        if selectMuEl(event) is not None: continue
        if selectElMu(event) is not None: continue
    else:
        # event should not give valid mumu or elel pair
        if selectMuMu(event) is not None: continue
        if selectElEl(event) is not None: continue
    hBkgdCutflow.Fill(cuts["no 3rd lepton"])
    
    if event.njets >= 4: continue
    hBkgdCutflow.Fill(cuts["njets<4"])

    jMaxPt = 0
    for j in range(event.njets):
        if np.reshape(event.jet_pt,20)[j] > np.reshape(event.jet_pt,20)[jMaxPt]:
            jMaxPt = j
    dR_lep1_jet = deltaR(event, l1Flav, l1Index, "jet", jMaxPt)
    dR_lep2_jet = deltaR(event, l2Flav, l2Index, "jet", jMaxPt)

    # ********** Plotting. ***********
    for plotVar in plotVarArr:
        if plotVar == "cutflow": break
        hBkgd = hBkgdDict[plotVar]
        if plotVar[:4] == "lep1": 
            val = np.reshape(getattr(event, l1Flav+plotVar[4:]),20)[l1Index]
        elif plotVar[:4] == "lep2": 
            val = np.reshape(getattr(event, l2Flav+plotVar[4:]),20)[l2Index]
        elif plotVar[:3] == "jet":
            val = np.reshape(getattr(event, "jet"+plotVar[4:]),20)[jMaxPt]
        elif plotVar == "dR_lep1_jet": val = dR_lep1_jet
        elif plotVar == "dR_lep2_jet": val = dR_lep2_jet
        else: val = getattr(event, plotVar)
        if val <= xMax:
            hBkgd.Fill(val, 1)
        else: # overflow
            hBkgd.Fill(xMax - binwidth/2, 1)
for plotVar in plotVarArr:
    c = canvasDict[plotVar]
    c.cd()
    hBkgd = hBkgdDict[plotVar]
    hBkgd.Sumw2()
    title = plotVar + " ("+bkgdNtupleAdr[-18:-14]+")"
    title += ", "+bkgdNtupleAdr[-13:-5]
    hBkgd.SetTitle(title)
    hBkgd.GetXaxis().SetTitle(plotVar+" [GeV]")
    hBkgd.GetYaxis().SetTitle("Number of Events, norm to 3000 /fb")
    hBkgd.Scale(xsec*lumi/hBkgd.GetSumOfWeights())
    hBkgd.SetMinimum(1)
    hBkgd.SetMaximum(10**12)
    hBkgd.SetLineColor(1) # black
    legend = legendDict[plotVar]
    legend.AddEntry(hBkgd, plotVar+"_bkgd")
    hBkgd.Draw("hist")
print "Num surviving events after each cut from bkgd:" 
for i in range(0,nCuts):
    print hBkgdCutflow.GetXaxis().GetBinLabel(i+1),hBkgdCutflow.GetBinContent(i+1)
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
        if plotVar == "cutflow": break
        nBins = plotSettings[plotVar][0]
        # if not testMode and nBins > 20: nBins = nBins * 5
        xMin = plotSettings[plotVar][1]
        xMax = plotSettings[plotVar][2]
        binwidth = (xMax - xMin)/nBins
        hSigArr = hSigArrDict[plotVar]  # one hist for each signal file
        hSig = TH1F(plotVar + "_sig_" + filename, plotVar + "_sig_" + \
                filename[21:31], nBins + 1, xMin, xMax)
        hSig.SetDirectory(0)
        hSigArr.append(hSig)

    hSigCutflow = hSigArrDict["cutflow"][fileNum]
    for i, cut in enumerate(cuts, start=1):
        hSigCutflow.GetXaxis().SetBinLabel(i, cut)

    hSigCutflow.Fill(cuts["no cut"])

    # ********** Additional cuts. ***********
    if findingSameFlavor:
        if muPreference:
            lepIndices = selectMuMu(event, l1MinOkPt=20, maxOkIso=0.3)
        else: lepIndices = selectElEl(event, l1MinOkPt=20, maxOkIso=0.3)
        if lepIndices is None: continue
    else:
        lepIndices = selectMuEl(event, maxOkIso=0.3)
        l1Flav = "muon"
        l2Flav = "electron"
        if lepIndices is None:
            lepIndices = selectElMu(event, maxOkIso=0.3)
            if lepIndices is None: continue
            l1Flav = "electron"
            l2Flav = "muon"
    l1Index = lepIndices[0]
    l2Index = lepIndices[1]
    sigCutflowHist[fileNum].Fill(cuts["dilepton"])

    # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
    # sigCutflowHist[fileNum].Fill(cuts["deltaR(ll)>0.3"])

    if event.nbtag > 1: continue
    sigCutflowHist[fileNum].Fill(cuts["nbtag<2"])

    if event.met_pt < 80: continue
    sigCutflowHist[fileNum].Fill(cuts["MET>80"])

    # veto (3rd lepton) checks:
    if findingSameFlavor:
        # event should not give valid muel or elmu pair
        if selectMuEl(event) is not None: continue
        if selectElMu(event) is not None: continue
    else:
        # event should not give valid mumu or elel pair
        if selectMuMu(event) is not None: continue
        if selectElEl(event) is not None: continue
    sigCutflowHist[fileNum].Fill(cuts["no 3rd lepton"])
    
    if event.njets >= 4: continue
    sigCutflowHist[fileNum].Fill(cuts["njets<4"])

    jMaxPt = 0
    for j in range(event.njets):
        if np.reshape(event.jet_pt,20)[j] > np.reshape(event.jet_pt,20)[jMaxPt]:
            jMaxPt = j
    dR_lep1_jet = deltaR(event, l1Flav, l1Index, "jet", jMaxPt)
    dR_lep2_jet = deltaR(event, l2Flav, l2Index, "jet", jMaxPt)

    # ********** Plotting. ***********
    for plotVar in plotVarArr:
        if plotVar == "cutflow": break
        hSig = hBkgdDict[plotVar]
        if plotVar[:4] == "lep1": 
            val = np.reshape(getattr(event, l1Flav+plotVar[4:]),20)[l1Index]
        elif plotVar[:4] == "lep2": 
            val = np.reshape(getattr(event, l2Flav+plotVar[4:]),20)[l2Index]
        elif plotVar[:3] == "jet":
            val = np.reshape(getattr(event, "jet"+plotVar[4:]),20)[jMaxPt]
        elif plotVar == "dR_lep1_jet": val = dR_lep1_jet
        elif plotVar == "dR_lep2_jet": val = dR_lep2_jet
        else: val = getattr(event, plotVar)
        if val <= xMax:
            hSig.Fill(val, 1)
        else: # overflow
            hSig.Fill(xMax - binwidth/2, 1)

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

        hSig.Sumw2()
        hSig.Scale(xsec * lumi / hSig.GetSumOfWeights())
        hSig.SetMinimum(1)
        hSig.SetMaximum(10**12)
        legend = legendDict[plotVar]
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad

    print "Num surviving events after each cut from sig %s:" % filename 
    for i in range(0,nCuts):
        # print hBkgd.GetXaxis().GetBinLabel(i+1),"S/sqrt(B):",\
        #         hSigArr[fileNum].GetBinContent(i+1)/sqrt(hBkgd.GetBinContent(i+1))
        print hBkgd.GetXaxis().GetBinLabel(i+1),hSigCutflow.GetBinContent(i+1)
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
                "plots/v3CutSequence/"+plotVar+"_"+l1Flav[:2]+l2Flav[:2]+".png"
        print "Saving image", imgName
        img = TImage.Create()
        img.FromPad(canvasDict[plotVar])
        img.WriteImage(imgName)
    print "Done."

