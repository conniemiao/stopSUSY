#!/usr/bin/env python

# Draws 1D hist for data for some variable, for the summed bkgd data and for each 
# of the signal files.
# Uses the root file outputted by makeSusyBkgd+SigRoot.py
# Uses xsec info from sig_SingleStop_files

from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

plotVar = "met_pt" # **** change this line for different vars

# copy in the output name from running makeSusyBkgd+SigRoot.py:
allDataFile = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/stopCut_02Bkgd_TTDiLept_02Sig_muel_baseline.root"
print "Plotting from "+allDataFile

plotSettings = { # [nBins,xMin,xMax,listForm]
        "lep1_pt":[100,0,400,False], 
        "lep1_eta":[100,-3,3,False],
        "lep1_phi":[100,-4,4,False],
        "lep1_relIso":[100,0,0.1,False],
        "lep2_pt":[100,0,400,False],
        "lep2_eta":[100,-4,4,False],
        "lep2_phi":[100,-4,4,False],
        "lep2_relIso":[100,0,0.1,False],
        "njets":[15,0,15,False],
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
        "genweight":[100,2.980,2.995,False],
        }
numSigFiles = int(allDataFile[64:66])
testMode = True 
if numSigFiles > 10: testMode = False 
nBins = plotSettings[plotVar][0]
if not testMode and nBins > 20: nBins = nBins * 5
xMin = plotSettings[plotVar][1]
xMax = plotSettings[plotVar][2]
listForm = plotSettings[plotVar][3] # only for some of the jet variables 

binwidth = (xMax - xMin)/nBins # include overflow bin
hBkgd = TH1D(plotVar + "_bkgd", plotVar + "_bkgd", nBins + 1, \
        xMin, xMax + binwidth)
lumi = 3000000 # luminosity = 3000 /pb = 3,000,000 /fb

c1 = TCanvas("c1","Plot",10,20,1000,700)
gStyle.SetOptStat(0) # don't show any stats

inFile = TFile.Open(allDataFile)

#--------------------------------------------------------------------------------#
# *************** Filling bkgd data summed together  ************
print "Plotting", plotVar, "from background."
xsec = 67.75 # production cross section
inTree = inFile.Get("tBkgd")

nentries = inTree.GetEntries()
print("nentries={0:d}".format(nentries))
assert nentries > 0, "You have no events in your tree..."

for count, event in enumerate(inTree):
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

title = plotVar + " ("+allDataFile[70:74]+", normalized to 3000 /pb)"
if allDataFile[-13:-5] == "baseline": title += ", baseline"
else: title += ", with cuts"
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

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVar + " from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
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

    inTree = inFile.Get("tSig"+str(fileNum))
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    assert nentries > 0, "You have no events in your tree..."

    hSigArr.append(TH1D(plotVar + "_sig_" + filename, plotVar + "_sig_" + \
            filename[19:24], nBins + 1, xMin, xMax + binwidth))
    
    for count, event in enumerate(inTree):
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

# if not testMode:
#     print "Saving image."
#     gSystem.ProcessEvents()
#     img = TImage.Create()
#     img.FromPad(c1)
#     img.WriteImage(plotVar + "_Bkgd+Sig.png")
print "Done. Press enter to finish."
raw_input()

