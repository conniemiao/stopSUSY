#!/usr/bin/env python

# Draws hist for data corresponding to some variable for the summed bkgd data and
# for each of the signal files.
# Uses the root file outputted by makeSusyBkgd+SigRoot.py
# Uses xsec info from sig_SingleStop_files

from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

plotVar = "mtmu" # **** change this line for different vars
allDataFile = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/stopCut_02Bkgd_TTDiLept_02Sig.root"

plotSettings = { #[nBins,xMin,xMax]]
        "muon_px":[100,-300,300],
        "muon_py":[100,-300,300],
        "muon_pz":[100,-700,700],
        "muon_pt":[100,0,400], 
        "muon_eta":[100,-3,3],
        "muon_phi":[100,-4,4],
        "pfjet_px":[100,-300,300],
        "pfjet_py":[100,-300,300],
        "pfjet_pz":[100,-700,700],
        "pfjet_pt":[100,0,400], 
        "pfjet_eta":[100,-3,3],
        "pfjet_phi":[100,-4,4],
        "pfjet_btag":[10,0,10],
        "electron_px":[100,-300,300],
        "electron_py":[100,-300,300],
        "electron_pz":[100,-700,700],
        "electron_pt":[100,0,400],
        "electron_eta":[100,-4,4],
        "electron_phi":[100,-4,4],
        "mtmu":[100,0,500],
        "mtel":[100,0,500],
        "pfmet_pt":[100,0,500],
        "pfmet_ex":[100,-350,500],
        "pfmet_ey":[100,-250,350],
        "pfmet_ez":[100,-250,350],
        "genweight":[100,2.980,2.995],
        }
numSigFiles = int(allDataFile[64:66])
testMode = True 
if numSigFiles > 10: testMode = False 
nBins = plotSettings[plotVar][0]
if not testMode and nBins > 20: nBins = nBins * 5
xMin = plotSettings[plotVar][1]
xMax = plotSettings[plotVar][2]

binwidth = (xMax - xMin)/nBins # include overflow bin
histBkgd = TH1D(plotVar + "_bkgd", plotVar + "_bkgd", nBins + 1, \
        xMin, xMax + binwidth)
lumi = 150000 # luminosity = 150 /pb = 150,000 /fb

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

for count, entry in enumerate(inTree):
    if count % 500000 == 0: print("count = {0:d}".format(count))
    # *************** CHANGE BELOW FOR DIFFERENT PLOT VARS ************
    val = entry.mtmu # <----- CHANGE HERE
    if val <= xMax:
        histBkgd.Fill(val, 1)
    else: # overflow
        histBkgd.Fill(xMax + binwidth/2, 1)
histBkgd.Sumw2()

# rebinning
# print histBkgd.GetBinError((xMax+xMin)/2)/histBkgd.GetSumOfWeights()
# while histBkgd.GetBinError((xMax+xMin)/2) > 0.1*histBkgd.GetSumOfWeights():
#     histBkgd.Rebin(2)
#     print rebinned

histBkgd.SetTitle(plotVar)
histBkgd.GetXaxis().SetTitle(plotVar + " (normalized to 150000 /pb)")
histBkgd.GetYaxis().SetTitle("Number of Events")
histBkgd.Scale(xsec*lumi/histBkgd.GetSumOfWeights())
histBkgd.SetLineColor(1) # black
# histBkgd.SetStats(0)
histBkgd.Draw("hist")

c1.Update()

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVar + " from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

histSigArr = [] # one hist for each signal file 
for fileNum, line in enumerate(sigDataListFile):
    if fileNum + 1 > numSigFiles: break

    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename

    inTree = inFile.Get("tSig"+str(fileNum))
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))

    histSigArr.append(TH1D(plotVar + "_sig_" + filename, plotVar + "_sig_" + \
            filename[19:24], nBins + 1, xMin, xMax + binwidth))
    
    for count, entry in enumerate(inTree):
        if count % 500000 == 0: print("count={0:d}".format(count))
        # *************** CHANGE BELOW FOR DIFFERENT PLOT VARS ************
        val = entry.mtmu # <----- CHANGE HERE
        if val <= xMax:
            histSigArr[fileNum].Fill(val, 1)
        else: # overflow
            histSigArr[fileNum].Fill(xMax + binwidth/2, 1)

    histSigArr[fileNum].SetDirectory(0) # necessary to keep hist from closing

    histcolor = (fileNum + 1) % 3 + 2 # use colors 2 (red), 3 (green), 4 (blue)
    histSigArr[fileNum].Sumw2()
    histSigArr[fileNum].Scale(xsec * lumi /
            histSigArr[fileNum].GetSumOfWeights())
    histSigArr[fileNum].SetLineColor(histcolor) 
    histSigArr[fileNum].Draw("hist same") # same pad
    c1.Update()


#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************

legend = TLegend(.75,.80,.95,.95)
legend.AddEntry(histBkgd, histBkgd.GetTitle() + "_bkgd");
for fileNum, hist in enumerate(histSigArr):
    legend.AddEntry(histSigArr[fileNum], histSigArr[fileNum].GetTitle());
legend.Draw("same")

c1.SetLogy();
c1.Update();

if not testMode:
    print "Saving image."
    gSystem.ProcessEvents()
    img = TImage.Create()
    img.FromPad(c1)
    img.WriteImage(plotVar + "_Bkgd+Sig.png")
    print "Done."
else:
    print "Done. Press enter to finish."
    raw_input()

