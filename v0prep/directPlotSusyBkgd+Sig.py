#!/usr/bin/env python

# (OUTDATED) draws hist for data corresponding to some variable for the 
# summed bkgd data and for each of the signal files
# uses bkgd data files listed in bkgd_TTDiLept_files
# uses signal data files listed in sig_SingleStop_files

from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

plotVar = "pfmet_pt" # **** change this line for different vars
testMode = True 

numBkgdFiles = 100 # max number of files to process
numSigFiles = 100
if testMode: 
    numBkgdFiles = 3 
    numSigFiles = 3 

plotSettings = { #[nBins,xMin,xMax]
        "muon_pt":[100,0,400, True], # <- list val
        "muon_px":[100,-300,300, True], # <- list val
        "muon_py":[100,-300,300, True], # <- list val
        "muon_pz":[100,-700,700, True], # <- list val
        "muon_phi":[100,-4,4, True],
        "muon_eta":[100,-3,3, True],
        "muon_charge":[100,-1.5,1.5, True],
        "muon_relIso":[100,0,0.35, True],
        "electron_pt":[100,0,400, True], # <- list val
        "electron_px":[100,-300,300, True], # <- list val
        "electron_py":[100,-300,300, True], # <- list val
        "electron_pz":[100,-700,700, True], # <- list val
        "electron_phi":[100,-4,4, True],
        "electron_eta":[100,-4,4, True],
        "electron_charge":[100,-1,1, True],
        "electron_relIso":[100,0,0.35, True],
        "pfmet_pt":[100,0,500, False], # <- float value
        "pfmet_ex":[100,-350,500, False],
        "pfmet_ey":[100,-250,350, False],
        "pfmet_ez":[100,-250,350, False],
        "genweight":[100,2.980,2.995, False],
        }

nBins = plotSettings[plotVar][0]
if (not testMode): nBins = nBins * 5
xMin = plotSettings[plotVar][1]
xMax = plotSettings[plotVar][2]
listFormat = plotSettings[plotVar][3]

binwidth = (xMax - xMin)/nBins # include overflow bin
histBkgd = TH1D(plotVar + "_bkgd", plotVar + "_bkgd", nBins + 1, \
        xMin, xMax + binwidth)
lumi = 150000 # luminosity = 150 /pb = 150,000 /fb

c1 = TCanvas("c1","Plot",10,20,1000,700)
gStyle.SetOptStat(0) # don't show any stats


# *************** Filling bkgd data summed together  ************
print "Plotting", plotVar, "from background."
bkgdDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/DESY_pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/"
bkgdDataListFile = open("bkgd_TTDiLept_files")

xsec = 67.75 # production cross section
filecount = 1
for line in bkgdDataListFile:
    if filecount > numBkgdFiles: break
    filecount += 1
    filename = line.rstrip()
    print filename
    inFile = TFile.Open(bkgdDataDir + filename)

    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    nMax = nentries
    if testMode: nMax = 5000
    count = 0
    
    for entry in inTree :
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))
        count += 1
        # *************** CHANGE BELOW FOR DIFFERENT PLOT VARS ************
        # ******* V1: STORED IN LIST FORMAT *****************
        if (listFormat):
            val = list(entry.muon_pt) # <----- CHANGE HERE
            if len(val)>0:
                if val[0] <= xMax:
                    histBkgd.Fill(val[0], 1)
                else: # overflow
                    histBkgd.Fill(xMax + binwidth/2, 1)
        # ******* V2: STORED IN FLOAT FORMAT *****************
        else:
            val = entry.pfmet_pt # <----- OR CHANGE HERE
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




# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVar + " from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

histSigArr = [] # one hist for each signal file 
filecount = 1
for line in sigDataListFile:
    if filecount > numSigFiles: break

    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename
    inFile = TFile.Open(sigDataDir + filename)

    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    nMax = nentries
    if testMode: nMax = 5000
    count = 0

    histSigArr.append(TH1D(plotVar + "_sig_" + filename, \
            plotVar + "_sig_" + filename[19:24], nBins + 1, xMin, xMax + binwidth))
    
    for entry in inTree :
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))
        count += 1
        # *************** CHANGE BELOW FOR DIFFERENT PLOT VARS ************
        # ******* V1: STORED IN LIST FORMAT *****************
        if (listFormat):
            val = list(entry.muon_pt) # <----- CHANGE HERE
            if len(val)>0:
                if val[0] <= xMax:
                    histSigArr[filecount-1].Fill(val[0], 1)
                else: # overflow
                    histSigArr[filecount-1].Fill(xMax + binwidth/2, 1)
        # ******* V2: STORED IN FLOAT FORMAT *****************
        else:
            val = entry.pfmet_pt # <----- OR CHANGE HERE
            if val <= xMax:
                histSigArr[filecount-1].Fill(val[0], 1)
            else: # overflow
                histSigArr[filecount-1].Fill(xMax + binwidth/2, 1)

    histSigArr[filecount-1].SetDirectory(0) # necessary to keep hist from closing

    histcolor = filecount % 3 + 2 # use colors 2 (red), 3 (green), 4 (blue)
    histSigArr[filecount-1].Sumw2()
    histSigArr[filecount-1].Scale(xsec * lumi /
            histSigArr[filecount-1].GetSumOfWeights())
    histSigArr[filecount-1].SetLineColor(histcolor) 
    histSigArr[filecount-1].Draw("hist same") # same pad
    c1.Update()

    filecount += 1



# *************** Wrap up. *******************

legend = TLegend(.75,.80,.95,.95)
legend.AddEntry(histBkgd, histBkgd.GetTitle() + "_bkgd");
for filecount, hist in enumerate(histSigArr):
    legend.AddEntry(histSigArr[filecount-1], histSigArr[filecount-1].GetTitle());
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

