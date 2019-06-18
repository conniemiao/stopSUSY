#!/usr/bin/env python

# Draws 2D hist for data for some variable, for the summed bkgd data and for each 
# of the signal files (e.g. met vs. pt(l), met vs. mt(l)).
# Uses the root file outputted by makeSusyBkgd+SigRoot.py
# Uses xsec info from sig_SingleStop_files

from ROOT import TFile, TTree, TH2F, TCanvas, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

plotVarY = "met_pt" # y-axis variable
plotVarX = "lep1_pt" # x-axis variable

# copy in the output name from running makeSusyBkgd+SigRoot.py:
allDataFile = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/stopCut_02Bkgd_TTDiLept_02Sig_mumu.root"
print "Plotting from "+allDataFile

plotSettings = { # [nBins,xMin,xMax]
        "lep1_pt":[500,0,400], 
        "lep2_pt":[500,0,400],
        "mtlep1":[500,0,500],
        "mtlep2":[500,0,500],
        "met_pt":[500,0,500],
        }
numSigFiles = int(allDataFile[64:66])
testMode = True 
if numSigFiles > 10: testMode = False 
nBinsX = plotSettings[plotVarX][0]
if not testMode: nBinsX = nBinsX * 5
xMin = plotSettings[plotVarX][1]
xMax = plotSettings[plotVarX][2]
binwidthX = (xMax - xMin)/nBinsX # include overflow bin
nBinsY = plotSettings[plotVarY][0]
if not testMode: nBinsY = nBinsY * 5
yMin = plotSettings[plotVarY][1]
yMax = plotSettings[plotVarY][2]
binwidthY = (yMax - yMin)/nBinsY # include overflow bin

hBkgd = TH2F(plotVarY+"_"+plotVarX+"_bkgd", plotVarY+"_"+plotVarX+"_bkgd", \
        nBinsX + 1, xMin, xMax + binwidthX, nBinsY + 1, yMin, yMax + binwidthY)
lumi = 3000000 # luminosity = 3000 /pb = 3,000,000 /fb

c1 = TCanvas("c1","Plot",10,20,1000,700)
gStyle.SetOptStat(0) # don't show any stats

inFile = TFile.Open(allDataFile)

#--------------------------------------------------------------------------------#
# *************** Filling bkgd data summed together  ************
print "Plotting " + plotVarY + " vs. " + plotVarX + " from background."
xsec = 67.75 # production cross section
inTree = inFile.Get("tBkgd")

nentries = inTree.GetEntries()
print("nentries={0:d}".format(nentries))
assert nentries > 0, "You have no events in your tree..."

for count, event in enumerate(inTree):
    if count % 500000 == 0: print("count={0:d}".format(count))
    valX = getattr(event, plotVarX)
    valY = getattr(event,plotVarY)
    if valX > xMax:
        if valY > yMax: # x and y overflow
            hBkgd.Fill(xMax + binwidthX/2, yMax + binwidthY/2, 1)
        else: # x overflow, y in range
            hBkgd.Fill(xMax + binwidthX/2, valY, 1)
    else:
        if valY > yMax: # x in range, y overflow
            hBkgd.Fill(valX, yMax + binwidthY/2, 1)
        else: # x in range, y in range
            hBkgd.Fill(valX, valY, 1)
hBkgd.Sumw2()

# rebinning
# print hBkgd.GetBinError((xMax+xMin)/2)/hBkgd.GetSumOfWeights()
# while hBkgd.GetBinError((xMax+xMin)/2) > 0.1*hBkgd.GetSumOfWeights():
#     hBkgd.Rebin(2)
#     print rebinned

title = plotVarY + " v. "+plotVarX+" ("+allDataFile[70:74]+\
        ", normalized to 3000 /pb)"
if allDataFile[-13:-5] == "baseline": title += ", baseline"
else: title += ", with cuts"
hBkgd.SetTitle(title)
hBkgd.GetXaxis().SetTitle(plotVarX + " [GeV]")
hBkgd.GetYaxis().SetTitle(plotVarY + " [GeV]")
hBkgd.Scale(xsec*lumi/hBkgd.GetSumOfWeights())
hBkgd.SetLineColor(1) # black
# hBkgd.SetStats(0)
hBkgd.Draw("hist")
c1.Update()

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate hist  ************
print "Plotting " + plotVarY + " vs. " + plotVarX + " from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

coloropts = [2,4,3,6,7,9,28,46] # some good colors for lines

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

    hSigArr.append(TH2F(plotVarY+"_"+plotVarX+"_sig_"+filename[19:24], \
            plotVarY+"_"+plotVarX +"_sig_"+filename[19:24], nBinsX + 1, \
            xMin, xMax + binwidthX, nBinsY + 1, yMin, yMax + binwidthY))
    
    for count, event in enumerate(inTree):
        if count % 500000 == 0: print("count={0:d}".format(count))
        valX = getattr(event, plotVarX)
        valY = getattr(event,plotVarY)
        if valX > xMax:
            if valY > yMax: # x and y overflow
                hSigArr[fileNum].Fill(xMax + binwidthX/2, yMax + binwidthY/2, 1)
            else: # x overflow, y in range
                hSigArr[fileNum].Fill(xMax + binwidthX/2, valY, 1)
        else:
            if valY > yMax: # x in range, y overflow
                hSigArr[fileNum].Fill(valX, yMax + binwidthY/2, 1)
            else: # x in range, y in range
                hSigArr[fileNum].Fill(valX, valY, 1)

    hSigArr[fileNum].SetDirectory(0) # necessary to keep hist from closing

    hcolor = coloropts[fileNum % len(coloropts)]
    hSigArr[fileNum].SetLineColor(hcolor) 
    hSigArr[fileNum].SetMarkerColor(hcolor)

    hSigArr[fileNum].Sumw2()
    hSigArr[fileNum].Scale(xsec * lumi / hSigArr[fileNum].GetSumOfWeights())
    hSigArr[fileNum].Draw("hist same") # same pad, draw marker
    c1.Update()


#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************

legend = TLegend(.65,.75,.90,.90)
legend.AddEntry(hBkgd, plotVarY+"_"+plotVarX+"_bkgd")
legend.SetTextSize(0.02)
for fileNum, h in enumerate(hSigArr):
    legend.AddEntry(hSigArr[fileNum], hSigArr[fileNum].GetTitle())
legend.Draw("same")

c1.Update()

# if not testMode:
#     print "Saving image."
#     gSystem.ProcessEvents()
#     img = TImage.Create()
#     img.FromPad(c1)
#     img.WriteImage(plotVar + "_Bkgd+Sig.png")
print "Done. Press enter to finish."
raw_input()

