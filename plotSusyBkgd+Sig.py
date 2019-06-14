#!/usr/bin/env python

# Draws hist for data corresponding to some variable for the summed bkgd data and
# for each of the signal files.
# Uses the root file outputted by makeSusyBkgd+SigRoot.py
# Uses xsec info from sig_SingleStop_files

from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

plotVar = "pfmet_pt" # **** change this line for different vars

# copy in the output name from running makeSusyBkgd+SigRoot.py:
allDataFile = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/stopCut_02Bkgd_TTDiLept_02Sig_muel_baseline.root"
print "Plotting from "+allDataFile

plotSettings = { #[nBins,xMin,xMax]]
        "lep1_px":[100,-300,300],
        "lep1_py":[100,-300,300],
        "lep1_pz":[100,-700,700],
        "lep1_pt":[100,0,400], 
        "lep1_eta":[100,-3,3],
        "lep1_phi":[100,-4,4],
        "pfjet_px":[100,-300,300],
        "pfjet_py":[100,-300,300],
        "pfjet_pz":[100,-700,700],
        "pfjet_pt":[100,0,400], 
        "pfjet_eta":[100,-3,3],
        "pfjet_phi":[100,-4,4],
        "pfjet_btag":[10,0,10],
        "lep2_px":[100,-300,300],
        "lep2_py":[100,-300,300],
        "lep2_pz":[100,-700,700],
        "lep2_pt":[100,0,400],
        "lep2_eta":[100,-4,4],
        "lep2_phi":[100,-4,4],
        "mtlep1":[100,0,500],
        "mtlep2":[100,0,500],
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

for count, event in enumerate(inTree):
    if count % 500000 == 0: print("count={0:d}".format(count))
    val = getattr(event, plotVar)
    if plotVar[:5] == "pfjet" and not plotVar[:11] == "pfjet_count":
        numjets = event.pfjet_count
        val = np.reshape(val, 20)
        for j in range(numjets):
            pfjetval = val[j]
            if pfjetval <= xMax:
                hBkgd.Fill(pfjetval, 1)
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

hBkgd.SetTitle(plotVar + " ("+allDataFile[70:74]+", normalized to 3000 /pb)")
hBkgd.GetXaxis().SetTitle(plotVar+" [GeV]")
hBkgd.GetYaxis().SetTitle("Number of Events")
hBkgd.Scale(xsec*lumi/hBkgd.GetSumOfWeights())
hBkgd.SetLineColor(1) # black
# hBkgd.SetStats(0)
hBkgd.Draw("h")

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

    hSigArr.append(TH1D(plotVar + "_sig_" + filename, plotVar + "_sig_" + \
            filename[19:24], nBins + 1, xMin, xMax + binwidth))
    
    for count, event in enumerate(inTree):
        if count % 500000 == 0: print("count={0:d}".format(count))
        val = getattr(event, plotVar)
        if plotVar[:5] == "pfjet" and not plotVar[:11] == "pfjet_count":
            numjets = event.pfjet_count
            val = np.reshape(val, 20)
            for j in range(numjets):
                pfjetval = val[j]
                if pfjetval <= xMax:
                    hSigArr[fileNum].Fill(pfjetval, 1)
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
    hSigArr[fileNum].Draw("hist same") # same pad, draw marker
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

