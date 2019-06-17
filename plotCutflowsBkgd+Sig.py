#!/usr/bin/env python

# Overlays the cutflow hists for bkgd and each of the signal files.
# Uses the root file outputted by makeSusyBkgd+SigRoot.py
# Uses xsec info from sig_SingleStop_files

from ROOT import TFile, TTree, TH1F, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np

# copy in the output name from running makeSusyBkgd+SigRoot.py:
allDataFile = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/stopCut_02Bkgd_TTDiLept_02Sig_mumu.root"
print "Plotting from "+allDataFile

numSigFiles = int(allDataFile[64:66])

lumi = 3000000 # luminosity = 3000 /pb = 3,000,000 /fb
c1 = TCanvas("c1","Plot",10,20,1000,700)
gStyle.SetOptStat(0) # don't show any stats

inFile = TFile.Open(allDataFile, "READ")
xsec = 67.75 # production cross section

hBkgd = inFile.Get("bkgd_cutflow")
hBkgd.Sumw2()
hBkgd.SetTitle("Cutflow ("+allDataFile[70:74]+", norm to 3000 /pb)")
hBkgd.GetYaxis().SetTitle("Number of Events")
hBkgd.Scale(xsec*lumi/hBkgd.GetSumOfWeights())
hBkgd.SetLineColor(1) # black
hBkgd.SetMinimum(1)
hBkgd.SetMaximum(10**12)
hBkgd.Draw("hist")
c1.Update()


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

    hSigArr.append(inFile.Get("sig_"+filename[19:24]+"_cutflow"))
    hSigArr[fileNum].SetDirectory(0)

    hcolor = coloropts[fileNum % len(coloropts)]
    hSigArr[fileNum].SetLineColor(hcolor) 
    hmarkerstyle = markeropts[(fileNum/len(coloropts)) % len(markeropts)]
    hSigArr[fileNum].SetMarkerStyle(hmarkerstyle)
    hSigArr[fileNum].SetMarkerColor(hcolor)
    hlinestyle = linestyleopts[(fileNum/len(coloropts)/len(markeropts)) % \
            len(linestyleopts)]
    hSigArr[fileNum].SetLineStyle(hlinestyle)

    hSigArr[fileNum].Sumw2()
    hSigArr[fileNum].Scale(xsec * lumi /
            hSigArr[fileNum].GetSumOfWeights())
    hSigArr[fileNum].SetMinimum(1)
    hSigArr[fileNum].SetMaximum(10**12)
    hSigArr[fileNum].Draw("hist same") # same pad, draw marker
    c1.Update()

legend = TLegend(.70,.75,.90,.90)
legend.AddEntry(hBkgd, "bkgd_cutflow")
legend.SetTextSize(0.02)
for fileNum, h in enumerate(hSigArr):
    legend.AddEntry(hSigArr[fileNum], hSigArr[fileNum].GetTitle())
legend.Draw("same")

c1.SetLogy()
c1.Update()

print "Done. Press enter to finish."
raw_input()
