# plots a variable from a root file outputted by makeSusySignalRoot.py

from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector
import numpy as np

# filename = "Stop_175_LSP1_small"
# inFile = TFile.Open("selectedMuEl_" + filename + ".root")
# plotVar = "nbtag" # **** change this line

filename = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/single-stop14TeV_R_220_A.root"
inFile = TFile.Open(filename)
plotVar = "muon_px" # **** change this line

plotSettings = { #[nBins,xMin,xMax]
        "muon_pt":[100,0,3],
        "muon_phi":[100,-4,4],
        "muon_eta":[100,-3,3],
        "muon_px":[100,-150,150],
        "muon_py":[100,-150,150],
        "muon_pz":[100,-450,350],
        "muon_charge":[100,-1.5,1.5],
        "muon_relIso":[100,0,0.35],
        "electron_pt":[100,0,130],
        "electron_phi":[100,-4,4],
        "electron_eta":[100,-4,4],
        "electron_px":[100,-120,120],
        "electron_pz":[100,-300,200],
        "electron_charge":[100,-1,1],
        "electron_relIso":[100,0,0.35],
        # "njets":[10,0,10],
        # "nbtag":[3,0,3],
        # "nbtagTight":[2,0,2],
        "pfmet_corr_x":[100,-350,500],
        "pfmet_corr_y":[100,-250,350],
        "primvertex_count":[35,0,35],
        "genweight":[100,2.980,2.995],
        "mtmu":[100,0,350],
        "mtel":[100,0,300]
        }

nBins = plotSettings[plotVar][0]
xMin = plotSettings[plotVar][1]
xMax = plotSettings[plotVar][2]
print "Plotting", plotVar, "from", filename

# inTree = inFile.Get("t")
inTree = inFile.Get("AC1B")
nentries = inTree.GetEntries()
print("nentries={0:d}".format(nentries))
nMax = nentries
# nMax = 5000
count = 0
verbose = 0

hist = TH1D(plotVar, plotVar, nBins, xMin, xMax)
for entry in inTree :
    if count > nMax : break
    if verbose > 1 or (count % 5000 == 0) : print("count={0:d}".format(count))
    count += 1

    # hist.Fill(entry.muon_px) # **** change this line

    val = list(entry.muon_px)
    if len(val)>0: hist.Fill(val[0],1) # **** change this line

c1 = TCanvas("c1","Plot",200,50,1000,700)
hist.GetXaxis().SetTitle(plotVar + " (GeV)")
hist.GetYaxis().SetTitle("Number of Events")
hist.Draw()
c1.Update()

print "Press enter to finish."
raw_input()

