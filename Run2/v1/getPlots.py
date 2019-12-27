# NOTE: NEEDS 1 CMD LINE ARGS:
# file name of the root file containing the plots (just the file name, not the full
# path) and each plotVar to save
#
# Saves all the plots in the canvases in the file as .png images.

print "Importing modules."
import sys
from ROOT import TFile, TCanvas, TImage
from ROOT import gSystem, gStyle, gROOT, kTRUE
print "Beginning execution of", sys.argv

# location where the plots will be saved and from where plots file is retrieved
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/plot1D"

assert len(sys.argv) > 2, " need 1 arg with name of file containing plots and at least 1 plotVar"

histFileName = sys.argv[1]
histFileAdr = imgDir + "/" + histFileName
print "Retrieving plots from", histFileAdr
try:
    histFile = TFile.Open(histFileAdr, "READ")
except:
    print "Invalid hist file address."
    exit()

# plotVars = {"lep1_pt", "lep1_eta", "lep1_phi", "lep1_relIso", "lep1_mt", 
# "lep2_pt", "lep2_eta", "lep2_phi", "lep2_relIso", "lep2_mt", "nJet", "Jet_pt", 
# "Jet_eta", "Jet_phi", "Jet_ht", "nbtag", "nbtagLoose", "nbtagTight", "dR_lep1_jet",
# "dR_lep2_jet", "MET_pt", "mt2", "mll", "mt_tot", "mt_sum", "m_eff", "Jet_ht_div_sqrt_MET",
# "mt_tot_div_sqrt_MET", "m_eff_div_sqrt_MET", "fakeSort"}
plotVars = sys.argv[2:]

# do the thing
gROOT.SetBatch(kTRUE) # prevent displaying canvases
for plotVar in plotVars:
    gSystem.ProcessEvents()
    if not histFile.GetListOfKeys().Contains("c_"+plotVar): continue
    c = histFile.Get("c_"+plotVar)
    c.Draw()
    imgName = imgDir+"/"+histFileName[:-5]+"_"+plotVar+".png"
    img = TImage.Create()
    img.FromPad(c)
    img.WriteImage(imgName)
    print "Saved image", imgName
print "Done."
