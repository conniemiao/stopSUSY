# NOTE: NEEDS 6 CMD LINE ARGS with values:
# testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut,
# qcdType {qcdMC, qcdData}, region (only needed for qcdType=qcdMC) {A, B, C, D, any}
#
# Retreives the plots drawn by plot1D.py and either displays them all or saves them
# all to .png files.
#
# Uses the root files outputted by hadding the output from makeNtuple.py
# Possible lastcuts: see below.

print "Importing modules."
import sys
from ROOT import TFile, TCanvas, TImage
from ROOT import gSystem, gStyle, gROOT, kTRUE
from collections import OrderedDict
print "Beginning execution of", sys.argv

# location where the root file with all the plots was saved
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/plot1D"

assert len(sys.argv) >= 6, "need at least 5 command line args: testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut, qcdType {qcdMC, qcdData}, region (only needed for qcdType=qcdMC) {A, B, C, D, any}"

if sys.argv[1] == "test": testMode = True
elif sys.argv[1] == "all": testMode = False
else: assert False, "invalid test mode, need {test, all}"

if sys.argv[2] == "show": displayMode = True
elif sys.argv[2] == "save": displayMode = False
else: assert False, "invalid display mode, need {show, save}"

# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
channel = sys.argv[3]
if channel == "mumu": channel = "MuMu"
elif channel == "muel": channel = "MuEl"
elif channel == "elel": channel = "ElEl"
else: assert False, "invalid channel, need {mumu, elel, muel}"

cuts = OrderedDict([("baseline",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
        ("MET>80",4),("nJet<4",5)])
lastcut = sys.argv[4]
assert lastcut in cuts, "invalid last cut %s" % lastcut

qcdType = sys.argv[5]
assert qcdType == "qcdMC" or qcdType == "qcdData", \
        "invalid qcd type, need {qcdMC, qcdData}"

if qcdType == "qcdMC":
    region = sys.argv[6]
    assert region == "any" or region == "A" or region == "B" or region == "C" \
            or region == "D" or region == "any", \
            "invalid region, need {any, A, B, C, D, any}"
 
# assemble hist file adr
if qcdType == "qcdMC": histFileAdr = imgDir+"/QCDMC_plot1D_"
else: histFileAdr = imgDir+"/QCDData_plot1D_" 
if testMode: histFileAdr += "test_"
else: histFileAdr += "all_"
histFileAdr += channel+"_"+lastcut
if qcdType == "qcdMC": histFileAdr += "_"+region
histFileAdr += ".root"
print "Retrieving plots from", histFileAdr
histFile = TFile.Open(histFileAdr, "READ")

plotVars = {"lep1_pt", "lep1_eta", "lep1_phi", "lep1_relIso", "lep1_mt", 
"lep2_pt", "lep2_eta", "lep2_phi", "lep2_relIso", "lep2_mt", "nJet", "Jet_pt", 
"Jet_eta", "Jet_phi", "Jet_ht", "nbtag", "nbtagLoose", "nbtagTight", "dR_lep1_jet",
"dR_lep2_jet", "MET_pt", "mt_tot", "mt_sum", "m_eff", "Jet_ht_div_sqrt_MET",
"mt_tot_div_sqrt_MET", "m_eff_div_sqrt_MET", "fakeSort"}

# do the thing
if displayMode:
    for plotVar in plotVars:
        histFile.Get("c_"+plotVar).Draw()
    print "Done. Press enter to finish (plots not saved)."
    raw_input()
else:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases
    for plotVar in plotVars:
        gSystem.ProcessEvents()
        if not histFile.GetListOfKeys().Contains("c_"+plotVar): continue
        c = histFile.Get("c_"+plotVar)
        c.Draw()
        imgName = imgDir+"/"+plotVar+"_"+channel+"_"+lastcut
        if qcdType == "qcdMC": imgName += "_qcdMC_"+region+".png"
        else: imgName += "_qcdData.png"
        img = TImage.Create()
        img.FromPad(c)
        img.WriteImage(imgName)
        print "Saved image", imgName
    print "Done."
