# NOTE: NEEDS 4 CMD LINE ARGS with values:
# testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut
#
# Retreives the plots drawn by plot1D.py and either displays them all or saves them
# all to .png files.
#
# Uses the root files outputted by hadding the output from makeNtuple.py
# Possible lastcuts: see below.

print "Importing modules."
import sys
from ROOT import gSystem, TFile, TH1D, TCanvas, TImage
from collections import OrderedDict
print "Beginning execution of", sys.argv

assert len(sys.argv) == 5, "need 4 command line args: testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut"

if sys.argv[1] == "test": testMode = True
elif sys.argv[1] == "all": testMode = False
else: assert False, "invalid test mode, need {test, all}"

if sys.argv[2] == "show": displayMode = True
elif sys.argv[2] == "save": displayMode = False
else: assert False, "invalid display mode, need {show, save}"
if not displayMode:
    gROOT.SetBatch(kTRUE) # prevent displaying canvases

# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
channel = sys.argv[3]
if channel == "mumu": channel = "MuMu"
elif channel == "muel": channel = "MuEl"
elif channel == "elel": channel = "ElEl"
else: assert False, "invalid channel, need {mumu, elel, muel}"

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
        ("MET>80",4),("nJet<4",5)])
lastcut = sys.argv[4]
assert lastcut in cuts, "invalid last cut %s" % lastcut

plotSettings = { # [nBins,xMin,xMax,units]
        "lep1_pt":[100,0,400,"[Gev]"],
        "lep1_eta":[100,-4,4,""],
        "lep1_phi":[100,-4,4,""],
        "lep1_relIso":[100,0,0.2,""],
        "lep1_mt":[100,0,500,"[GeV]"],
        "lep2_pt":[100,0,400,"[GeV]"],
        "lep2_eta":[100,-4,4,""],
        "lep2_phi":[100,-4,4,""],
        "lep2_relIso":[100,0,0.2,""],
        "lep2_mt":[100,0,500,"[GeV]"],
        "nJet":[10,0,10,""],
        "Jet_pt":[100,0,400,"[GeV]"], 
        "Jet_eta":[100,-3,3,""],
        "Jet_phi":[100,-4,4,""],
        "Jet_ht":[100,0,800,"[GeV]"],
        "nbtag":[5,0,5,""],
        "nbtagLoose":[10,0,10,""],
        "nbtagTight":[5,0,5,""],
        "dR_lep1_jet":[100,0,7,""],
        "dR_lep2_jet":[100,0,7,""],
        "MET_pt":[100,0,500,"[GeV]"],
        "mt_tot":[100,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        "mt_sum":[100,0,1000,"[GeV]"], # mt1 + mt2
        "m_eff":[100,0,1000,"[GeV]"], # ht + MET + pt1 + pt2
        "Jet_ht_div_sqrt_MET":[100,0,100,""],
        "mt_tot_div_sqrt_MET":[100,0,100,""],
        "m_eff_div_sqrt_MET":[100,0,100,""]
        }

# assemble hist file adr
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
        "plots/Run2/v1/plot1D"
histFileAdr = imgDir+"/plot1D_"
if testMode: histFileAdr += "test_"
else: histFileAdr += "all_"
histFileAdr += channel+"_"+lastcut
print histFileAdr

#--------------------------------------------------------------------------------#

# A = SS, nominal rel iso
# B = OS, nominal rel iso (signal)
# C = OS, inverted rel iso
# D = SS, inverted rel iso
try: histFileA = TFile.Open(histFileAdr+"_A.root")
except:
    sys.stderr.write("WARNING: couldn't open region A hists file")
    exit()
try: histFileC = TFile.Open(histFileAdr+"_C.root")
except:
    sys.stderr.write("WARNING: couldn't open region C hists file")
    exit()
try: histFileD = TFile.Open(histFileAdr+"_D.root")
except:
    sys.stderr.write("WARNING: couldn't open region D hists file")
    exit()

# the hMCPlotVarDicts map a plotVar to the total (non QCD) MC hist for that plotVar,
# for 1 region
hMCPlotVarDictA = {}
hMCPlotVarDictC = {}
hMCPlotVarDictD = {}
# the hDataPlotVarDicts map a plotVar to the data hist for that plotVar, for 1 region
hDataPlotVarDictA = {}
hDataPlotVarDictC = {}
hDataPlotVarDictD = {}
# hQCDPlotVarDict maps the plotVar to the QCD hist for that plotVar
hQCDPlotVarDict = {}
canvasDict = {}
for i, plotVar in enumerate(plotSettings):
    if testMode:
        if i >= 2: break
    if displayMode:
        c = TCanvas("c_"+plotVar,"Plot",10,20,1000,700)
        canvasDict.update({plotVar:c})
    nBins = plotSettings[plotVar][0]
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    hMCPlotVarDictA.update({plotVar:TH1D("MC_"+plotVar+"_A", "MC_"+plotVar+"_A", \
            nBins, xMin, xMax)})
    hMCPlotVarDictC.update({plotVar:TH1D("MC_"+plotVar+"_C", "MC_"+plotVar+"_C", \
            nBins, xMin, xMax)})
    hMCPlotVarDictD.update({plotVar:TH1D("MC_"+plotVar+"_D", "MC_"+plotVar+"_D", \
            nBins, xMin, xMax)})
    hDataPlotVarDictA.update({plotVar:TH1D("data_"+plotVar+"_A", \
            "data_"+plotVar+"_A", nBins, xMin, xMax)})
    hDataPlotVarDictC.update({plotVar:TH1D("data_"+plotVar+"_C", \
            "data_"+plotVar+"_C", nBins, xMin, xMax)})
    hDataPlotVarDictD.update({plotVar:TH1D("data_"+plotVar+"_D", \
            "data_"+plotVar+"_D", nBins, xMin, xMax)})

subprocesses = []
with open("bkgd_fileRedirector") as bkgdSubprocessesListFile:
    for subprocessLine in bkgdSubprocessesListFile:
        subprocessLine = subprocessLine.rstrip('\n').split(" ")
        subprocess = subprocessLine[0]
        if subprocess[0] == "#": continue

        if subprocess == "WJetsToLNu" or subprocess == "DYJetsToLL_M-50":
            for i in range(5):
                subprocesses.append(subprocess+"_"+str(i)+"Parton")
        else: subprocesses.append(subprocess)

#--------------------------------------------------------------------------------#
print
print "----------- Calculating QCD. -----------"

for i, plotVar in enumerate(plotSettings):
    if testMode:
        if i >= 2: break
    hMC_A = hMCPlotVarDictA[plotVar]
    hMC_C = hMCPlotVarDictC[plotVar]
    hMC_D = hMCPlotVarDictD[plotVar]
    hData_A = hDataPlotVarDictA[plotVar]
    hData_C = hDataPlotVarDictC[plotVar]
    hData_D = hDataPlotVarDictD[plotVar]
    for subprocess in subprocesses:
        hBkgd_A = histFileA.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_A.GetSumOfWeights() <= 0: hBkgd_A.Scale(0)
        elif i == 0: print "Adding", subprocess, "region A"
        hMC_A.Add(hBkgd_A)

        hBkgd_C = histFileC.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_C.GetSumOfWeights() <= 0: hBkgd_C.Scale(0)
        elif i == 0: print "Adding", subprocess, "region C"
        hMC_C.Add(hBkgd_C)

        hBkgd_D = histFileD.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_D.GetSumOfWeights() <= 0: hBkgd_D.Scale(0)
        elif i == 0: print "Adding", subprocess, "region D"
        hMC_D.Add(hBkgd_D)

    print "mc a", hMC_A.GetSumOfWeights()
    print "mc c", hMC_C.GetSumOfWeights()
    print "mc d", hMC_D.GetSumOfWeights()

    hData_A.Add(histFileA.Get("data_"+plotVar))
    hData_C.Add(histFileC.Get("data_"+plotVar))
    hData_D.Add(histFileD.Get("data_"+plotVar))

    print "data a", hData_A.GetSumOfWeights()
    print "data c", hData_C.GetSumOfWeights()
    print "data d", hData_D.GetSumOfWeights()

    hDiffA = hData_A - hMC_A
    hDiffD = hData_D - hMC_D
    hDiffCDivDiffD = hData_C - hMC_C
    hDiffCDivDiffD.Divide(hDiffD)

    hQCD = hDiffA * hDiffCDivDiffD

    if displayMode:
        c = canvasDict[plotVar]
        c.cd()
        hQCD.Draw()

    hQCDPlotVarDict.update({plotVar:hQCD})

histFileA.Close()
histFileC.Close()
histFileD.Close()

#--------------------------------------------------------------------------------#
if displayMode:
    for i, plotVar in enumerate(plotSettings):
        if testMode:
            if i >= 2: break
        c = canvasDict[plotVar]
        c.cd()
        c.SetLogy()
        c.Update()
    print "Done. Press enter to finish (plots not saved)."
    raw_input()
else:
    gSystem.ProcessEvents()
    outHistFileAdr = imgDir+"/QCDFromData_"
    if testMode: outHistFileAdr += "test_"
    else: outHistFileAdr += "all_"
    outHistFileAdr += channel+"_"+lastcut+"_.root"
    outHistFile = TFile(outHistFileAdr, "recreate")
    for plotVar in plotSettings:
        hQCDPlotVarDict[plotVar].Write()
    outHistFile.Close()
    print "Saved hists in", outHistFileAdr
    print "Done."
