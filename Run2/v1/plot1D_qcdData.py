#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values:
# testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut
# Possible lastcuts: see below.
#
# Using the plots drawn by plot1D_qcdMC, estimates qcd from data for each of the
# plotVars. Then, draws each of the plotVars (signal region only) onto canvases using
# the data-estimated QCD and saves these and all of the individual hists into another
# output root.
#
# Uses the A, B, C, D root files outputted from plot1D_qcdMC.py
# Uses bkgd_fileRedirector
# Uses sig_fileRedirector

print "Importing modules."
import sys
from ROOT import TFile, TTree, TH1D, TCanvas, TImage, TLegend, TText, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE
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

# color for plotting : bkgd process name 
colorWJets = 38 # dark blue
colorDY = 46 # red
colorTTBar = 835 # teal 
colorSingleTop = 832  
colorTTX = 831 
colorDiboson = 806 #orange
colorQCD = 868 # light blue
processes = OrderedDict([(colorWJets,"W-Jets"), (colorDY,"Drell-Yan"), \
        (colorTTBar,"TTBar"), (colorSingleTop,"Single-Top"), (colorTTX, "TT+X"), \
        (colorDiboson,"Diboson"), (colorQCD, "QCD")])

# assemble hist file adr
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
        "plots/Run2/v1/plot1D"
histFileAdr = imgDir+"/QCDMC_plot1D_"
if testMode: histFileAdr += "test_"
else: histFileAdr += "all_"
histFileAdr += channel+"_"+lastcut
print histFileAdr

#--------------------------------------------------------------------------------#

# A = SS, nominal rel iso
# B = OS, nominal rel iso (signal)
# C = OS, inverted rel iso
# D = SS, inverted rel iso
try: histFileA = TFile.Open(histFileAdr+"_A.root", "READ")
except:
    sys.stderr.write("WARNING: couldn't open region A hists file")
    exit()
try: histFileB = TFile.Open(histFileAdr+"_B.root", "READ")
except:
    sys.stderr.write("WARNING: couldn't open region B (signal) hists file")
    exit()
try: histFileC = TFile.Open(histFileAdr+"_C.root", "READ")
except:
    sys.stderr.write("WARNING: couldn't open region C hists file")
    exit()
try: histFileD = TFile.Open(histFileAdr+"_D.root", "READ")
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
legendDict = {}
# hBkgdStacksDict maps plotVar to the stack of background
hBkgdStacksDict = {}
for plotVarNum, plotVar in enumerate(plotSettings):
    if testMode:
        if plotVarNum >= 2: break
    canvasDict.update({plotVar:TCanvas("c_"+plotVar,"c_"+plotVar,10,20,1000,700)})
    legendDict.update({plotVar:TLegend(.45,.75,.90,.90)})
    title = plotVar+" ("+channel+", cuts to "+lastcut+")"
    hBkgdStacksDict.update({plotVar:THStack(plotVar+"_bkgdStack", title)})

    nBins = plotSettings[plotVar][0]
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    hMCPlotVarDictA.update({plotVar:TH1D("MC_"+plotVar+"_A", "MC_A", \
            nBins, xMin, xMax)})
    hMCPlotVarDictC.update({plotVar:TH1D("MC_"+plotVar+"_C", "MC_C", \
            nBins, xMin, xMax)})
    hMCPlotVarDictD.update({plotVar:TH1D("MC_"+plotVar+"_D", "MC_D", \
            nBins, xMin, xMax)})
    hDataPlotVarDictA.update({plotVar:TH1D("data_"+plotVar+"_A", \
            "data_A", nBins, xMin, xMax)})
    hDataPlotVarDictC.update({plotVar:TH1D("data_"+plotVar+"_C", \
            "data_C", nBins, xMin, xMax)})
    hDataPlotVarDictD.update({plotVar:TH1D("data_"+plotVar+"_D", \
            "data_D", nBins, xMin, xMax)})
    hQCDPlotVarDict.update({plotVar:TH1D("bkgd_qcdData_"+plotVar, \
            "qcdData", nBins, xMin, xMax)})

# hBkgdSubprocessesPlotVarDict maps each bkgd subprocess to another dictionary,
# which maps each plotVar to a hist.
hBkgdSubprocessesPlotVarDict = {}
# list of "subprocesses" to loop on (ie separates nPartons for Wjets):
bkgdSubprocesses = [] 
with open("bkgd_fileRedirector") as bkgd_redirector:
    for subprocessLine in bkgd_redirector:
        subprocessLine = subprocessLine.rstrip('\n').split(" ")
        subprocess = subprocessLine[0]
        if subprocess[0] == "#": continue

        process = subprocessLine[1]
        if process == "QCD": continue

        if subprocess == "WJetsToLNu" or subprocess == "DYJetsToLL_M-50":
            for i in range(5):
                name = subprocess+"_"+str(i)+"Parton"
                
                # make sure this process was actually run in plot1D_qcdMC
                histname = "bkgd_"+name+"_m_eff"
                if not histFileB.GetListOfKeys().Contains(histname): continue

                bkgdSubprocesses.append(name)
                hBkgdSubprocessesPlotVarDict.update({name:{}})
        else: 
            # make sure this process was actually run in plot1D_qcdMC
            histname = "bkgd_"+histname+"_m_eff"
            if not histFileB.GetListOfKeys().Contains(histname): continue

            bkgdSubprocesses.append(subprocess)
            hBkgdSubprocessesPlotVarDict.update({subprocess:{}})

# hSigSubprocessesPlotVarDict maps the subprocess (signal type) to another dictionary
# which maps a plotVar to a hist
hSigSubprocessesPlotVarDict = {}
sigSubprocesses = []
with open("sig_fileRedirector") as sig_redirector:
    for subprocessLine in sig_redirector:
        subprocessLine = subprocessLine.rstrip('\n').split(" ")
        subprocess = subprocessLine[0]
        if subprocess[0] == "#": continue

        # make sure this process was actually run in plot1D_qcdMC
        histname = "sig_"+subprocess[10:27]+"_m_eff"
        if not histFileB.GetListOfKeys().Contains(histname): continue

        sigSubprocesses.append(subprocess)
        hSigSubprocessesPlotVarDict.update({subprocess:{}})

# hDataPlotVarDict maps each plotVar to a hist (which contains all the data from the
# process)
hDataPlotVarDict = {}
#--------------------------------------------------------------------------------#
print
print "----------- Calculating QCD. -----------"

prevBkgdColor = -1
for plotVarNum, plotVar in enumerate(plotSettings):
    if testMode:
        if plotVarNum >= 2: break
    hMC_A = hMCPlotVarDictA[plotVar]
    hMC_C = hMCPlotVarDictC[plotVar]
    hMC_D = hMCPlotVarDictD[plotVar]
    hData_A = hDataPlotVarDictA[plotVar]
    hData_C = hDataPlotVarDictC[plotVar]
    hData_D = hDataPlotVarDictD[plotVar]

    hBkgdStack = hBkgdStacksDict[plotVar]
    c = canvasDict[plotVar]
    c.cd()
    legend = legendDict[plotVar]

    for subprocess in bkgdSubprocesses:
        # for qcd estimation:
        hBkgd_A = histFileA.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_A.GetSumOfWeights() <= 0: hBkgd_A.Scale(0)
        elif plotVarNum == 0: print "Adding", subprocess, "region A"
        hMC_A.Add(hBkgd_A)

        hBkgd_C = histFileC.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_C.GetSumOfWeights() <= 0: hBkgd_C.Scale(0)
        elif plotVarNum == 0: print "Adding", subprocess, "region C"
        hMC_C.Add(hBkgd_C)

        hBkgd_D = histFileD.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_D.GetSumOfWeights() <= 0: hBkgd_D.Scale(0)
        elif plotVarNum == 0: print "Adding", subprocess, "region D"
        hMC_D.Add(hBkgd_D)

        # for final plotting:
        histname = "bkgd_"+subprocess+"_"+plotVar
        hBkgd_B = histFileB.Get(histname)
        hBkgd_B.SetDirectory(0)
        hBkgdSubprocessesPlotVarDict[subprocess].update({plotVar:hBkgd_B})
        hBkgdStack.Add(hBkgd_B) 
        hBkgdColor = hBkgd_B.GetFillColor() # colors/styles determined in plot1D_qcdMC
        if hBkgdColor == 0: continue # probably file not looped on in plot1D_qcdMC
        if hBkgdColor != prevBkgdColor:
            legend.AddEntry(hBkgd_B, processes[hBkgdColor])
        prevBkgdColor = hBkgdColor

    hData_A.Add(histFileA.Get("data_"+plotVar))
    hData_C.Add(histFileC.Get("data_"+plotVar))
    hData_D.Add(histFileD.Get("data_"+plotVar))

    # ********** hQCD calculation. ***********
    hQCD = hQCDPlotVarDict[plotVar]
    hQCD.Add(hData_C, hMC_C, 1, -1) # hQCD = hDataC - hMC_C
    if plotVarNum == 0 and hQCD.GetSumOfWeights() < 0:
        sys.stderr.write("WARNING: hData_C - hMC_C < 0!\n")
    hDiffD = hData_D.Clone()
    hDiffD.Add(hMC_D, -1) # hDiffD = hData_D - hMC_D
    if plotVarNum == 0 and hDiffD.GetSumOfWeights() < 0:
        sys.stderr.write("WARNING: hData_D - hMC_D < 0!\n")
    hQCD.Divide(hDiffD) # hQCD = hQCD/hDiffD
    hDiffA = hData_A.Clone()
    if plotVarNum == 0 and hDiffA.GetSumOfWeights() < 0:
        sys.stderr.write("WARNING: hData_A - hMC_A < 0!\n")
    hDiffA.Add(hMC_A, -1) # hDiffA = hData_A - hMC_A
    hQCD.Multiply(hDiffA) # hQCD = hQCD * hDiffA

    hQCD.SetFillColor(colorQCD)
    hQCD.SetLineColor(colorQCD)
    hQCD.SetMarkerColor(colorQCD)
    hQCD.SetDirectory(0)
    legend = legendDict[plotVar]
    legend.AddEntry(hQCD, hQCD.GetTitle())
    hQCDPlotVarDict.update({plotVar:hQCD})
    hBkgdStack.Add(hQCD)

    # ********** Drawing. ***********
    hBkgdStack.Draw("hist")
    unitsLabel = plotSettings[plotVar][3]
    hBkgdStack.GetXaxis().SetTitle(plotVar+" "+unitsLabel)
    hBkgdStack.GetYaxis().SetTitle("Number of Events, norm to 35921 /pb")
    hBkgdStack.SetMinimum(1)
    hBkgdStack.SetMaximum(10**12)

    for subprocess in sigSubprocesses:
        histname = "sig_"+subprocess[10:27]+"_"+plotVar
        hSig = histFileB.Get(histname)
        hSig.SetDirectory(0)
        hSigSubprocessesPlotVarDict[subprocess].update({plotVar:hSig})
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad

    histname = "data_"+plotVar
    if not histFileB.GetListOfKeys().Contains(histname):
        assert False, "Something's wrong, you have no data hist in your hists file!"
    hData = histFileB.Get(histname)
    if hData.GetLineColor() == 0: continue # something went wrong in plot1D_qcdMC
    hData.SetDirectory(0)
    hDataPlotVarDict.update({plotVar:hData})
    legend.AddEntry(hData, hData.GetTitle())
    hData.Draw("* hist same") # same pad

    legend.SetTextSize(0.017)
    legend.SetNColumns(3)
    legend.Draw("same")
    c.SetLogy()
    c.Update()

histFileA.Close()
histFileB.Close()
histFileC.Close()
histFileD.Close()

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print
print "Drawing."
if displayMode:
    print "Done. Press enter to finish (plots not saved)."
    raw_input()
else:
    gSystem.ProcessEvents()
    outHistFileAdr = imgDir+"/QCDData_plot1D_"
    if testMode: outHistFileAdr += "test_"
    else: outHistFileAdr += "all_"
    outHistFileAdr += channel+"_"+lastcut+".root"
    outHistFile = TFile(outHistFileAdr, "recreate")
    for plotVarNum, plotVar in enumerate(plotSettings):
        if testMode:
            if plotVarNum >= 2: break
        canvasDict[plotVar].Write()
        for subprocess in hBkgdSubprocessesPlotVarDict:
            hBkgdSubprocessesPlotVarDict[subprocess][plotVar].Write()
        hQCDPlotVarDict[plotVar].Write()
        for subprocess in hSigSubprocessesPlotVarDict:
            hSigSubprocessesPlotVarDict[subprocess][plotVar].Write()
        hDataPlotVarDict[plotVar].Write()
    outHistFile.Close()
    print "Saved hists in", outHistFileAdr
    print "Done."
