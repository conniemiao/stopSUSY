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
from ROOT import TFile, TTree, TH1D, TCanvas, TPad, TLegend, TText, THStack, TLine
from ROOT import gSystem, gStyle, gROOT, kTRUE
from collections import OrderedDict
print "Beginning execution of", sys.argv

# location where the root file with all the qcdMC plots will be taken from, and also
# where the new qcdData control plots will be saved
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/plot1D"

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

cuts = OrderedDict([("baseline",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
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
        "nJet":[10,0.5,10.5,""],
        "Jet_pt":[100,0,400,"[GeV]"], 
        "Jet_eta":[100,-3,3,""],
        "Jet_phi":[100,-4,4,""],
        "Jet_ht":[100,0,800,"[GeV]"],
        "nbtag":[5,0.5,5.5,""],
        "nbtagLoose":[10,0.5,10.5,""],
        "nbtagTight":[5,0.5,5.5,""],
        "dR_lep1_jet":[100,0,7,""],
        "dR_lep2_jet":[100,0,7,""],
        "mt2":[100,0,150,"[GeV]"],
        "MET_pt":[100,0,500,"[GeV]"], 
        "mt_tot":[100,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        "mt_sum":[100,0,1000,"[GeV]"], # mt1 + mt2
        "m_eff":[100,0,1000,"[GeV]"], # ht + MET + pt1 + pt2
        "Jet_ht_div_sqrt_MET":[100,0,200,""],
        "mt_tot_div_sqrt_MET":[100,0,200,""],
        "m_eff_div_sqrt_MET":[100,0,200,""]
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

# the hBkgdMCPlotVarDicts map a plotVar to the total (non QCD) MC hist for that 
# plotVar, for 1 region
hBkgdMCPlotVarDictA = {}
hBkgdMCPlotVarDictC = {}
hBkgdMCPlotVarDictD = {}
# the hDataPlotVarDicts map a plotVar to the data hist for that plotVar, for 1 region
hDataPlotVarDictA = {}
hDataPlotVarDictC = {}
hDataPlotVarDictD = {}
# hQCDPlotVarDict maps the plotVar to the QCD hist for that plotVar
hQCDPlotVarDict = {}

canvasDict = {}
hRatioDict = {} # maps each plotVar to the ratio histogram
plotPadDict = {}
ratioPadDict = {}
ratioLineDict = {}
legendDict = {}
# hBkgdStacksDict maps plotVar to the stack of background
hBkgdStacksDict = {}
for plotVarNum, plotVar in enumerate(plotSettings):
    if testMode:
        if plotVarNum >= 2: break
    canvasDict.update({plotVar:TCanvas("c_"+plotVar,"c_"+plotVar,10,20,1000,700)})
    legendDict.update({plotVar:TLegend(.5,.75,.95,.90)})
    title = plotVar+" ("+channel+", cuts to "+lastcut+", region B)"
    hBkgdStacksDict.update({plotVar:THStack(plotVar+"_bkgdStack", title)})

    nBins = plotSettings[plotVar][0]
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    
    # need to initialize these as new histograms because want to add everything
    # into 1 hist
    hBkgdMCPlotVarDictA.update({plotVar:TH1D("MC_"+plotVar+"_A", "MC_A", \
            nBins, xMin, xMax)})
    hBkgdMCPlotVarDictC.update({plotVar:TH1D("MC_"+plotVar+"_C", "MC_C", \
            nBins, xMin, xMax)})
    hBkgdMCPlotVarDictD.update({plotVar:TH1D("MC_"+plotVar+"_D", "MC_D", \
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
            useNPartons = True # true if npartons were all handled in qcd mc
            for i in range(5):
                name = subprocess+"_"+str(i)+"Parton"
                # make sure all the npartons were actually run in plot1D_qcdMC
                histname = "bkgd_"+name+"_mt2"
                useNPartons = histFileB.GetListOfKeys().Contains(histname) and \
                        useNPartons
            if useNPartons:
                for i in range(5):
                    name = subprocess+"_"+str(i)+"Parton"
                    bkgdSubprocesses.append(name)
                    hBkgdSubprocessesPlotVarDict.update({name:{}})
            else: # W/DY Jets incl only
                # make sure this process was actually run in plot1D_qcdMC
                histname = "bkgd_"+subprocess+"_mt2"
                if not histFileB.GetListOfKeys().Contains(histname): continue
                bkgdSubprocesses.append(subprocess)
                hBkgdSubprocessesPlotVarDict.update({subprocess:{}})
        else: 
            # make sure this process was actually run in plot1D_qcdMC
            histname = "bkgd_"+subprocess+"_mt2"
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
        histname = "sig_"+subprocess[10:27]+"_mt2"
        if not histFileB.GetListOfKeys().Contains(histname): continue

        sigSubprocesses.append(subprocess)
        hSigSubprocessesPlotVarDict.update({subprocess:{}})

# hDataPlotVarDictB maps each plotVar to its region B (signal) data hist
hDataPlotVarDictB = {}
#--------------------------------------------------------------------------------#
print
print "----------- Calculating QCD. -----------"

prevBkgdColor = -1
for plotVarNum, plotVar in enumerate(plotSettings):
    if testMode:
        if plotVarNum >= 2: break
    hBkgdMC_A = hBkgdMCPlotVarDictA[plotVar]
    hBkgdMC_C = hBkgdMCPlotVarDictC[plotVar]
    hBkgdMC_D = hBkgdMCPlotVarDictD[plotVar]
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
        # hBkgd_A.Scale(0.5) # lol this is definitely not ok (!!!)
        hBkgdMC_A.Add(hBkgd_A)

        hBkgd_C = histFileC.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_C.GetSumOfWeights() <= 0: hBkgd_C.Scale(0)
        elif plotVarNum == 0: print "Adding", subprocess, "region C"
        # hBkgd_C.Scale(0.5) # lol this is definitely not ok (!!!)
        hBkgdMC_C.Add(hBkgd_C)

        hBkgd_D = histFileD.Get("bkgd_"+subprocess+"_"+plotVar)
        if hBkgd_D.GetSumOfWeights() <= 0: hBkgd_D.Scale(0)
        elif plotVarNum == 0: print "Adding", subprocess, "region D"
        # hBkgd_D.Scale(0.5) # lol this is definitely not ok (!!!)
        hBkgdMC_D.Add(hBkgd_D)

        # for final plotting (using plots from MC_B):
        histname = "bkgd_"+subprocess+"_"+plotVar
        hBkgd_B = histFileB.Get(histname)
        hBkgd_B.SetDirectory(0)
        # hBkgd_B.Scale(0.5) # lol this is definitely not ok (!!!)
        hBkgdSubprocessesPlotVarDict[subprocess].update({plotVar:hBkgd_B})
        if hBkgd_B.GetSumOfWeights() <= 0:
            if plotVarNum == 0: # just print once
                sys.stderr.write("WARNING: Region B sum of weights "+subprocess+\
                        " <= 0! ("+str(hBkgd_B.GetSumOfWeights())+")\n")
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
    hQCD.Add(hData_A, hBkgdMC_A, 1, -1) # hQCD = hDataA - hBkgdMC_A
    if hQCD.GetSumOfWeights() < 0:
        if plotVarNum == 0: # just print once
            sys.stderr.write("WARNING: Sum of weights hData_A - hBkgdMC_A < 0! ("+\
                    str(hQCD.GetSumOfWeights())+")\n")

    hDiffD = hData_D.Clone()
    hDiffD.Add(hBkgdMC_D, -1) # hDiffD = hData_D - hBkgdMC_D
    if hDiffD.GetSumOfWeights() < 0:
        if plotVarNum == 0: # just print once
            sys.stderr.write("WARNING: Sum of weights hData_D - hBkgdMC_D < 0! ("+\
                    str(hDiffD.GetSumOfWeights())+")\n")
    hQCD.Divide(hDiffD) # hQCD = hQCD/hDiffD

    hDiffC = hData_C.Clone()
    hDiffC.Add(hBkgdMC_C, -1) # hDiffC = hData_C - hBkgdMC_C
    if hDiffC.GetSumOfWeights() < 0:
        if plotVarNum == 0: # don't keep printing this
            sys.stderr.write("WARNING: Sum of weights hData_C - hBkgdMC_C < 0! ("+\
                    str(hDiffC.GetSumOfWeights())+")\n")
    hQCD.Multiply(hDiffC) # hQCD = hQCD * hDiffC

    hQCD.SetFillColor(colorQCD)
    hQCD.SetLineColor(colorQCD)
    hQCD.SetMarkerColor(colorQCD)
    hQCD.SetDirectory(0)
    legend = legendDict[plotVar]
    legend.AddEntry(hQCD, hQCD.GetTitle())
    hQCDPlotVarDict.update({plotVar:hQCD})
    hBkgdStack.Add(hQCD)

    # ********** Drawing. ***********
    plotPad = TPad("p_"+plotVar,"p_"+plotVar, 0.0, 0.3, 1.0, 1.0)
    plotPadDict[plotVar] = plotPad
    plotPad.SetNumber(0)
    plotPad.SetTicks(0, 0)
    plotPad.SetBottomMargin(0)
    plotPad.SetLeftMargin(0.1)
    plotPad.SetRightMargin(0.05)
    plotPad.SetFillColor(4000) # transparent
    c.cd()
    ratioPad = TPad("pRatio_"+plotVar,"pRatio_"+plotVar, 0.0, 0.0, 1.0, 0.3)
    ratioPadDict[plotVar] = ratioPad
    ratioPad.SetNumber(1)
    ratioPad.SetTopMargin(0.01)
    ratioPad.SetBottomMargin(0.25)
    ratioPad.SetLeftMargin(0.1)
    ratioPad.SetRightMargin(0.05)
    ratioPad.SetFillColor(4000) # transparent
    ratioPad.Draw()
    plotPad.Draw()

    plotPad.cd()
    plotPad.SetLogy()

    # ********** Background. ***********
    hBkgdStack.Draw("hist")
    unitsLabel = plotSettings[plotVar][3]
    hBkgdStack.GetXaxis().SetTitle(plotVar+" "+unitsLabel)
    hBkgdStack.GetYaxis().SetTitle("Number of Events, norm to 35921 /pb")
    hBkgdStack.SetMinimum(1)
    hBkgdStack.SetMaximum(10**8)
    nBins = plotSettings[plotVar][0]
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    # for ratio canvas:
    hMC = TH1D("allMC_"+plotVar, "allMC_"+plotVar, nBins, xMin, xMax)
    for hBkgdPlotVarDict in hBkgdSubprocessesPlotVarDict.values():
        hMC.Add(hBkgdPlotVarDict[plotVar])
    hMC.Add(hQCD)

    # ********** Signal. ***********
    for subprocess in sigSubprocesses:
        histname = "sig_"+subprocess[10:27]+"_"+plotVar
        hSig = histFileB.Get(histname)
        hSig.SetDirectory(0)
        hSigSubprocessesPlotVarDict[subprocess].update({plotVar:hSig})
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad

    # ********** Data. ***********
    histname = "data_"+plotVar
    if not histFileB.GetListOfKeys().Contains(histname):
        assert False, "Something's wrong, you have no data hist in your hists file!"
    hData = histFileB.Get(histname)
    if hData.GetLineColor() == 0: continue # something went wrong in plot1D_qcdMC
    hData.SetDirectory(0)
    hDataPlotVarDictB.update({plotVar:hData})
    legend.AddEntry(hData, hData.GetTitle())
    hData.Draw("* hist same") # same pad

    legend.SetTextSize(0.025)
    legend.SetNColumns(3)
    legend.Draw("same")

    # ********** Ratio canvas. ***********
    ratioPad.cd()
    ratioPad.SetGridy(1)
    hRatio = hData.Clone()
    hRatio.SetDirectory(0)
    hRatioDict[plotVar] = hRatio
    hRatio.Divide(hRatio, hMC)
    hRatio.SetMarkerStyle(20)
    hRatio.SetTitle("")
    hRatio.SetLabelSize(0.08,"Y")
    hRatio.SetLabelSize(0.08,"X")
    hRatio.GetYaxis().SetRangeUser(0,3.5)
    hRatio.GetYaxis().SetNdivisions(204)
    hRatio.GetYaxis().SetTitle("Obs/Exp    ")
    hRatio.SetTitle(";"+plotVar+" "+unitsLabel)
    hRatio.SetTitleSize(0.08,"Y")
    hRatio.SetTitleOffset(0.45,"Y")
    hRatio.SetTitleSize(0.08,"X")
    hRatio.SetTitleOffset(0.8,"X")
    line = TLine(xMin, 1.0, xMax, 1.0)
    ratioLineDict[plotVar] = line
    line.SetLineWidth(2)
    line.SetLineColor(2) # red
    hRatio.Draw("P")
    line.Draw()

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
        hDataPlotVarDictB[plotVar].Write()
    outHistFile.Close()
    print "Saved hists in", outHistFileAdr
    print "Done."
