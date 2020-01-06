#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values:
# testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}
#
# Using the plots drawn by plot1D_fakeRegions, rescales everything in the SR by the
# ratio of mc to data from the mt2 plot in CR1B.
# 
# Uses the cr1b and sr root files outputted from plot1D_fakeRegions.py
# Uses bkgd_fileRedirector
# Uses sig_fileRedirector

print "Importing modules."
import sys
assert len(sys.argv) == 4, "need 3 command line args: testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}"

from ROOT import TFile, TTree, TH1D, TCanvas, TPad, TLegend, TText, THStack, TLine
from ROOT import gSystem, gStyle, gROOT, kTRUE
from collections import OrderedDict
print "Beginning execution of", sys.argv

# location where the cutflow stats were saved
statsDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/cutflow_stats"
# location where the root file with all the fakeRegions plots will be taken from, and also
# where the new  plots will be saved
imgDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/plot1D"

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

cuts = OrderedDict([("baseline",0), ("allCuts",1)])
lastcut = "allCuts"
assert lastcut in cuts, "invalid last cut %s" % lastcut
nCuts = cuts[lastcut]+1

plotSettings = { # [nBins,xMin,xMax,units]
        #"lep1_pt":[100,0,400,"[Gev]"],
        #"lep1_eta":[100,-4,4,""],
        #"lep1_phi":[100,-4,4,""],
        #"lep1_relIso":[100,0,0.2,""],
        #"lep1_mt":[100,0,500,"[GeV]"],
        #"lep2_pt":[100,0,400,"[GeV]"],
        #"lep2_eta":[100,-4,4,""],
        #"lep2_phi":[100,-4,4,""],
        #"lep2_relIso":[100,0,0.2,""],
        #"lep2_mt":[100,0,500,"[GeV]"],
        #"nJet":[10,0.5,10.5,""],
        #"Jet_pt":[100,0,400,"[GeV]"], 
        #"Jet_eta":[100,-3,3,""],
        #"Jet_phi":[100,-4,4,""],
        #"Jet_ht":[100,0,800,"[GeV]"],
        #"nbtag":[5,0.5,5.5,""],
        #"nbtagLoose":[10,0.5,10.5,""],
        #"nbtagTight":[5,0.5,5.5,""],
        #"dR_lep1_jet":[100,0,7,""],
        #"dR_lep2_jet":[100,0,7,""],
        "mt2":[100,0,150,"[GeV]"],
        #"MET_pt":[100,0,500,"[GeV]"], 
        #"mt_tot":[100,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        #"mt_sum":[100,0,1000,"[GeV]"], # mt1 + mt2
        #"m_eff":[100,0,1000,"[GeV]"], # ht + MET + pt1 + pt2
        #"Jet_ht_div_sqrt_MET":[100,0,200,""],
        #"mt_tot_div_sqrt_MET":[100,0,200,""],
        #"m_eff_div_sqrt_MET":[100,0,200,""]
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
histFileAdr = imgDir+"/fakeRegions_plot1D_"
if testMode: histFileAdr += "test_"
else: histFileAdr += "all_"
histFileAdr += channel+"_"+lastcut
print histFileAdr

#--------------------------------------------------------------------------------#

try: histFileSR = TFile.Open(histFileAdr+"_sr.root", "READ")
except:
    sys.stderr.write("WARNING: couldn't open sr hists file")
    exit()
try: histFileCR1B = TFile.Open(histFileAdr+"_cr1b.root", "READ")
except:
    sys.stderr.write("WARNING: couldn't open cr1b hists file")
    exit()

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
    title = plotVar+" ("+channel+", cuts to "+lastcut+", sr with correction from CR1b)"
    hBkgdStacksDict.update({plotVar:THStack(plotVar+"_bkgdStack", title)})

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
        if subprocess == "WJetsToLNu" or subprocess == "DYJetsToLL_M-50":
            useNPartons = True # true if npartons were all handled in fakeRegions
            for i in range(5):
                name = subprocess+"_"+str(i)+"Parton"
                # make sure all the npartons were actually run in plot1D_fakeRegions
                histname = "bkgd_"+name+"_mll"
                useNPartons = histFileSR.GetListOfKeys().Contains(histname) and \
                        useNPartons
            if useNPartons:
                for i in range(5):
                    name = subprocess+"_"+str(i)+"Parton"
                    bkgdSubprocesses.append(name)
                    hBkgdSubprocessesPlotVarDict.update({name:{}})
            else: # W/DY Jets incl only
                # make sure this process was actually run in plot1D_fakeRegions
                histname = "bkgd_"+subprocess+"_mll"
                if not histFileSR.GetListOfKeys().Contains(histname): continue
                bkgdSubprocesses.append(subprocess)
                hBkgdSubprocessesPlotVarDict.update({subprocess:{}})
        else: 
            # make sure this process was actually run in plot1D_fakeRegions
            histname = "bkgd_"+subprocess+"_mll"
            if not histFileSR.GetListOfKeys().Contains(histname): continue
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

        # make sure this process was actually run in plot1D_fakeRegions
        histname = "sig_"+subprocess[10:27]+"_mll"
        if not histFileSR.GetListOfKeys().Contains(histname): continue

        sigSubprocesses.append(subprocess)
        hSigSubprocessesPlotVarDict.update({subprocess:{}})

# hDataPlotVarDict maps each plotVar to its region SR (signal) data hist
hDataPlotVarDict = {}
#--------------------------------------------------------------------------------#
print
print "----------- Rescaling by CR1b correction. -----------"

# Calculate the correction factor from mt2 in CR1B
totBkgdSumWeights = 0
for subprocess in bkgdSubprocesses:
    histname = "bkgd_"+subprocess+"_mt2"
    hBkgd_cr1b = histFileCR1B.Get(histname)
    hBkgd_cr1b.SetDirectory(0)
    totBkgdSumWeights += hBkgd_cr1b.GetSumOfWeights()
totDataSumWeights = histFileCR1B.Get("data_mt2").GetSumOfWeights()
cr1bScaleFactor = totDataSumWeights/totBkgdSumWeights
print "Scaling by", cr1bScaleFactor

prevBkgdColor = -1
for plotVarNum, plotVar in enumerate(plotSettings):
    if testMode:
        if plotVarNum >= 2: break

    if not histFileSR.GetListOfKeys().Contains("bkgd_"+bkgdSubprocesses[0]+"_"+\
            plotVar): continue
    hBkgdStack = hBkgdStacksDict[plotVar]
    c = canvasDict[plotVar]
    c.cd()
    legend = legendDict[plotVar]

    for subprocess in bkgdSubprocesses:
        histname = "bkgd_"+subprocess+"_"+plotVar
        if not histFileSR.GetListOfKeys().Contains(histname): continue
        hBkgd = histFileSR.Get(histname)
        hBkgd_new = hBkgd.Clone()
        hBkgd_new.Scale(cr1bScaleFactor)
        hBkgd = hBkgd_new
        hBkgd.SetDirectory(0)
        hBkgdSubprocessesPlotVarDict[subprocess].update({plotVar:hBkgd})
        if plotVarNum == 0: print "Adding", subprocess
        hBkgdStack.Add(hBkgd) 
        hBkgdColor = hBkgd.GetFillColor() # colors/styles determined in plot1D_fakeRegions
        if hBkgdColor == 0: continue # probably file not looped on in plot1D_fakeRegions
        if hBkgdColor != prevBkgdColor:
            legend.AddEntry(hBkgd, processes[hBkgdColor])
        prevBkgdColor = hBkgdColor

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
    hMC = TH1D("allMC_"+plotVar, "allMC_"+plotVar, nBins, xMin, xMax)
    for hBkgdPlotVarDict in hBkgdSubprocessesPlotVarDict.values():
        hMC.Add(hBkgdPlotVarDict[plotVar]) # for ratio canvas

    # ********** Signal. ***********
    for subprocess in sigSubprocesses:
        histname = "sig_"+subprocess[10:27]+"_"+plotVar
        hSig = histFileSR.Get(histname)
        hSig.SetDirectory(0)
        hSigSubprocessesPlotVarDict[subprocess].update({plotVar:hSig})
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad

    # ********** Data. ***********
    histname = "data_"+plotVar
    if not histFileSR.GetListOfKeys().Contains(histname):
        assert False, "Something's wrong, you have no data hist in your hists file!"
    hData = histFileSR.Get(histname)
    if hData.GetLineColor() == 0: continue # something went wrong in plot1D_fakeRegions
    hData.SetDirectory(0)
    hDataPlotVarDict.update({plotVar:hData})
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

histFileSR.Close()
histFileCR1B.Close()

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print
if displayMode:
    print "Done. Press enter to finish (plots not saved)."
    raw_input()
else:
    gSystem.ProcessEvents()
    outHistFileAdr = imgDir+"/fakeRegionsCR1bScaled_plot1D_"
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
        for subprocess in hSigSubprocessesPlotVarDict:
            hSigSubprocessesPlotVarDict[subprocess][plotVar].Write()
        hDataPlotVarDict[plotVar].Write()
    outHistFile.Close()
    print "Saved hists in", outHistFileAdr
    print "Done."
