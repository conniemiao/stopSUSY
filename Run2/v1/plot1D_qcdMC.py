#!/usr/bin/env python

# NOTE: NEEDS 5 CMD LINE ARGS with values:
# testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, 
# lastcut, region {any, A, B, C, D}
# A = SS, nominal rel iso
# B = OS, nominal rel iso (signal)
# C = OS, inverted rel iso
# D = SS, inverted rel iso
# Possible lastcuts: listed in cuts below.
#
# Implements additional cuts and then draws 1D hist for data for each variable in
# plotSettings for the summed bkgd data and for each of the signal files. Saves
# them in a histogram root file. Also saves both a .txt and a .hdf (pandas dataframe)
# version of the cutflow if lastcut is nJet<4. Uses MC for QCD.
#
# Uses the root files outputted by hadding the output from makeNtuple.py
# Uses bkgd_fileRedirector
# Uses sig_fileRedirector
# Uses data_fileRedirector
# Uses stopSelection.py

import sys, os
assert len(sys.argv) == 6, "need 5 command line args: testMode {test, all}, displayMode {show, save}, channel {mumu, elel, muel}, lastcut, region {A, B, C, D, any}"

print "Importing modules."
from ROOT import TFile, TTree, TH1D, TCanvas, TImage, TLegend, TText, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR
from stopSelection import isRegionA, isRegionB, isRegionC, isRegionD
from collections import OrderedDict
from math import sqrt, cos
import numpy as np
import pandas as pd
import time
print "Beginning execution of", sys.argv

# location of input ntuples
# myDataDir = "/eos/user/a/alkaloge/Connie/"
myDataDir = "/eos/user/c/cmiao/private/myDataSusy/Run2"
# location where cutflow stats will be saved
# statsDir = "/afs/cern.ch/user/a/alkaloge/work/Connie/CMSSW_10_2_9/src/stopSUSY/Run2/v1/cutflow_stats"
statsDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/cutflow_stats"
# location where the root file with all the control plots will be saved
# imgDir = "/afs/cern.ch/user/a/alkaloge/work/Connie/CMSSW_10_2_9/src/stopSUSY/Run2/v1/plot1D"
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
if sys.argv[3] == "mumu":
    findingSameFlavor = True
    muPref = True
    l1Flav = "Muon"
    l2Flav = "Muon"
    dataProcess = "DoubleMuon"
elif sys.argv[3] == "elel":
    findingSameFlavor = True
    muPref = False
    l1Flav = "Electron"
    l2Flav = "Electron"
    dataProcess = "DoubleEG"
elif sys.argv[3] == "muel":
    findingSameFlavor = False
    muPref = False
    l1Flav = "Muon"
    l2Flav = "Electron"
    dataProcess = "MuonEG"
else: assert False, "invalid channel, need {mumu, elel, muel}"
channel = l1Flav[:2] + l2Flav[:2]

# cuts = OrderedDict([("baseline",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
#         ("MET>80",4),("nJet<4",5)])
cuts = OrderedDict([("baseline",0), ("no3rdlept",1), ("nbtag<2",2), ("MET>80",3), \
        ("nJet<4",4)])
lastcut = sys.argv[4]
assert lastcut in cuts, "invalid last cut %s" % lastcut
nCuts = cuts[lastcut]+1

region = sys.argv[5]
assert region == "any" or region == "A" or region == "B" or region == "C" or region == "D", "invalid region, need {any, A, B, C, D}"

#--------------------------------------------------------------------------------#

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
        "MET_pt":[100,0,500,"[GeV]"], 
        "mt_tot":[100,0,1000,"[GeV]"], # sqrt(mt1^2 + mt2^2)
        "mt_sum":[100,0,1000,"[GeV]"], # mt1 + mt2
        "m_eff":[100,0,1000,"[GeV]"], # ht + MET + pt1 + pt2
        "Jet_ht_div_sqrt_MET":[100,0,200,""],
        "mt_tot_div_sqrt_MET":[100,0,200,""],
        "m_eff_div_sqrt_MET":[100,0,200,""]
        }

# produced particle -> labeled particle
# heavy quarks (c,b,t), light quarks (d,u,s), g: gluon
# tau->mu,e, e->e, and mu->mu are allowed, others are fake in some way.
# fakeTypes_idDict maps the genPartFlav tag number to the label
fakeTypes_idDict = OrderedDict([(1, "e->e/mu->mu"), (15, "tau->e/mu"), \
        (22, "photon->e"), (5, "b->e/mu"), (4, "c->e/mu"), \
         (3, "d/u/s/unknown->e/mu"), (0, "unmatched")])
# fakeTypes_orderDict maps id of the fake type to the index in the histogram
fakeTypes_orderDict = OrderedDict()
for i, fakeTypeId in enumerate(fakeTypes_idDict):
    fakeTypes_orderDict[fakeTypeId] = i

nFakeTypes = len(fakeTypes_idDict)

# bkgd process name : color for plotting
colorWJets = 38 # dark blue
colorDY = 46 # red
colorTTBar = 835 # teal 
colorSingleTop = 832  
colorTTX = 831 
colorDiboson = 806 # orange
colorQCD = 868 # light blue
processes = OrderedDict([("W-Jets",colorWJets), ("Drell-Yan",colorDY), \
        ("TTBar",colorTTBar), ("Single-Top",colorSingleTop), ("TT+X",colorTTX), \
        ("Diboson",colorDiboson), ("QCD", colorQCD)])

canvasDict = {}
legendDict = {}
hBkgdStacksDict = {} # maps plotVar to the stack of background
for plotVar in plotSettings: # add an entry to the plotVar:hist dictionary
    canvasDict.update({plotVar:TCanvas("c_"+plotVar,"c_"+plotVar,10,20,1000,700)})
    legendDict.update({plotVar:TLegend(.45,.75,.90,.90)})
    title = plotVar+" ("+channel+", cuts to "+lastcut+", region "+region+")"
    hBkgdStacksDict.update({plotVar:THStack(plotVar+"_bkgdStack", title)})
title = "cutflow ("+channel+", region "+region+")"
hBkgdCutflowStack = THStack("cutflow_bkgdStack", title)
c_fakeSort = TCanvas("c_fakeSort","c_fakeSort",10,20,1000,700)
legend_fakeSorting = TLegend(.45,.75,.90,.90)
title = "fakeSorting ("+channel+", cuts to "+lastcut+", region "+region+")"
hFakeSortingStack = THStack("fakeSorting", title)

# limit the number of files to process (other than what is commented out in the file
# redirector)
numBkgdFiles = float("inf")  # note: must loop over all files to have correct xsec
numSigFiles = 3

lumi = 35921 # 2016 lumi in /pb
WJets_kfactor = 1.221
DYJets_kfactor = 1.1637

gStyle.SetOptStat(0) # don't show any stats

#--------------------------------------------------------------------------------#
# *************** Filling bkgd mc summed together  ************
print
print "----------- Plotting from background. -----------"

# hBkgdCutflowDict maps every subprocess to an hBkgdCutflow that contains the cutflow
# from all the ntuples for that subprocess.
hBkgdCutflowDict = {}
# hBkgdCutsCountDict maps every process to an array of size nCuts that keeps track
# of the num evts remaining after each cut for that process
hBkgdCutsCountDict = {}
# hFakeSorting maps every subprocess to an hFakeSorting that contains the histogram
# that sorts the events into different fake processes
hFakeSortingDict = {}
# hBkgdSubprocessesPlotVarDict maps each bkgd subprocess to another dictionary,
# which maps each plotVar to a hist.
hBkgdSubprocessesPlotVarDict = {}
WNJetsXsecs = [47297.3] # first entry: W0Jets xsec (not in file redirector)
DYNJetsXsecs = [4263.5]  # first entry: DY0Jets xsec (not in file redirector)
with open("bkgd_fileRedirector") as bkgd_redirector:
    for subprocessLine in bkgd_redirector:
        subprocessLine = subprocessLine.rstrip('\n').split(" ")
        subprocess = subprocessLine[0]
        if subprocess[0] == "#": continue

        process = subprocessLine[1]

        if subprocess[0] == "W" and subprocess[2:] == "JetsToLNu":
            WNJetsXsecs.append(float(subprocessLine[2]))
        elif subprocess[:2] == "DY" and subprocess[3:] == "JetsToLL_M-50":
            DYNJetsXsecs.append(float(subprocessLine[2]))

        if subprocess != "WJetsToLNu" and subprocess != "DYJetsToLL_M-50":
            hBkgdSubprocessesPlotVarDict.update({subprocess:{}})
            hBkgdCutflow = TH1D("cutflow_"+subprocess+"_bkgd", subprocess, nCuts, 0,\
                    nCuts)
            hBkgdCutflowDict.update({subprocess:hBkgdCutflow})
            hFakeSorting = TH1D("fakeSorting_"+subprocess, subprocess, nFakeTypes, \
                    0, nFakeTypes)
            hFakeSortingDict.update({subprocess:hFakeSorting})

# deal with W/DYJets special cases:
subprocess = "WJetsToLNu"
WJetsIncl_only = False # default: running on both WJets incl and WnJets 
if len(WNJetsXsecs) != 5:
    WJetsIncl_only = True # only running on WJets inclusive
    hBkgdSubprocessesPlotVarDict.update({subprocess:{}})
    hBkgdCutflow = TH1D("cutflow_"+subprocess+"_bkgd", subprocess, nCuts, 0, nCuts)
    hBkgdCutflowDict.update({subprocess:hBkgdCutflow})
    hFakeSorting = TH1D("fakeSorting_"+subprocess, subprocess, nFakeTypes, 0, \
            nFakeTypes)
    hFakeSortingDict.update({subprocess:hFakeSorting})
else:
    for i in range(5):
        name = subprocess+"_"+str(i)+"Parton"
        hBkgdSubprocessesPlotVarDict.update({name:{}})
        hBkgdCutflow = TH1D("cutflow_"+name+"_bkgd", name, nCuts, 0, nCuts)
        hBkgdCutflowDict.update({name:hBkgdCutflow})
        hFakeSorting = TH1D("fakeSorting_"+name, name, nFakeTypes, 0, nFakeTypes)
        hFakeSortingDict.update({name:hFakeSorting})
subprocess = "DYJetsToLL_M-50"
DYJetsIncl_only = False # default: running on both DYJets incl and DYnJets 
if len(DYNJetsXsecs) != 5:
    DYJetsIncl_only = True # only running on DYJets inclusive
    hBkgdSubprocessesPlotVarDict.update({subprocess:{}})
    hBkgdCutflow = TH1D("cutflow_"+subprocess+"_bkgd", subprocess, nCuts, 0, nCuts)
    hBkgdCutflowDict.update({subprocess:hBkgdCutflow})
    hFakeSorting = TH1D("fakeSorting_"+subprocess, subprocess, nFakeTypes, 0, \
            nFakeTypes)
    hFakeSortingDict.update({subprocess:hFakeSorting})
else:
    for i in range(5):
        name = subprocess+"_"+str(i)+"Parton"
        hBkgdSubprocessesPlotVarDict.update({name:{}})
        hBkgdCutflow = TH1D("cutflow_"+name+"_bkgd", "cutflow_"+name, nCuts, 0, \
                nCuts)
        hBkgdCutflowDict.update({name:hBkgdCutflow})
        hFakeSorting = TH1D("fakeSorting_"+name, "fakeSorting_"+name, nFakeTypes, 0,\
                nFakeTypes)
        hFakeSortingDict.update({name:hFakeSorting})

# make/update the histgorams for all the subprocesses
for name in hBkgdSubprocessesPlotVarDict: # same as names in hBkgdCutflowDict
    hBkgdCutflow = hBkgdCutflowDict[name]
    hFakeSorting = hFakeSortingDict[name]
    for c, cut in enumerate(cuts, start=1):
        if c>nCuts: break
        hBkgdCutflow.GetXaxis().SetBinLabel(c, cut)
    hBkgdCutflow.SetDirectory(0) # necessary to keep hist from closing
    hBkgdCutflow.SetDefaultSumw2() # automatically sum w^2 while filling

    for i, fakeTypeId in enumerate(fakeTypes_idDict, start=1):
        hFakeSorting.GetXaxis().SetBinLabel(i, fakeTypes_idDict[fakeTypeId])
    hFakeSorting.SetDirectory(0)
    hFakeSorting.SetDefaultSumw2()

    for plotVar in plotSettings:
        nBins = plotSettings[plotVar][0]
        xMin = plotSettings[plotVar][1]
        xMax = plotSettings[plotVar][2]
        binwidth = (xMax - xMin)/nBins
        hBkgd = TH1D("bkgd_"+name+"_"+plotVar, name, nBins, xMin, xMax)
        hBkgd.SetDirectory(0) # necessary to keep hist from closing
        hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
        hBkgdSubprocessesPlotVarDict[name].update({plotVar:hBkgd})

# ********** Looping over each subprocess. ***********
start_time = time.time()

prevProcess = "" # to determine when you got to the next process
processNum = 0
bkgd_redirector = open("bkgd_fileRedirector")
WIncl_totgenwt = 0
DYIncl_totgenwt = 0
for process in processes:
    hBkgdCutsCountDict.update({process:[0]*nCuts})

for subprocessLine in bkgd_redirector:
    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    if subprocess[0] == "#": continue
    process = subprocessLine[1]
    if not process in processes: continue
    xsec = float(subprocessLine[2])

    # assemble the bkgdNtupleAdr
    bkgdNtupleAdr = myDataDir+"/bkgd/"+process+"/"+subprocess+"/"+subprocess+"_"
    # if testMode: bkgdNtupleAdr += "test_"
    # else: bkgdNtupleAdr += "all_"
    bkgdNtupleAdr += "all_"  
    bkgdNtupleAdr += channel+".root"
    print "Plotting from", bkgdNtupleAdr

    try:
        bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
        tBkgd = bkgdFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+bkgdNtupleAdr+\
                ", skipping\n")
        continue
    
    try: nentries = tBkgd.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+bkgdNtupleAdr+\
                ", skipping\n")
        continue
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+bkgdNtupleAdr+" has no entries!"+\
                " Skipping.\n")
        continue
    print("nentries={0:d}".format(nentries))

    hBkgdGenweights = bkgdFile.Get("genWeights")
    bkgdTotGenweight = hBkgdGenweights.GetSumOfWeights() # tot for this subprocess
    if subprocess == "WJetsToLNu":
        WxGenweightsArr = []
        for i in range(5):
            WxGenweightsArr.append(bkgdFile.Get("W"+str(i)+"genWeights")\
                    .GetSumOfWeights())
    elif subprocess == "DYJetsToLL_M-50":
        DYxGenweightsArr = []
        for i in range(5):
            DYxGenweightsArr.append(bkgdFile.Get("DY"+str(i)+"genWeights").\
                    GetSumOfWeights())

    # don't define this if this subprocess is W/DY incl and we're including W/DYnJets,
    # because then need to separate it out into the partons:
    if (subprocess[:4] != "WJet" or WJetsIncl_only) and \
            (subprocess != "DYJetsToLL_M-50" or DYJetsIncl_only):
        hBkgdPlotVarDict = hBkgdSubprocessesPlotVarDict[subprocess]
        hBkgdCutflow = hBkgdCutflowDict[subprocess] 
        hFakeSorting = hFakeSortingDict[subprocess]

    nMax = nentries
    if testMode: nMax = 10000

    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        if count > nMax : break
        if count % 100000 == 0: print "count =", count
        genwt = event.genWeight
        puwt = event.puWeight
        evtwt = genwt*puwt
        # if evtwt < 0: continue

        nPartons = event.LHE_Njets

        if (not WJetsIncl_only and subprocess[:4] == "WJet") or \
                (not DYJetsIncl_only and subprocess == "DYJetsToLL_M-50"):
            nPartons = event.LHE_Njets
            hBkgdCutflow = hBkgdCutflowDict[subprocess+"_"+str(nPartons)+"Parton"]
            hFakeSorting = hFakeSortingDict[subprocess+"_"+str(nPartons)+"Parton"]
            hBkgdPlotVarDict = hBkgdSubprocessesPlotVarDict[subprocess+"_"+\
                    str(nPartons)+"Parton"]
        else: pass # already defined above
    
        # ********** Additional cuts. ***********

        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index

        l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
        l2Charge = list(getattr(event, l2Flav+"_charge"))[l2Index]
        l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
        l2RelIso = list(getattr(event, l2Flav+"_relIso"))[l2Index]

        if not (event.PV_npvsGood and event.PV_npvs > 1): continue

        if region == "any": pass
        elif region == "A":
            if not isRegionA(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "B":
            if not isRegionB(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "C":
            if not isRegionC(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "D":
            if not isRegionD(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        hBkgdCutflow.Fill(cuts["baseline"], evtwt)

        # if nCuts > cuts["dilepton"]: # not doing this; doing ABCD regions instead
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

        if nCuts > cuts["no3rdlept"]:
            if (event.nExtraMuon > 0 or event.nExtraMuon > 0): continue
            # if event.found3rdLept: continue
            hBkgdCutflow.Fill(cuts["no3rdlept"], evtwt)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hBkgdCutflow.Fill(cuts["nbtag<2"], evtwt)

        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
            hBkgdCutflow.Fill(cuts["MET>80"], evtwt)
        
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
            hBkgdCutflow.Fill(cuts["nJet<4"], evtwt)

        # ********** Sorting fake types. ***********
        # fill twice for each event (once for each lepton)
        l1_genPartFlav = getattr(event, l1Flav+"_genPartFlav")[l1Index]
        l2_genPartFlav = getattr(event, l2Flav+"_genPartFlav")[l2Index]
        hFakeSorting.Fill(fakeTypes_orderDict[l1_genPartFlav])
        hFakeSorting.Fill(fakeTypes_orderDict[l2_genPartFlav])

        # ********** Filling plotVars. ***********
        if event.nJet > 0:
            Jet_pt_arr = np.reshape(event.Jet_pt, 20)
            jMaxPt = 0
            for j in range(event.nJet):
                if Jet_pt_arr[j] > Jet_pt_arr[jMaxPt]: jMaxPt = j

        for plotVar in plotSettings:
            hBkgd = hBkgdPlotVarDict[plotVar]
            nBins = plotSettings[plotVar][0]
            xMin = plotSettings[plotVar][1]
            xMax = plotSettings[plotVar][2]
            binwidth = (xMax - xMin)/nBins

            # Figure out what value to plot.
            # flag to divide by sqrt MET before plotting
            div_sqrt_MET = False
            if "div_sqrt_MET" in plotVar: 
                div_sqrt_MET = True
                plotVar = plotVar[:-len("_div_sqrt_MET")]

            # plotting a jet related var
            if "Jet" in plotVar and plotVar != "nJet":
                if event.nJet == 0: continue # to the next plotVar
                # just want to plot the jet var for the jet with max pt
                if plotVar[:3] == "Jet" and plotVar != "Jet_ht":
                    val = np.reshape(getattr(event, "Jet_"+plotVar[4:]),20)[jMaxPt]
                # dR_lep1/2_jet, Jet_ht
                else: val = getattr(event, plotVar)
            elif plotVar[:4] == "lep1":
                val = list(getattr(event, l1Flav+plotVar[4:]))[l1Index]
            elif plotVar[:4] == "lep2":
                val = list(getattr(event, l2Flav+plotVar[4:]))[l2Index]
            # everything else
            else: val = getattr(event, plotVar)

            # extra processing for these vars
            if div_sqrt_MET:
                val /= sqrt(event.MET_pt)

            # Fill.
            if val <= xMax: hBkgd.Fill(val, evtwt)
            else: hBkgd.Fill(xMax - binwidth/2, evtwt) # overflow 
    
    # ********** Drawing. ***********
    newProcess = False
    if not prevProcess == process:
        processNum += 1
        newProcess = True
    prevProcess = process

    # default:
    norm = xsec*lumi/bkgdTotGenweight

    # special processing for WJets and DYJets inclusive if including W/DYxJets:
    if not WJetsIncl_only and subprocess == "WJetsToLNu":
        WIncl_totgenwt = bkgdTotGenweight # will be used later for WxJets
        WIncl_xsec = xsec
        for i in range(5):
            if i == 0:
                norm = WNJetsXsecs[0]*lumi/WxGenweightsArr[0]
            else:
                norm = lumi/(WIncl_totgenwt/WIncl_xsec + \
                        WxGenweightsArr[i]/(WNJetsXsecs[i]*WJets_kfactor))
            for plotVar in plotSettings:
                hBkgd = hBkgdSubprocessesPlotVarDict[subprocess+"_"+str(i)\
                            +"Parton"][plotVar]
                hBkgd.Scale(norm)
                hBkgd.SetFillColor(processes[process])
                hBkgd.SetLineColor(processes[process])
                hBkgdStacksDict[plotVar].Add(hBkgd)
                if i == 0:
                    legend = legendDict[plotVar]
                    legend.AddEntry(hBkgd, process)
            name = subprocess+"_"+str(i)+"Parton"
            print name
            hBkgdCutflow = hBkgdCutflowDict[name]
            hBkgdCutflow.Scale(norm)
            for c, cut in enumerate(cuts):
                if c >= nCuts: break
                hBkgdCutsCountDict[process][c]+=int(hBkgdCutflow.GetBinContent(c+1))
            hFakeSorting = hFakeSortingDict[name]
            hFakeSorting.Scale(norm)
            hFakeSorting.SetFillColor(processes[process])
            hFakeSorting.SetLineColor(processes[process])
            hFakeSortingStack.Add(hFakeSorting)
            if i == 0:
                legend_fakeSorting.AddEntry(hFakeSorting, process)
    elif not DYJetsIncl_only and subprocess == "DYJetsToLL_M-50":
        DYIncl_totgenwt = bkgdTotGenweight # will be used later for DYxJets
        DYIncl_xsec = xsec
        for i in range(5):
            if i == 0:
                norm = DYNJetsXsecs[0]*lumi/DYxGenweightsArr[0]
            else:
                norm = lumi/(DYIncl_totgenwt/DYIncl_xsec + \
                        DYxGenweightsArr[i]/(DYNJetsXsecs[i]*DYJets_kfactor))
            for plotVar in plotSettings:
                hBkgd = hBkgdSubprocessesPlotVarDict[subprocess+"_"+str(i)\
                            +"Parton"][plotVar]
                hBkgd.Scale(norm)
                hBkgd.SetFillColor(processes[process])
                hBkgd.SetLineColor(processes[process])
                hBkgdStacksDict[plotVar].Add(hBkgd)
                if i == 0:
                    legend = legendDict[plotVar]
                    legend.AddEntry(hBkgd, process)
            name = subprocess+"_"+str(i)+"Parton"
            print name
            hBkgdCutflow = hBkgdCutflowDict[name]
            hBkgdCutflow.Scale(norm)
            for c, cut in enumerate(cuts):
                if c >= nCuts: break
                hBkgdCutsCountDict[process][c] += int(hBkgdCutflow.GetBinContent(c+1))
            hFakeSorting = hFakeSortingDict[name]
            hFakeSorting.Scale(norm)
            hFakeSorting.SetFillColor(processes[process])
            hFakeSorting.SetLineColor(processes[process])
            hFakeSortingStack.Add(hFakeSorting)
            if i == 0:
                legend_fakeSorting.AddEntry(hFakeSorting, process)

    # if running only on W/DYJets incl, then all subprocesses; otherwise all
    # subprocesses other than WJets incl and DYJets incl:
    else:
        # special processing for WnJets and DYnJets (will never happen if running
        # on only W/DYJets incl):
        if subprocess[0] == "W" and subprocess[2:] == "JetsToLNu":
            if WIncl_totgenwt > 0: # if not, then missed WIncl; use the default norm
                norm = lumi/(WIncl_totgenwt/WIncl_xsec + bkgdTotGenweight/\
                        (WNJetsXsecs[int(subprocess[1])-1]*WJets_kfactor))
        elif subprocess[0] == "DY" and subprocess[2:] == "Jets_M-50":
            if DYIncl_totgenwt > 0: # if not, then missed DYIncl; use the default norm
                norm = lumi/(DYIncl_totgenwt/DYIncl_xsec + bkgdTotGenweight/\
                        (DYNJetsXsecs[int(subprocess[1])-1]*DYJets_kfactor))

        # if running only on W/DYJets incl, then all subprocesses; otherwise all
        # subprocesses other than WJets incl and DYJets incl:
        for plotVar in plotSettings:
            hBkgd = hBkgdPlotVarDict[plotVar]
            # hBkgd.Sumw2() # already summed while filling
            hBkgd.Scale(norm)
            hBkgd.SetFillColor(processes[process])
            hBkgd.SetLineColor(processes[process])
            hBkgdStacksDict[plotVar].Add(hBkgd)
            if newProcess:
                legend = legendDict[plotVar]
                legend.AddEntry(hBkgd, process)
        hBkgdCutflow.Scale(norm)
        for c, cut in enumerate(cuts):
            if c >= nCuts: break
            hBkgdCutsCountDict[process][c] += int(hBkgdCutflow.GetBinContent(c+1))
        hFakeSorting.Scale(norm)
        hFakeSorting.SetFillColor(processes[process])
        hFakeSorting.SetLineColor(processes[process])
        hFakeSortingStack.Add(hFakeSorting)
        if newProcess:
            legend_fakeSorting.AddEntry(hFakeSorting, process)

    # all subprocesses:
    bkgdFile.Close()

for plotVar in plotSettings:
    c = canvasDict[plotVar]
    c.cd()
    hBkgdStack = hBkgdStacksDict[plotVar]
    hBkgdStack.Draw("hist")
    unitsLabel = plotSettings[plotVar][3]
    try:
        hBkgdStack.GetXaxis().SetTitle(plotVar+" "+unitsLabel)
        hBkgdStack.GetYaxis().SetTitle("Number of Events, norm to 35921 /pb")
        hBkgdStack.SetMinimum(1)
        hBkgdStack.SetMaximum(10**8)
    except:
        sys.stderr.write("WARNING: no hBkgds were filled!\n")
        continue

c_fakeSort.cd()
hFakeSortingStack.Draw("hist")
hFakeSortingStack.SetMinimum(1)
hFakeSortingStack.SetMaximum(10**8)

#--------------------------------------------------------------------------------#
# *************** Filling each signal in a separate hist ***************
print
print "----------- Plotting from signal. -----------"

sig_redirector = open("sig_fileRedirector")

# hSigCutflowDict maps the subprocess to its cutflow hist
hSigCutflowDict = {} 
# hSigCutsCountDict maps the subprocess (signal type) to an array containing num 
# evts remaining after cut i
hSigCutsCountDict = {}

coloropts = [2,4,3,6,7,9,28,46] # some good colors for lines
markeropts = [1,20,21,22,23] # some good marker styles for lines
linestyleopts = [1,2,3,7,9] # some good styles for lines

# hSigSubprocessesPlotVarDict maps the subprocess (signal type) to another dictionary
# which maps a plotVar to a hist
hSigSubprocessesPlotVarDict = {}

for fileNum, subprocessLine in enumerate(sig_redirector):
    if fileNum + 1 > numSigFiles: break

    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    process = subprocessLine[1]
    xsec = float(subprocessLine[2])
    print subprocess 
    
    # assemble the sigNtupleAdr
    sigNtupleAdr = myDataDir+"/sig/"+process+"/"+subprocess+"/"+subprocess+"_"
    # if testMode: sigNtupleAdr += "test_"
    # else: sigNtupleAdr += "all_"
    sigNtupleAdr += "all_"
    sigNtupleAdr += channel+".root"

    try:
        sigFile = TFile.Open(sigNtupleAdr, "READ")
        tSig = sigFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+sigNtupleAdr+\
                ", skipping\n")
        continue
    try:
        nentries = tSig.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+sigNtupleAdr+\
                ", skipping\n")
        continue
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+sigNtupleAdr+" has no entries!"+\
                " Skipping.\n")
        continue
    print("nentries={0:d}".format(nentries))

    hSigSubprocessesPlotVarDict.update({subprocess:{}})
    hSigPlotVarDict = hSigSubprocessesPlotVarDict[subprocess]
    for plotVar in plotSettings:
        nBins = plotSettings[plotVar][0]
        xMin = plotSettings[plotVar][1]
        xMax = plotSettings[plotVar][2]
        hSig = TH1D("sig_"+subprocess[10:27]+"_"+plotVar, \
                subprocess[10:27], nBins, xMin, xMax)
        hSig.SetDirectory(0)
        hSig.SetDefaultSumw2() # automatically sum w^2 while filling
        hSigPlotVarDict.update({plotVar:hSig})

    hSigCutflow = TH1D("sig_" + subprocess, "sig_" + subprocess[10:27], nCuts, 0, \
            nCuts)
    hSigCutflow.SetDirectory(0)
    hSigCutflow.SetDefaultSumw2() # automatically sum w^2 while filling
    hSigCutflowDict.update({subprocess:hSigCutflow})
    hSigCutsCountDict.update({subprocess:[0]*nCuts})
    for c, cut in enumerate(cuts, start=1):
        if c>nCuts: break
        hSigCutflow.GetXaxis().SetBinLabel(c, cut)

    hSigGenweights = sigFile.Get("genWeights")
    sigTotGenweight = hSigGenweights.GetSumOfWeights()

    nMax = nentries
    # if testMode: nMax = 1000

    # ********** Looping over events. ***********
    for count, event in enumerate(tSig):
        if count > nMax : break
        if count % 100000 == 0: print "count =", count
        genwt = event.genWeight
        puwt = event.puWeight
        evtwt = genwt*puwt

        # ********** Additional cuts. ***********

        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
    
        l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
        l2Charge = list(getattr(event, l2Flav+"_charge"))[l2Index]
        l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
        l2RelIso = list(getattr(event, l2Flav+"_relIso"))[l2Index]

        if not (event.PV_npvsGood and event.PV_npvs > 1): continue

        if region == "any": pass
        elif region == "A":
            if not isRegionA(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "B":
            if not isRegionB(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "C":
            if not isRegionC(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "D":
            if not isRegionD(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        hSigCutflow.Fill(cuts["baseline"], evtwt)

        # if nCuts > cuts["dilepton"]: # not doing this; doing ABCD regions instead
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

        if nCuts > cuts["no3rdlept"]:
            if (event.nExtraMuon > 0 or event.nExtraMuon > 0): continue
            # if event.found3rdLept: continue
            hSigCutflow.Fill(cuts["no3rdlept"], evtwt)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hSigCutflow.Fill(cuts["nbtag<2"], evtwt)

        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
            hSigCutflow.Fill(cuts["MET>80"], evtwt)
        
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
            hSigCutflow.Fill(cuts["nJet<4"], evtwt)

        # ********** Filling. ***********
        if event.nJet > 0:
            Jet_pt_arr = np.reshape(event.Jet_pt, 20)
            jMaxPt = 0
            for j in range(event.nJet):
                if Jet_pt_arr[j] > Jet_pt_arr[jMaxPt]: jMaxPt = j

        for plotVar in plotSettings:
            hSig = hSigPlotVarDict[plotVar]
            nBins = plotSettings[plotVar][0]
            xMin = plotSettings[plotVar][1]
            xMax = plotSettings[plotVar][2]
            binwidth = (xMax - xMin)/nBins

            # Figure out what value to plot.

            # flag to divide by sqrt MET before plotting
            div_sqrt_MET = False
            if "div_sqrt_MET" in plotVar: 
                div_sqrt_MET = True
                plotVar = plotVar[:-13]

            # plotting a jet related var
            if "Jet" in plotVar and plotVar != "nJet":
                if event.nJet == 0: continue # to the next plotVar
                # just want to plot the jet var for the jet with max pt
                if plotVar[:3] == "Jet" and plotVar != "Jet_ht":
                    val = np.reshape(getattr(event, "Jet_"+plotVar[4:]),20)[jMaxPt]
                # dR_lep1/2_jet, Jet_ht
                else: val = getattr(event, plotVar)
            elif plotVar[:4] == "lep1":
                val = list(getattr(event, l1Flav+plotVar[4:]))[l1Index]
            elif plotVar[:4] == "lep2":
                val = list(getattr(event, l2Flav+plotVar[4:]))[l2Index]
            # everything else
            else: val = getattr(event, plotVar)

            # extra processing for these vars
            if div_sqrt_MET:
                val /= sqrt(event.MET_pt)

            # Fill.
            if val <= xMax: hSig.Fill(val, evtwt)
            else: hSig.Fill(xMax - binwidth/2, evtwt) # overflow 

    # ********** Drawing. ***********
    hcolor = coloropts[fileNum % len(coloropts)]
    hmarkerstyle = markeropts[(fileNum/len(coloropts)) % len(markeropts)]

    norm = xsec*lumi/sigTotGenweight

    for plotVar in plotSettings:
        c = canvasDict[plotVar]
        c.cd()
        hSig = hSigPlotVarDict[plotVar]
        hSig.SetLineColor(hcolor) 
        hSig.SetMarkerStyle(hmarkerstyle)
        hSig.SetMarkerColor(hcolor)
        hlinestyle = linestyleopts[(fileNum/len(coloropts)/len(markeropts)) % \
                len(linestyleopts)]
        hSig.SetLineStyle(hlinestyle)

        hSig.Scale(norm)
        legend = legendDict[plotVar]
        legend.AddEntry(hSig, hSig.GetTitle())
        hSig.Draw("hist same") # same pad

    hSigCutflow.Scale(norm)
    for c, cut in enumerate(cuts):
        if c >= nCuts: break
        hSigCutsCountDict[subprocess][c] = int(hSigCutflow.GetBinContent(c+1))

    sigFile.Close()

#--------------------------------------------------------------------------------#
# *************** Filling data in a separate hist  ************
print
print "----------- Plotting from data. -----------"
data_redirector = open("data_fileRedirector")

# hDataPlotVarDict maps each plotVar to a hist (which contains all the data from the
# process)
hDataPlotVarDict = {}
for plotVar in plotSettings:
    nBins = plotSettings[plotVar][0]
    xMin = plotSettings[plotVar][1]
    xMax = plotSettings[plotVar][2]
    hData = TH1D("data_"+plotVar, "data", nBins, xMin, xMax)
    hData.SetDirectory(0)
    hData.SetDefaultSumw2() # automatically sum w^2 while filling
    hDataPlotVarDict.update({plotVar:hData})

hDataCutflow = TH1D("data", "data", nCuts, 0, nCuts)
hDataCutflow.SetDirectory(0)
hDataCutflow.SetDefaultSumw2() # automatically sum w^2 while filling
hDataCutCountArr = [0]*nCuts

for fileNum, subprocessLine in enumerate(data_redirector):
    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    if subprocess[0] == "#": continue
    process = subprocessLine[1]
    if process != dataProcess: continue

    # assemble the dataNtupleAdr
    dataNtupleAdr = myDataDir+"/data/"+process+"/"+subprocess+"/"+subprocess+"_"
    # if testMode: dataNtupleAdr += "test_"
    # else: dataNtupleAdr += "all_"
    dataNtupleAdr += "all_"
    dataNtupleAdr += channel+".root"
    print dataNtupleAdr
    
    try:
        dataFile = TFile.Open(dataNtupleAdr, "READ")
        tData = dataFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+dataNtupleAdr+\
                ", skipping\n")
        continue
    try:
        nentries = tData.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+dataNtupleAdr+\
                ", skipping\n")
        continue
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+dataNtupleAdr+" has no entries!"+\
                " Skipping.\n")
        continue
    print("nentries={0:d}".format(nentries))
    
    nMax = nentries
    if testMode: nMax = 1000

    # ********** Looping over events. ***********
    for count, event in enumerate(tData):
        if count > nMax : break
        if count % 100000 == 0: print "count =", count
        evtwt = 1 # just in case I forgot to change it somewhere
    
        # ********** Additional cuts. ***********
    
        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
    
        l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
        l2Charge = list(getattr(event, l2Flav+"_charge"))[l2Index]
        l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
        l2RelIso = list(getattr(event, l2Flav+"_relIso"))[l2Index]

        if not (event.PV_npvsGood and event.PV_npvs > 1): continue

        if region == "any": pass
        elif region == "A":
            if not isRegionA(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "B":
            if not isRegionB(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "C":
            if not isRegionC(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        elif region == "D":
            if not isRegionD(l1Charge, l2Charge, l1RelIso, l2RelIso, \
                    findingSameFlavor): continue
        hDataCutflow.Fill(cuts["baseline"], 1)

        # if nCuts > cuts["dilepton"]: # not doing this; doing ABCD regions instead
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

        if nCuts > cuts["no3rdlept"]:
            if (event.nExtraMuon > 0 or event.nExtraMuon > 0): continue
            # if event.found3rdLept: continue
            hDataCutflow.Fill(cuts["no3rdlept"], 1)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hDataCutflow.Fill(cuts["nbtag<2"], 1)

        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
            hDataCutflow.Fill(cuts["MET>80"], 1)

        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
            hDataCutflow.Fill(cuts["nJet<4"], 1)

        # ********** Filling. ***********
        if event.nJet > 0:
            Jet_pt_arr = np.reshape(event.Jet_pt, 20)
            jMaxPt = 0
            for j in range(event.nJet):
                if Jet_pt_arr[j] > Jet_pt_arr[jMaxPt]: jMaxPt = j

        for plotVar in plotSettings:
            hData = hDataPlotVarDict[plotVar]
            nBins = plotSettings[plotVar][0]
            xMin = plotSettings[plotVar][1]
            xMax = plotSettings[plotVar][2]
            binwidth = (xMax - xMin)/nBins

            # Figure out what value to plot.
            # flag to divide by sqrt MET before plotting
            div_sqrt_MET = False
            if "div_sqrt_MET" in plotVar: 
                div_sqrt_MET = True
                plotVar = plotVar[:-len("_div_sqrt_MET")]

            # plotting a jet related var
            if "Jet" in plotVar and plotVar != "nJet":
                if event.nJet == 0: continue # to the next plotVar
                # just want to plot the jet var for the jet with max pt
                if plotVar[:3] == "Jet" and plotVar != "Jet_ht":
                    val = np.reshape(getattr(event, "Jet_"+plotVar[4:]),20)[jMaxPt]
                # dR_lep1/2_jet, Jet_ht
                else: val = getattr(event, plotVar)
            elif plotVar[:4] == "lep1":
                val = list(getattr(event, l1Flav+plotVar[4:]))[l1Index]
            elif plotVar[:4] == "lep2":
                val = list(getattr(event, l2Flav+plotVar[4:]))[l2Index]
            # everything else
            else: val = getattr(event, plotVar)
 
            # extra processing for these vars
            if div_sqrt_MET:
                val /= sqrt(event.MET_pt)

            # Fill.
            if val <= xMax: hData.Fill(val, 1)
            else: hData.Fill(xMax - binwidth/2, 1) # overflow 

    dataFile.Close()

for c, cut in enumerate(cuts):
    if c >= nCuts: break
    hDataCutCountArr[c] += int(hDataCutflow.GetBinContent(c+1))

# ********** Drawing. ***********
hcolor = 1 # black
hmarkerstyle = 3 # asterisk (to match with the *H draw option)

for plotVar in plotSettings:
    c = canvasDict[plotVar]
    c.cd()
    hData = hDataPlotVarDict[plotVar]
    hData.SetLineColor(hcolor) 
    hData.SetMarkerStyle(hmarkerstyle)
    hData.SetMarkerColor(hcolor)

    legend = legendDict[plotVar]
    legend.AddEntry(hData, hData.GetTitle())
    hData.Draw("* hist same") # same pad

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print
print int(time.time()-start_time), "secs of processing."
print "Drawing."

for plotVar in plotSettings:
    c = canvasDict[plotVar]
    c.cd()
    legend = legendDict[plotVar]
    legend.SetTextSize(0.017)
    legend.SetNColumns(3)
    legend.Draw("same")
    c.SetLogy()
    c.Update()
c_fakeSort.cd()
legend_fakeSorting.SetTextSize(0.017)
legend_fakeSorting.SetNColumns(3)
legend_fakeSorting.Draw("same")
c_fakeSort.SetLogy()
c_fakeSort.Update()

# making cutflow pandas dataframe
statsStack = {}
thisCuts = list(cuts)[:nCuts]
statsStack.update({"Cut_name":thisCuts})
statsNamesList = []
for process in processes:
    statsStack.update({process:hBkgdCutsCountDict[process]})
    statsNamesList.append(process)
for subprocess in hSigCutflowDict:
    sigName = hSigCutflowDict[subprocess].GetTitle()[4:]
    statsStack.update({sigName:hSigCutsCountDict[subprocess]})
    statsNamesList.append(sigName)
statsStack.update({"data":hDataCutCountArr})
statsNamesList.append("data")
statsDF = pd.DataFrame(statsStack)
statsDF.set_index('Cut_name', inplace = True) # keep the Cut_name as first column
statsDF = statsDF[statsNamesList] # reorder columns
if not os.path.exists(statsDir): os.makedirs(statsDir)
statsFileName = statsDir+"/cutflow_stats_"+channel+"_"+region
# if experimental: statsFileName += "_experimental"

if displayMode:
    print statsDF
    print "Done. Press enter to finish (plots not saved)."
    raw_input()
else:
    gSystem.ProcessEvents()
    if not os.path.exists(imgDir): os.makedirs(imgDir) 
    outHistFileAdr = imgDir+"/QCDMC_plot1D_"
    if testMode: outHistFileAdr += "test_"
    else: outHistFileAdr += "all_"
    outHistFileAdr += channel+"_"+lastcut+"_"+region+".root"
    outHistFile = TFile(outHistFileAdr, "recreate")
    for plotVar in plotSettings:
        canvasDict[plotVar].Write()
        for subprocess in hBkgdSubprocessesPlotVarDict:
            hBkgdSubprocessesPlotVarDict[subprocess][plotVar].Write()
        for subprocess in hSigSubprocessesPlotVarDict:
            hSigSubprocessesPlotVarDict[subprocess][plotVar].Write()
        hDataPlotVarDict[plotVar].Write()
    c_fakeSort.Write()
    for subprocess in hFakeSortingDict:
        hFakeSortingDict[subprocess].Write()
    outHistFile.Close()
    print "Saved hists in", outHistFileAdr

    if lastcut == "nJet<4": # don't save the cutflow unless made all cuts
        statsFile = open(statsFileName+".txt", "w")
        statsFile.write(statsDF.to_string()+"\n")
        statsFile.close()
        print "Saved file", statsFileName+".txt"
        
        statsDF.to_hdf(statsFileName+".hdf", key="statsDF")
        print "Saved file", statsFileName+".hdf"
    else: print statsDF

    print "Done."
