#!/usr/bin/env python

# NOTE: NEEDS 2 CMD LINE ARGS with values:
# testMode {test, all}, channel {mumu, elel, muel}
#
# Implements additional cuts and then writes the counts of each sig/bkgd type
# after each cut to an output file of the form cutflow_stats_[channel].txt in
# the plots folder.
#
# Uses the root files outputted by hadding the output from makeNtuple.py
# Uses xsec info from bkgd_fileRedirector
# Uses xsec info from sig_fileRedirector

print "Importing modules."
import sys, os
from ROOT import TFile, TTree, TH1D, TCanvas, TImage, TLegend, TText, THStack
from ROOT import gSystem, gStyle, gROOT, kTRUE
from stopSelection import deltaR
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from collections import OrderedDict
import numpy as np
import time
print "Beginning execution of", sys.argv

assert len(sys.argv) == 3, "needs 2 command line args: testMode{0,1}, channel {mumu, elel, muel}"

if sys.argv[1] == "test": testMode = True
elif sys.argv[1] == "all": testMode = False
else: assert False, "invalid test mode, need {test, all}"

# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
if sys.argv[2] == "mumu":
    findingSameFlavor = True
    muPreference = True
    l1Flav = "Muon"
    l2Flav = "Muon"
    dataProcess = "DoubleMuon"
elif sys.argv[2] == "elel":
    findingSameFlavor = True
    muPreference = False
    l1Flav = "Electron"
    l2Flav = "Electron"
    dataProcess = "DoubleEG"
elif sys.argv[2] == "muel":
    findingSameFlavor = False
    muPreference = False
    l1Flav = "Muon"
    l2Flav = "Electron"
    dataProcess = "MuonEG"
else: assert False, "invalid channel, need {mumu, elel, muel}"
channelName = l1Flav[:2] + l2Flav[:2]

cuts = OrderedDict([("nocut",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3), \
        ("MET>80",4),("nJet<4",5)])
# if experimental:
#     cuts = OrderedDict([("nocut",0), ("dilepton",1), ("no3rdlept",2), ("nbtag<2",3),\
#             ("MET>20",4), ("MET>50",5), ("MET>80",6), ("MET>110",7), ("nJet<4",8)])
nCuts = len(cuts)


# bkgd process name : color for plotting
processes = OrderedDict([("W-Jets",38), ("Drell-Yan",46), ("TTBar",30), \
        ("Diboson",41), ("Single-Top",40), ("TT+X",7)])

myDataDir = "/eos/user/c/cmiao/private/myDataSusy/Run2/"
# number of files to process
numBkgdFiles = float("inf")  # note: must loop over all files to have correct xsec
numSigFiles = 3 # max 25

lumi = 35921 # 2016 lumi in /pb
WJets_kfactor = 1.221
DYJets_kfactor = 1.1637

gStyle.SetOptStat(0) # don't show any stats

#--------------------------------------------------------------------------------#
start_time = time.time()
# *************** Filling bkgd data summed together  ************
print
print "----------- Plotting from background. -----------"

# hBkgdDict maps every subprocess to an hBkgd which contains data from all the 
# ntuples for that subprocess.
hBkgdDict = {}
# hBkgdCutsCountDict maps every process to an array of size nCuts that keeps track
# of the num evts remaining after each cut for that process
hBkgdCutsCountDict = {}
WNJetsXsecs = [47297.3] # first entry: W0Jets xsec
DYNJetsXsecs = [4263.5]  # first entry: DY0Jets xsec
with open("bkgd_fileRedirector") as bkgdSubprocessesListFile:
    for subprocessLine in bkgdSubprocessesListFile:
        subprocessLine = subprocessLine.rstrip('\n').split(" ")
        subprocess = subprocessLine[0]
        if subprocess[0] == "#": continue # problematic input files

        if subprocess[0] == "W" and subprocess[2:] == "JetsToLNu":
            WNJetsXsecs.append(float(subprocessLine[2]))
        if subprocess[0] == "DY" and subprocess[2:] == "JetsToLL_M-50":
            DYNJetsXsecs.append(float(subprocessLine[2]))

        if subprocess == "WJetsToLNu" or subprocess == "DYJetsToLL_M-50":
            for i in range(5):
                name = subprocess+"_"+str(i)+"Parton"
                hBkgd = TH1D("cutflow_"+name+"_bkgd", \
                        "cutflow_"+name+"_bkgd", nCuts, 0, nCuts)
                for i, cut in enumerate(cuts, start=1):
                    if i>nCuts: break
                    hBkgd.GetXaxis().SetBinLabel(i, cut)
                hBkgd.SetDirectory(0) # necessary to keep hist from closing
                hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
                hBkgdDict.update({name:hBkgd})
        else:
            hBkgd = TH1D("cutflow_"+subprocess+"_bkgd", \
                    "cutflow_"+subprocess+"_bkgd", nCuts, 0, nCuts)
            for i, cut in enumerate(cuts, start=1):
                if i>nCuts: break
                hBkgd.GetXaxis().SetBinLabel(i, cut)
            hBkgd.SetDirectory(0) # necessary to keep hist from closing
            hBkgd.SetDefaultSumw2() # automatically sum w^2 while filling
            hBkgdDict.update({subprocess:hBkgd})
for process in processes:
    hBkgdCutsCountDict.update({process:[0]*nCuts})
title = "cutflow ("+channelName+")"
hBkgdStack = THStack("cutflow_bkgdStack", title)

# ********** Looping over each subprocess. ***********
prevProcess = "" # to determine when you got to the next process
processNum = 0
bkgdSubprocessesListFile = open("bkgd_fileRedirector")
for subprocessLine in bkgdSubprocessesListFile:
    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    if subprocess[0] == "#": continue
    process = subprocessLine[1]
    if not process in processes: continue
    xsec = float(subprocessLine[2])

    # assemble the bkgdNtupleAdr
    bkgdNtupleAdr = myDataDir+"bkgd/"+process+"/"+subprocess+"/"+subprocess+"_"
    if testMode: bkgdNtupleAdr += "test_"
    else: bkgdNtupleAdr += "all_"
    bkgdNtupleAdr += channelName+".root"
    print "Plotting from", bkgdNtupleAdr

    try:
        bkgdFile = TFile.Open(bkgdNtupleAdr, "READ")
        tBkgd = bkgdFile.Get("Events")
    except:
        sys.stderr.write("WARNING: nonexistent or corrupted file "+bkgdNtupleAdr+\
                ", skipping\n")
        continue
    
    try:
        nentries = tBkgd.GetEntries()
    except:
        sys.stderr.write("WARNING: unable to get entries from "+bkgdNtupleAdr+\
                ", skipping\n")
        continue
    print("nentries={0:d}".format(nentries))
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+bkgdNtupleAdr+" has no entries!"+\
                " Skipping\n")
        continue

    if subprocess[:4] != "WJet" and subprocess != "DYJetsToLL_M-50":
        hBkgd = hBkgdDict[subprocess] 
    
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
    
    nMax = nentries
    if testMode: nMax = 10000

    # ********** Looping over events. ***********
    for count, event in enumerate(tBkgd):
        if count > nMax : break
        if count % 100000 == 0: print "count =", count
        genwt = event.genWeight
        puwt = event.puWeight
        evtwt = genwt*puwt

        if subprocess[:4] == "WJet" or subprocess == "DYJetsToLL_M-50":
            nPartons = event.LHE_Njets
            if nPartons < 1 or nPartons > 4: continue
            hBkgd = hBkgdDict[subprocess+"_"+str(nPartons)+"Parton"]
    
        # ********** Additional cuts. ***********
        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
        hBkgd.Fill(cuts["nocut"], evtwt)
    
        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if list(getattr(event, l1Flav+"_charge"))[l1Index] * \
                    list(getattr(event, l2Flav+"_charge"))[l2Index] >= 0: continue
            if findingSameFlavor:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.1: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.1: continue
            else:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.2: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.2: continue
            hBkgd.Fill(cuts["dilepton"], evtwt)
    
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hBkgd.Fill(cuts["deltaR(ll)>0.3"], evtwt)
    
        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hBkgd.Fill(cuts["no3rdlept"], evtwt)
    
        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hBkgd.Fill(cuts["nbtag<2"], evtwt)

        # if experimental:
        #     if nCuts > cuts["MET>20"]:
        #         if event.MET_pt < 20: continue
        #         hBkgd.Fill(cuts["MET>20"], evtwt)
        #     if nCuts > cuts["MET>50"]:
        #         if event.MET_pt < 50: continue
        #         hBkgd.Fill(cuts["MET>50"], evtwt)

        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
            hBkgd.Fill(cuts["MET>80"], evtwt)

        # if experimental:
        #     if nCuts > cuts["MET>110"]:
        #         if event.MET_pt < 110: continue
        #         hBkgd.Fill(cuts["MET>110"], evtwt)
            
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
            hBkgd.Fill(cuts["nJet<4"], evtwt)
    
    # default:
    norm = xsec*lumi/bkgdTotGenweight

    # special processing for WJets and DYJets inclusive:
    if subprocess == "WJetsToLNu":
        WIncl_totgenwt = bkgdTotGenweight # will be used later for WxJets
        WIncl_xsec = xsec
        for i in range(5):
            if i == 0:
                norm = lumi/(WxGenweightsArr[i]/(WNJetsXsecs[i]))
            else:
                norm = lumi/(WIncl_totgenwt/WIncl_xsec + \
                        WxGenweightsArr[i]/(WNJetsXsecs[i]*WJets_kfactor))
            hBkgd = hBkgdDict[subprocess+"_"+str(i)+"Parton"]
            hBkgd.Scale(norm)
            print subprocess+"_"+str(i)+"Parton"
            for i, cut in enumerate(cuts):
                hBkgdCutsCountDict[process][i] += int(hBkgd.GetBinContent(i+1))
                print i, hBkgdCutsCountDict[process][i]
    elif subprocess == "DYJetsToLL_M-50":
        DYIncl_totgenwt = bkgdTotGenweight # will be used later for DYxJets
        DYIncl_xsec = xsec
        for i in range(5):
            if i == 0:
                norm = lumi/(DYxGenweightsArr[i]/(DYNJetsXsecs[i]))
            else:
                norm = lumi/(DYIncl_totgenwt/DYIncl_xsec + \
                        DYxGenweightsArr[i]/(DYNJetsXsecs[i]*DYJets_kfactor))
            hBkgd = hBkgdDict[subprocess+"_"+str(i)+"Parton"]
            hBkgd.Scale(norm)
            print subprocess+"_"+str(i)+"Parton"
            for i, cut in enumerate(cuts):
                hBkgdCutsCountDict[process][i] += int(hBkgd.GetBinContent(i+1))
                print i, hBkgdCutsCountDict[process][i]

    # all subprocesses other than WJets incl and DYJets incl:
    else:
        # special processing for WnJets and DYnJets:
        if subprocess[0] == "W" and subprocess[2:] == "JetsToLNU":
            if WIncl_totgenwt > 0: # if missed WIncl, use the default norm
                norm = lumi/(WIncl_totgenwt/WIncl_xsec + bkgdTotGenweight/\
                        (WNJetsXsecs[int(subprocess[1])-1]*WJets_kfactor))
        elif subprocess[0] == "DY" and subprocess[2:] == "Jets_M-50":
            if DYIncl_totgenwt > 0: # if missed DYIncl, use the default norm
                norm = lumi/(DYIncl_totgenwt/DYIncl_xsec + bkgdTotGenweight/\
                        (DYNJetsXsecs[int(subprocess[1])-1]*DYJets_kfactor))
        # all subprocesses other than WJets incl and DYJets incl:
        hBkgd.Scale(norm)
        for i, cut in enumerate(cuts):
            hBkgdCutsCountDict[process][i] += int(hBkgd.GetBinContent(i+1))
            print i, hBkgdCutsCountDict[process][i]

    # all subprocesses:
    bkgdFile.Close()


#--------------------------------------------------------------------------------#
# *************** Filling each signal in a separate hist  ************
print
print "----------- Plotting from signal. -----------"

sig_redirector = open("sig_fileRedirector")

# hSigDict maps the subprocess to its cutflow hist
hSigDict = {} 

# hSigCutsCountDict maps the subprocess (signal type) to an array containing num 
# evts remaining after cut i
hSigCutsCountDict = {}

for fileNum, subprocessLine in enumerate(sig_redirector):
    if fileNum + 1 > numSigFiles: break

    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    process = subprocessLine[1]
    xsec = float(subprocessLine[2])
    print subprocess 
    
    # assemble the sigNtupleAdr
    sigNtupleAdr = myDataDir+"sig/"+process+"/"+subprocess+"/"+subprocess+"_"
    if testMode: sigNtupleAdr += "test_"
    else: sigNtupleAdr += "all_"
    sigNtupleAdr += channelName+".root"

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
    print("nentries={0:d}".format(nentries))
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+sigNtupleAdr+" has no entries!"+\
                " Skipping\n")
        continue

    hSig = TH1D("sig_" + subprocess, "sig_" + subprocess[10:27], nCuts, 0, nCuts)
    hSig.SetDirectory(0)
    hSig.SetDefaultSumw2() # automatically sum w^2 while filling
    hSigDict.update({subprocess:hSig})
    hSigCutsCountDict.update({subprocess:[0]*nCuts})

    for i, cut in enumerate(cuts, start=1):
        if i>nCuts: break
        hSig.GetXaxis().SetBinLabel(i, cut)

    hSigGenweights = sigFile.Get("genWeights")
    sigTotGenweight = hSigGenweights.GetSumOfWeights()

    # ********** Looping over events. ***********
    for count, event in enumerate(tSig):
        if count % 100000 == 0: print("count={0:d}".format(count))
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
        hSig.Fill(cuts["nocut"], evtwt)

        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if list(getattr(event, l1Flav+"_charge"))[l1Index] * \
                    list(getattr(event, l2Flav+"_charge"))[l2Index] >= 0: continue
            if findingSameFlavor:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.1: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.1: continue
            else:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.2: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.2: continue
            hSig.Fill(cuts["dilepton"], evtwt)


        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        # hSig.Fill(cuts["deltaR(ll)>0.3"], evtwt)

        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hSig.Fill(cuts["no3rdlept"], evtwt)

        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hSig.Fill(cuts["nbtag<2"], evtwt)

        # if experimental:
        #     if nCuts > cuts["MET>20"]:
        #         if event.MET_pt < 20: continue
        #         hSig.Fill(cuts["MET>20"], evtwt)
        #     if nCuts > cuts["MET>50"]:
        #         if event.MET_pt < 50: continue
        #         hSig.Fill(cuts["MET>50"], evtwt)

        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
            hSig.Fill(cuts["MET>80"], evtwt)

        # if experimental:
        #     if nCuts > cuts["MET>110"]:
        #         if event.MET_pt < 110: continue
        #         hSig.Fill(cuts["MET>110"], evtwt)
        
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
            hSig.Fill(cuts["nJet<4"], evtwt)

    hSig.Scale(xsec*lumi/sigTotGenweight)

    for i, cut in enumerate(cuts):
        hSigCutsCountDict[subprocess][i] = int(hSig.GetBinContent(i+1))

    sigFile.Close()

#--------------------------------------------------------------------------------#
# *************** Filling data in a separate hist  ************
print
print "----------- Plotting from data. -----------"
data_redirector = open("data_fileRedirector")

hData = TH1D("data", "data", nCuts, 0, nCuts)
hData.SetDirectory(0)
hData.SetDefaultSumw2() # automatically sum w^2 while filling

hDataCutCountArr = [0]*nCuts

for fileNum, subprocessLine in enumerate(data_redirector):
    subprocessLine = subprocessLine.rstrip('\n').split(" ")
    subprocess = subprocessLine[0]
    if subprocess[0] == "#": continue
    process = subprocessLine[1]
    if process != dataProcess: continue

    # assemble the dataNtupleAdr
    dataNtupleAdr = myDataDir+"data/"+process+"/"+subprocess+"/"+subprocess+"_"
    if testMode: dataNtupleAdr += "test_"
    else: dataNtupleAdr += "all_"
    dataNtupleAdr += channelName+".root"
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
    print("nentries={0:d}".format(nentries))
    if nentries == 0:
        sys.stderr.write("WARNING: tree in "+dataNtupleAdr+" has no entries!"+\
                " Skipping\n")
        continue
    
    # ********** Looping over events. ***********
    for count, event in enumerate(tData):
        if count % 100000 == 0: print "count =", count
    
        # ********** Additional cuts. ***********
    
        # if findingSameFlavor, l1/l2Flav set at runtime
        if not findingSameFlavor: 
            if event.lep1_isMu: l1Flav = "Muon"
            else: l1Flav = "Electron"
            if event.lep2_isMu: l1Flav = "Muon"
            else: l2Flav = "Electron"
        l1Index = event.lep1_index
        l2Index = event.lep2_index
        hData.Fill(cuts["nocut"], evtwt)
    
        if nCuts > cuts["dilepton"]: # currently just tighter relIso cuts
            if list(getattr(event, l1Flav+"_charge"))[l1Index] * \
                    list(getattr(event, l2Flav+"_charge"))[l2Index] >= 0: continue
            if findingSameFlavor:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.1: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.1: continue
            else:
                if list(getattr(event, l1Flav+"_relIso"))[l1Index] >= 0.2: continue
                if list(getattr(event, l2Flav+"_relIso"))[l2Index] >= 0.2: continue
            hData.Fill(cuts["dilepton"], evtwt)
    
        # if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue

        if nCuts > cuts["no3rdlept"]:
            if event.found3rdLept: continue
            hData.Fill(cuts["no3rdlept"], evtwt)
    
        if nCuts > cuts["nbtag<2"]:
            if event.nbtag > 1: continue
            hData.Fill(cuts["nbtag<2"], evtwt)
    
        if nCuts > cuts["MET>80"]:
            if event.MET_pt < 80: continue
            hData.Fill(cuts["MET>80"], evtwt)
        
        if nCuts > cuts["nJet<4"]:
            if event.nJet >= 4: continue
            hData.Fill(cuts["nJet<4"], evtwt)
    
    for i, cut in enumerate(cuts):
        hDataCutCountArr[i] += int(hData.GetBinContent(i+1))

    dataFile.Close()

#--------------------------------------------------------------------------------#
# *************** Wrap up. *******************
print
print int(time.time()-start_time), "secs of processing."

statsHeader = "Cut_name      "
statsStack = np.array(cuts.keys()).reshape(nCuts, 1)
for process in processes:
    statsHeader +=  process + ''.join([' ']*(17-len(process)))
    statsStack = np.append(statsStack, \
            np.array(hBkgdCutsCountDict[process]).reshape(nCuts, 1), axis=1)
for subprocess in hSigDict:
    sigName = hSigDict[subprocess].GetTitle()[4:]
    statsHeader += sigName + ''.join([' ']*(17-len(sigName)))
    statsStack = np.append(statsStack, \
            np.array(hSigCutsCountDict[subprocess]).reshape(nCuts, 1), axis=1)
statsHeader += "data" 
statsStack = np.append(statsStack,np.array(hDataCutCountArr).reshape(nCuts,1),axis=1)
statsDir = "/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/"+\
        "plots/Run2/v1/cutflow_stats"
if not os.path.exists(statsDir): os.makedirs(statsDir)
statsFileName = statsDir+"/cutflow_stats_"+channelName
# if experimental: statsFileName += "_experimental"
statsFileName += ".txt"
np.savetxt(statsFileName, statsStack, delimiter='    ', header=statsHeader, \
        fmt='%-13s')
print "Saved file", statsFileName

print "Done."

