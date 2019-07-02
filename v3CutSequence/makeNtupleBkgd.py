#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, findingSameFlavor, muPreference, process
# For each process listed in bkgd_files, outputs an ntuple located in ../myData/, 
# each containing 1 tree which contains the events from all ntuples listed in 
# the corresponding bkgdProcesses file(s) that have survived loose dilepton 
# selection cuts.
# The tree in each outputted ntuple, tBkgd, has branches for the same variables 
# as the tSig{i} outputted by makeNtupleSigs.py.
# Uses bkgd_files for xsec and process lists
# Uses files in bkgdProcesses dir for ntuple lists.
# Process options: TT+X Diboson W-Jets Drell-Yan Single-Top

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
import numpy as np
from math import sqrt, cos
from array import array

assert len(sys.argv) == 5, "need 4 command line args: testMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}, process"

# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
findingSameFlavor = bool(int(sys.argv[2]))
print "Finding same flavor:", findingSameFlavor
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
muPreference = bool(int(sys.argv[3]))
print "Mu preference:", muPreference

# name of the eventual process that this output root file will be hstacked into
process = sys.argv[4]
processes = {"TT+X", "Diboson", "W-Jets", "Drell-Yan", "Single-Top"}
assert process in processes, "invalid process %s" % process
print "Process:", process

if findingSameFlavor:
    if muPreference: 
        l1Flav = "muon"
        l2Flav = "muon"
    else: 
        l1Flav = "electron"
        l2Flav = "electron"
else: 
    # these 2 lines just matter for creating the outFile name; actual 
    # selection of leading/trailing flavors occurs when looping over events:
    l1Flav = "muon"
    l2Flav = "electron"
channelName = l1Flav[:2] + l2Flav[:2]

# number of ntuples to loop on for each background process
numBkgdFiles = float("inf")  # note: must loop over all files to have correct xsec
if testMode: 
    numBkgdFiles = 2 

outDir = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/"

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
muon_count = array('i',[0])
muon_pt = np.zeros(20, dtype=np.float32)
muon_eta = np.zeros(20, dtype=np.float32)
muon_phi = np.zeros(20, dtype=np.float32)
muon_relIso = np.zeros(20, dtype=np.float32)
muon_charge = np.zeros(20, dtype=np.float32)
muon_mt = np.zeros(20, dtype=np.float32)
electron_count = array('i',[0])
electron_pt = np.zeros(20, dtype=np.float32)
electron_eta = np.zeros(20, dtype=np.float32)
electron_phi = np.zeros(20, dtype=np.float32)
electron_relIso = np.zeros(20, dtype=np.float32)
electron_charge = np.zeros(20, dtype=np.float32)
electron_mt = np.zeros(20, dtype=np.float32)
found3rdLept = array('i',[0])
lep1_isMu = array('i',[0]) # True (1) if muon, False (0) if electron
lep1_index = array('i',[0])
lep2_isMu = array('i',[0]) # True (1) if muon, False (0) if electron
lep2_index = array('i',[0])
jet_pt = np.zeros(20, dtype=np.float32)
jet_eta = np.zeros(20, dtype=np.float32)
jet_phi = np.zeros(20, dtype=np.float32)
jet_ht = array('f',[0.])
# jet_flavour = array('f',[0])
njets = array('i',[0])
nbtag = array('i',[0])
nbtagLoose = array('i',[0])
nbtagTight = array('i',[0])
met_pt = array('f',[0.])
met_phi = array('f',[0.])
genweight = array('f',[0.])

#--------------------------------------------------------------------------------#
# ********************** Filling bkgd data  **********************
print "Storing variables from background."

bkgdSubprocessesListFile = open("bkgd_files")

for processLine in bkgdSubprocessesListFile:
    processLine = processLine.rstrip('\n')
    subProcessName, processName, xsec = processLine.split(" ")
    if subProcessName[0] == "#": continue # problematic input files
    if not processName == process: continue
    # all subProcesses with the same processName will be hstacked with the same
    # color during plotting.
    print
    print "Filling from", subProcessName, "for", processName
    
    # assemble the outName
    outName = outDir+"stopCut_"
    if testMode: outName += "test_"
    else: outName+="all_"
    outName += "Bkgd_"+subProcessName+"_"+channelName+".root"
    outFile = TFile(outName, "recreate")
    
    subProcessListFile = open("/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/v3CutSequence/bkgdProcesses/"+process+"/"+subProcessName)

    # SET UP THE OUTPUT TREE
    tBkgd = TTree("tBkgd", "SUSY stop cut events")
    tBkgd.Branch("muon_count", muon_count, "muon_count/i")
    tBkgd.Branch("muon_pt", muon_pt, "muon_pt[20]/F")
    tBkgd.Branch("muon_eta", muon_eta, "muon_eta[20]/F")
    tBkgd.Branch("muon_phi", muon_phi, "muon_phi[20]/F")
    tBkgd.Branch("muon_relIso", muon_relIso, "muon_relIso[20]/F")
    tBkgd.Branch("muon_charge", muon_charge, "muon_charge[20]/F")
    tBkgd.Branch("muon_mt", muon_mt, "muon_mt[20]/F")
    tBkgd.Branch("electron_count", electron_count, "electron_count/i")
    tBkgd.Branch("electron_pt", electron_pt, "electron_pt[20]/F")
    tBkgd.Branch("electron_eta", electron_eta, "electron_eta[20]/F")
    tBkgd.Branch("electron_phi", electron_phi, "electron_phi[20]/F")
    tBkgd.Branch("electron_relIso", electron_relIso, "electron_relIso[20]/F")
    tBkgd.Branch("electron_charge", electron_charge, "electron_charge[20]/F")
    tBkgd.Branch("electron_mt", electron_mt, "electron_mt[20]/F")
    tBkgd.Branch("found3rdLept", found3rdLept, "found3rdLept/i")
    tBkgd.Branch("lep1_isMu", lep1_isMu, "lep1_isMu/i")
    tBkgd.Branch("lep1_index", lep1_index, "lep1_index/i")
    tBkgd.Branch("lep2_isMu", lep2_isMu, "lep2_isMu/i")
    tBkgd.Branch("lep2_index", lep2_index, "lep2_index/i")
    tBkgd.Branch("njets", njets, "njets/i")
    tBkgd.Branch("jet_pt", jet_pt, "jet_pt[20]/F")
    tBkgd.Branch("jet_eta", jet_eta, "jet_eta[20]/F")
    tBkgd.Branch("jet_phi", jet_phi, "jet_phi[20]/F")
    tBkgd.Branch("jet_ht", jet_ht, "jet_ht/F")
    # tBkgd.Branch("jet_flavour", jet_flavour, "jet_flavour/F")
    tBkgd.Branch("nbtag", nbtag, "nbtag/i")
    tBkgd.Branch("nbtagLoose", nbtagLoose, "nbtagLoose/i")
    tBkgd.Branch("nbtagTight", nbtagTight, "nbtagTight/i")
    tBkgd.Branch("met_pt", met_pt, "met_pt/F")
    tBkgd.Branch("met_phi", met_phi, "met_phi/F")
    tBkgd.Branch("genweight", genweight, "genweight/F")
    
    hGenweights = TH1F("genweights","genweights",1,0,1)
    
    for fileNum, ntuplesLine in enumerate(subProcessListFile):
        if fileNum + 1 > numBkgdFiles: break
        filename = ntuplesLine.rstrip()
        print filename
    
        inFile = TFile.Open(filename, "READ")
        inTree = inFile.Get("AC1B")
        nentries = inTree.GetEntries()
        print("nentries={0:d}".format(nentries))
    
        nMax = nentries
        if testMode: nMax = 5000 
    
        # ***** EVERYTHING BELOW THIS LINE MUST MATCH makeNtupleSigs.py *****
        # ************ BEGIN LOOPING OVER EVENTS **********
        for count, event in enumerate(inTree):
            if count > nMax : break
            if count % 500000 == 0: print("count={0:d}".format(count))
    
            hGenweights.Fill(0,event.genweight)
            # ****** Loose selection of events with valid lep1, lep2, jets ******
            if findingSameFlavor:
                if muPreference:
                    lepIndices = selectMuMu(event, maxOkIso=0.3)
                else: lepIndices = selectElEl(event, maxOkIso=0.3)
                if lepIndices is None: continue
            else:
                lepIndices = selectMuEl(event, maxOkIso=0.3)
                l1Flav = "muon"
                l2Flav = "electron"
                if lepIndices is None:
                    lepIndices = selectElMu(event, maxOkIso=0.3)
                    if lepIndices is None: continue
                    l1Flav = "electron"
                    l2Flav = "muon"
    
            l1Index = lepIndices[0]
            l2Index = lepIndices[1]
    
            jets = findValidJets(event, l1Flav, l1Index, l2Flav, l2Index)
            numGoodJets = len(jets)
    
            numBtag = getNumBtag(event, jets)
            numBtagLoose = getNumBtag(event, jets, 0)
            numBtagTight = getNumBtag(event, jets, 2)
    
            # *********** STORE THE DATA. *************
            # Save all the leptons' and jets' info for this event if it could 
            # possibly contain a good lepton pair.
            assert l1Index > -1
            muon_count[0] = event.muon_count
            for i in range(event.muon_count):
                muon_pt[i] = list(getattr(event, "muon_pt"))[i]
                muon_eta[i] = list(getattr(event, "muon_eta"))[i]
                muon_phi[i] = list(getattr(event, "muon_phi"))[i]
                muon_relIso[i] = list(getattr(event, "muon_relIso"))[i]
                muon_charge[i] = list(getattr(event, "muon_charge"))[i]
                muon_mt[i] = sqrt(2 * muon_pt[i] * event.pfmet_pt * \
                        (1 - cos(muon_phi[i] - event.pfmet_phi)))
    
            assert l2Index > -1
            electron_count[0] = event.electron_count
            for i in range(event.electron_count):
                electron_pt[i] = list(getattr(event, "electron_pt"))[i]
                electron_eta[i] = list(getattr(event, "electron_eta"))[i]
                electron_phi[i] = list(getattr(event, "electron_phi"))[i]
                electron_relIso[i] = list(getattr(event, "electron_relIso"))[i]
                electron_charge[i] = list(getattr(event, "electron_charge"))[i]
                electron_mt[i] = sqrt(2 * electron_pt[i] * event.pfmet_pt * \
                        (1 - cos(electron_phi[i] - event.pfmet_phi)))
    
            lep1_isMu[0] = int(l1Flav == "muon")
            lep1_index[0] = l1Index
            lep2_isMu[0] = int(l2Flav == "muon") 
            lep2_index[0] = l2Index
    
            # veto (3rd lepton) checks:
            found3rdLept[0] = False
            if findingSameFlavor:
                # event should not give valid muel or elmu pair
                if selectMuEl(event) is not None: found3rdLept[0] = True
                if selectElMu(event) is not None: found3rdLept[0] = True
            else:
                # event should not give valid mumu or elel pair
                if selectMuMu(event) is not None: found3rdLept[0] = True
                if selectElEl(event) is not None: found3rdLept[0] = True
    
            njets[0] = numGoodJets
            nbtag[0] = numBtag
            nbtagLoose[0] = numBtagLoose
            nbtagTight[0] = numBtagTight
    
            if numGoodJets > 0:
                jet_ht[0] = 0
                iMaxPtJ = jets[0] 
                for j in range(numGoodJets):
                    jIndex = jets[j]
                    jet_pt[j] = list(event.pfjet_pt)[jIndex]
                    if jet_pt[j] > list(event.pfjet_pt)[iMaxPtJ]:
                        iMaxPtJ = jIndex
                    jet_eta[j] = list(event.pfjet_eta)[jIndex]
                    jet_phi[j] = list(event.pfjet_phi)[jIndex]
                    jet_ht[0] += jet_pt[j]
                    # jet_flavour[j] = list(event.pfjet_flavour)[jIndex]
    
            met_pt[0] = event.pfmet_pt
            met_phi[0] = event.pfmet_phi
            genweight[0] = event.genweight
    
            tBkgd.Fill()
    
    outFile.cd() # cd to outFile to write to it
    tBkgd.Write()
    hGenweights.Write()
    
    outFile.Close()
    print "Finished creating", outName

#--------------------------------------------------------------------------------#

# f = TFile.Open(outName, "READ")
# t = f.Get("tBkgd")
# for event in t:
#     for j in range(event.jet_count):
#         print event.jet_pt[j]
# h = f.Get("bkgd_cutflow")
# h.Sumw2()
# h.Draw()
# raw_input()

print "Done."
