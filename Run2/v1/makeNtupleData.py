#!/usr/bin/env python

# NOTE: NEEDS 3 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, findingSameFlavor, muPreference
#
# Depending on the selected channel, loops over the relevant datasets (SingleMuon 
# for mumu, SingleElectron for elel) and outputs one ntuple for every dataset, 
# located in ../myData/, each containing 1 tree with all events that have survived
# loose dilepton selection cuts.
# The tree in each outputted ntuple, tData, has branches for the same variables 
# as the trees outputted by makeNtupleBkgd.py and makeNtupleSigs.py.
# Uses data_datasets to get names of datasets
# Uses files in dataChannels dir for ntuple lists.

print "Importing modules."
import sys
from ROOT import TFile, TTree, TH1D, TCanvas, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np
from math import sqrt, cos
from array import array
import time
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from jsonChecker import jsonChecker
print "Beginning execution of", sys.argv

assert len(sys.argv) == 4, "need 3 command line args: testMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}"

# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
print "Test mode:", testMode
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
findingSameFlavor = bool(int(sys.argv[2]))
print "Finding same flavor:", findingSameFlavor
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
muPreference = bool(int(sys.argv[3]))
print "Mu preference:", muPreference

if findingSameFlavor:
    if muPreference: 
        l1Flav = "Muon"
        l2Flav = "Muon"
        expectedDatasetChannel = "SingleMuon"
    else: 
        l1Flav = "Electron"
        l2Flav = "Electron"
        expectedDatasetChannel = "SingleElectron"
else: 
    # these 2 lines just matter for creating the outFile name; actual 
    # selection of leading/trailing flavors occurs when looping over events:
    l1Flav = "Muon"
    l2Flav = "Electron"
    assert False, "not supported yet!"
channelName = l1Flav[:2] + l2Flav[:2]

# number of ntuples to loop on for each dataset
# note: must loop over all files and all datasets to have correct xsec
numDataFiles = float("inf") 
if testMode: 
    numDataFiles = 2 

outDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/Run2/"

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
nMuon = array('i',[0])
Muon_pt = np.zeros(20, dtype=np.float32)
Muon_eta = np.zeros(20, dtype=np.float32)
Muon_phi = np.zeros(20, dtype=np.float32)
Muon_relIso = np.zeros(20, dtype=np.float32)
Muon_charge = np.zeros(20, dtype=np.float32)
Muon_mt = np.zeros(20, dtype=np.float32)
nElectron = array('i',[0])
Electron_pt = np.zeros(20, dtype=np.float32)
Electron_eta = np.zeros(20, dtype=np.float32)
Electron_phi = np.zeros(20, dtype=np.float32)
Electron_relIso = np.zeros(20, dtype=np.float32)
Electron_charge = np.zeros(20, dtype=np.float32)
Electron_mt = np.zeros(20, dtype=np.float32)
found3rdLept = array('i',[0])
lep1_isMu = array('i',[0]) # True (1) if muon, False (0) if electron
lep1_index = array('i',[0])
lep1_pt = array('f',[0.])
lep1_eta = array('f',[0.])
lep1_phi = array('f',[0.])
lep1_relIso = array('f',[0.])
lep1_charge = array('f',[0.])
lep1_mt = array('f',[0.])
lep2_isMu = array('i',[0]) # True (1) if muon, False (0) if electron
lep2_index = array('i',[0])
lep2_pt = array('f',[0.])
lep2_eta = array('f',[0.])
lep2_phi = array('f',[0.])
lep2_relIso = array('f',[0.])
lep2_charge = array('f',[0.])
lep2_mt = array('f',[0.])
Jet_pt = np.zeros(20, dtype=np.float32)
Jet_eta = np.zeros(20, dtype=np.float32)
Jet_phi = np.zeros(20, dtype=np.float32)
Jet_ht = array('f',[0.])
# Jet_flavour = array('f',[0])
dR_lep1_jet = array('f',[0.])
dR_lep2_jet = array('f',[0.])
nJet = array('i',[0])
nbtag = array('i',[0])
nbtagLoose = array('i',[0])
nbtagTight = array('i',[0])
met_pt = array('f',[0.])
met_phi = array('f',[0.])
mt_tot = array('f',[0.])
mt_sum = array('f',[0.])
m_eff = array('f',[0.])
run = array('i',[0])
luminosityBlock = array('i',[0])

#--------------------------------------------------------------------------------#
start_time = time.time()
# ********************** Filling events **********************
print "Storing variables from background."

dataDatasetsListFile = open("data_datasets")
jc = jsonChecker()

for datasetLine in dataDatasetsListFile:
    datasetLine = datasetLine.rstrip('\n')
    datasetLine = datasetLine.split()
    datasetName = datasetLine[0]
    datasetChannel = datasetLine[1]
    if datasetName[0] == "#": continue
    if not datasetChannel == expectedDatasetChannel: continue
    print
    print "Filling from", datasetName, "for", datasetChannel 

    # assemble the outName
    outName = outDir+"stopCut_"
    if testMode: outName += "test_"
    else: outName+="all_"
    outName += "Data_"+datasetName+"_"+channelName+".root"
    outFile = TFile(outName, "recreate")
    
    dataNtuplesListFile = open("/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1/dataChannels/"+datasetChannel+"/"+datasetName)

    # SET UP THE OUTPUT TREE
    tData = TTree("tData", "SUSY stop cut events")
    tData.Branch("nMuon", nMuon, "nMuon/i")
    tData.Branch("Muon_pt", Muon_pt, "Muon_pt[20]/F")
    tData.Branch("Muon_eta", Muon_eta, "Muon_eta[20]/F")
    tData.Branch("Muon_phi", Muon_phi, "Muon_phi[20]/F")
    tData.Branch("Muon_relIso", Muon_relIso, "Muon_relIso[20]/F")
    tData.Branch("Muon_charge", Muon_charge, "Muon_charge[20]/F")
    tData.Branch("Muon_mt", Muon_mt, "Muon_mt[20]/F")
    tData.Branch("nElectron", nElectron, "nElectron/i")
    tData.Branch("Electron_pt", Electron_pt, "Electron_pt[20]/F")
    tData.Branch("Electron_eta", Electron_eta, "Electron_eta[20]/F")
    tData.Branch("Electron_phi", Electron_phi, "Electron_phi[20]/F")
    tData.Branch("Electron_relIso", Electron_relIso, "Electron_relIso[20]/F")
    tData.Branch("Electron_charge", Electron_charge, "Electron_charge[20]/F")
    tData.Branch("Electron_mt", Electron_mt, "Electron_mt[20]/F")
    tData.Branch("found3rdLept", found3rdLept, "found3rdLept/i")
    tData.Branch("lep1_isMu", lep1_isMu, "lep1_isMu/i")
    tData.Branch("lep1_index", lep1_index, "lep1_index/i")
    tData.Branch("lep1_pt", lep1_pt, "lep1_pt/F")
    tData.Branch("lep1_eta", lep1_eta, "lep1_eta/F")
    tData.Branch("lep1_phi", lep1_phi, "lep1_phi/F")
    tData.Branch("lep1_relIso", lep1_relIso, "lep1_relIso/F")
    tData.Branch("lep1_charge", lep1_charge, "lep1_charge/F")
    tData.Branch("lep1_mt", lep1_mt, "lep1_mt/F")
    tData.Branch("lep2_isMu", lep2_isMu, "lep2_isMu/i")
    tData.Branch("lep2_index", lep2_index, "lep2_index/i")
    tData.Branch("lep2_pt", lep2_pt, "lep2_pt/F")
    tData.Branch("lep2_eta", lep2_eta, "lep2_eta/F")
    tData.Branch("lep2_phi", lep2_phi, "lep2_phi/F")
    tData.Branch("lep2_relIso", lep2_relIso, "lep2_relIso/F")
    tData.Branch("lep2_charge", lep2_charge, "lep2_charge/F")
    tData.Branch("lep2_mt", lep2_mt, "lep2_mt/F")
    tData.Branch("nJet", nJet, "nJet/i")
    tData.Branch("Jet_pt", Jet_pt, "Jet_pt[20]/F")
    tData.Branch("Jet_eta", Jet_eta, "Jet_eta[20]/F")
    tData.Branch("Jet_phi", Jet_phi, "Jet_phi[20]/F")
    tData.Branch("Jet_ht", Jet_ht, "Jet_ht/F")
    # tData.Branch("Jet_flavour", Jet_flavour, "Jet_flavour/F")
    tData.Branch("dR_lep1_jet", dR_lep1_jet, "dR_lep1_jet/F")
    tData.Branch("dR_lep2_jet", dR_lep2_jet, "dR_lep2_jet/F")
    tData.Branch("nbtag", nbtag, "nbtag/i")
    tData.Branch("nbtagLoose", nbtagLoose, "nbtagLoose/i")
    tData.Branch("nbtagTight", nbtagTight, "nbtagTight/i")
    tData.Branch("met_pt", met_pt, "met_pt/F")
    tData.Branch("met_phi", met_phi, "met_phi/F")
    tData.Branch("mt_tot", mt_tot, "mt_tot/F")
    tData.Branch("mt_sum", mt_sum, "mt_sum/F")
    tData.Branch("m_eff", m_eff, "m_eff/F")
    tData.Branch("run", run, "run/i")
    tData.Branch("luminosityBlock", luminosityBlock, "luminosityBlock/i")
    
    for fileNum, ntuplesLine in enumerate(dataNtuplesListFile):
        if fileNum + 1 > numDataFiles: break
        filename = ntuplesLine.rstrip()
        print filename
    
        inFile = TFile.Open("root://cms-xrd-global.cern.ch//"+filename, "READ")
        inTree = inFile.Get("Events")
        nentries = inTree.GetEntries()
        print "nentries =", nentries
    
        nMax = nentries
        if testMode: nMax = 500 
    
        # ***** EVERYTHING BELOW THIS LINE MUST MATCH THE OTHER MAKENTUPLES *****
        # ************ BEGIN LOOPING OVER EVENTS **********
        for count, event in enumerate(inTree):
            if count > nMax : break
            if count % 500000 == 0: print "count =", count

            if not jc.checkJSON(event.luminosityBlock, event.run): continue

            # ****** Loose selection of events with valid lep1, lep2, jets ******
            if findingSameFlavor:
                if muPreference:
                    lepIndices = selectMuMu(event, maxOkIso=0.3)
                else: lepIndices = selectElEl(event, maxOkIso=0.3)
                if lepIndices is None: continue
            else:
                lepIndices = selectMuEl(event, maxOkIso=0.3)
                l1Flav = "Muon"
                l2Flav = "Electron"
                if lepIndices is None:
                    lepIndices = selectElMu(event, maxOkIso=0.3)
                    if lepIndices is None: continue
                    l1Flav = "Electron"
                    l2Flav = "Muon"
    
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
            nMuon[0] = event.nMuon
            for i in range(event.nMuon):
                Muon_pt[i] = list(event.Muon_pt)[i]
                Muon_eta[i] = list(event.Muon_eta)[i]
                Muon_phi[i] = list(event.Muon_phi)[i]
                Muon_relIso[i] = list(event.Muon_pfRelIso04_all)[i]
                Muon_charge[i] = list(event.Muon_charge)[i]
                Muon_mt[i] = sqrt(2 * Muon_pt[i] * event.MET_pt * \
                        (1 - cos(Muon_phi[i] - event.MET_phi)))
    
            assert l2Index > -1
            nElectron[0] = event.nElectron
            for i in range(event.nElectron):
                Electron_pt[i] = list(event.Electron_pt)[i]
                Electron_eta[i] = list(event.Electron_eta)[i]
                Electron_phi[i] = list(event.Electron_phi)[i]
                Electron_relIso[i] = list(event.Electron_pfRelIso03_all)[i]
                Electron_charge[i] = list(event.Electron_charge)[i]
                Electron_mt[i] = sqrt(2 * Electron_pt[i] * event.MET_pt * \
                        (1 - cos(Electron_phi[i] - event.MET_phi)))
    
            lep1_isMu[0] = int(l1Flav == "Muon")
            lep1_index[0] = l1Index
            lep1_pt[0] = list(getattr(event, l1Flav+"_pt"))[l1Index]
            lep1_eta[0] = list(getattr(event, l1Flav+"_eta"))[l1Index]
            lep1_phi[0] = list(getattr(event, l1Flav+"_phi"))[l1Index]
            if lep1_isMu: lep1_relIso[0] = Muon_relIso[l1Index]
            else: lep1_relIso[0] = Electron_relIso[l1Index]
            lep1_charge[0] = list(getattr(event, l1Flav+"_charge"))[l1Index]
            lep1_mt[0] = sqrt(2 * lep1_pt[0] * event.MET_pt * \
                        (1 - cos(lep1_phi[0] - event.MET_phi)))
            lep2_isMu[0] = int(l2Flav == "Muon") 
            lep2_index[0] = l2Index
            lep2_pt[0] = list(getattr(event, l2Flav+"_pt"))[l2Index]
            lep2_eta[0] = list(getattr(event, l2Flav+"_eta"))[l2Index]
            lep2_phi[0] = list(getattr(event, l2Flav+"_phi"))[l2Index]
            if lep2_isMu: lep2_relIso[0] = Muon_relIso[l2Index]
            else: lep2_relIso[0] = Electron_relIso[l2Index]
            lep2_charge[0] = list(getattr(event, l2Flav+"_charge"))[l2Index]
            lep2_mt[0] = sqrt(2 * lep2_pt[0] * event.MET_pt * \
                        (1 - cos(lep2_phi[0] - event.MET_phi)))
    
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
    
            nJet[0] = numGoodJets
            nbtag[0] = numBtag
            nbtagLoose[0] = numBtagLoose
            nbtagTight[0] = numBtagTight
    
            if numGoodJets > 0:
                Jet_ht[0] = 0
                iMaxPtJ = jets[0] # =index within the Events jets arr
                for j in range(numGoodJets): # =index within the tData jets arr
                    jIndex = jets[j] #  =index within the Events jets arr
                    Jet_pt[j] = list(event.Jet_pt)[jIndex]
                    if Jet_pt[j] > list(event.Jet_pt)[iMaxPtJ]:
                        iMaxPtJ = jIndex
                    Jet_eta[j] = list(event.Jet_eta)[jIndex]
                    Jet_phi[j] = list(event.Jet_phi)[jIndex]
                    Jet_ht[0] += Jet_pt[j]
                    # Jet_flavour[j] = list(event.Jet_flavour)[jIndex]
                dR_lep1_jet[0] = deltaR(event, l1Flav, l1Index, "Jet", iMaxPtJ)
                dR_lep2_jet[0] = deltaR(event, l2Flav, l2Index, "Jet", iMaxPtJ)
    
            met_pt[0] = event.MET_pt
            met_phi[0] = event.MET_phi
            mt_tot[0] = sqrt(lep1_mt[0]**2+ lep2_mt[0]**2)
            mt_sum[0] = lep1_mt[0] + lep2_mt[0]
            m_eff[0] = met_pt[0] + lep1_pt[0]+ lep2_pt[0]
            if nJet[0] > 0: m_eff[0] += Jet_ht[0]
    
            run[0] = event.run
            luminosityBlock[0] = event.luminosityBlock

            tData.Fill()
    
    outFile.cd() # cd to outFile to write to it
    tData.Write()
    
    outFile.Close()
    print "Finished creating", outName

#--------------------------------------------------------------------------------#
print
print int(time.time()-start_time), "secs of processing."

# f = TFile.Open(outName, "READ")
# t = f.Get("tData")
# for event in t:
#     for j in range(event.nJet):
#         print event.Jet_pt[j]
# h.Sumw2()
# h.Draw()
# raw_input()

print "Done."
