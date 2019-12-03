#!/usr/bin/env python

# NOTE: NEEDS 5-6 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# ntuple, subprocess, process (only required for bkgd)
#
# Process options: 
#   (optional arg) data: DoubleMuon, DoubleEG, MuonEG
#   bkgd: TTBar TT+X Diboson W-Jets Drell-Yan Single-Top
#   (optional arg) sig: Stop-Pair
# Subprocess: name of the dataset for data and sig, or name of the subprocess for bkgd
#
# Loops over events in the specified ntuple and outputs another ntuple located in 
# {myDataDir}/{input}/{process}/{subprocess}/
# which contains the Events tree with all events that have survived loose dilepton 
# selection cuts.
#
# Uses Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt (indirectly) 
# for data json checking, and Data_Pileup_2016_271036-284044_80bins.root and 
# MC_Moriond17_PU25ns_V1.root for puWeight calculation

print "Importing modules."
import sys, os
from ROOT import TFile, TTree, TH1D, TCanvas, TImage, TLegend
from ROOT import gSystem, gStyle
import numpy as np
from math import sqrt, cos
from array import array
from collections import OrderedDict
import time
from stopSelection import deltaR,  getBtagIndices, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
from jsonChecker import jsonChecker
print "Beginning execution of", sys.argv

# location of the files (must be in same directory) MC_Moriond17_PU25ns_V1.root, 
# Data_Pileup_2016_271036-284044_80bins.root, 
# Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
# myReferenceDataDir = "/afs/cern.ch/user/a/alkaloge/work/Connie/CMSSW_10_2_9/src/stopSUSY/Run2/"
myReferenceDataDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/Run2"

# location where the output ntuple will be placed
# myDataDir = "/eos/user/a/alkaloge/Connie/"
myDataDir = "/eos/user/c/cmiao/private/myDataSusy/Run2"

assert (len(sys.argv) == 6 or len(sys.argv) == 7), "need 5 or 6 command line args: testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, ntuple, subprocess, process (only required for bkgd)"

# limits the number of events and files to loop over
if sys.argv[1] == "test": testMode = True
elif sys.argv[1] == "all": testMode = False
else: assert False, "invalid test mode, need {test, all}"

# slightly different processing for data, bkgd, and sig
inputType = sys.argv[2]
if inputType == "data":
    isData = True
    isSig = False
elif inputType == "bkgd":
    isData = False
    isSig = False
elif inputType == "sig":
    isData = False
    isSig = True
else: assert False, "invalid type, need {data, bkgd, sig}"

# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
if sys.argv[3] == "mumu":
    findSameFlav = True
    muPref = True
    l1Flav = "Muon"
    l2Flav = "Muon"
    if isData: process = "DoubleMuon"
elif sys.argv[3] == "elel":
    findSameFlav = True
    muPref = False
    l1Flav = "Electron"
    l2Flav = "Electron"
    if isData: process = "DoubleEG"
elif sys.argv[3] == "muel":
    findSameFlav = False
    muPref = False
    l1Flav = "Muon"
    l2Flav = "Electron"
    if isData: process = "MuonEG"
else: assert False, "invalid channel, need {mumu, elel, muel}"
channelName = l1Flav[:2] + l2Flav[:2]

ntupleFileName = sys.argv[4]

subprocess = sys.argv[5]

# handle the parent folder ("process") for the subprocess you want to run
# (ignores sys.argv[6] for data and sig)
if inputType == "bkgd":
    assert len(sys.argv) == 7, "need a process arg for bkgd"
    process = sys.argv[6]
    colorWJets = 38 # dark blue
    colorDY = 46 # red
    colorTTBar = 835 # teal 
    colorSingleTop = 832  
    colorTTX = 831 
    colorDiboson = 806 #orange
    colorQCD = 868 # light blue
    processes = OrderedDict([("W-Jets",colorWJets), ("Drell-Yan",colorDY), 
        ("TTBar",colorTTBar), ("Single-Top",colorSingleTop), ("TT+X",colorTTX), \
        ("Diboson",colorDiboson), ("QCD", colorQCD)])
    assert process in processes, "invalid process %s" % process
elif inputType == "sig": process = "Stop-Pair"

if not isData:
    # dataPileupRoot = TFile.Open(myReferenceDataDir+"/data_pileup_2016.root", "READ")
    dataPileupRoot = TFile.Open(myReferenceDataDir+"/Data_Pileup_2016_271036-284044_80bins.root", "READ")
    dataPileupHist = dataPileupRoot.Get("pileup")
    dataPileupHist.Scale(1/dataPileupHist.Integral())
    # mcPileupRoot = TFile.Open(myReferenceDataDir+"/MC_2016.root", "READ")
    mcPileupRoot = TFile.Open(myReferenceDataDir+"/MC_Moriond17_PU25ns_V1.root", "READ")
    mcPileupHist = mcPileupRoot.Get("pileup")
    mcPileupHist.Scale(1/mcPileupHist.Integral())

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
# Note: the vector branches and the jet branches do not zero out the elements that
# were not filled in the event, so need to check the number of objects to loop over
# before reading these branches (i.e. nMuon, nElectron, nJet).
nMuon = array('i',[0])
Muon_pt = np.zeros(10, dtype=np.float32)
Muon_eta = np.zeros(10, dtype=np.float32)
Muon_phi = np.zeros(10, dtype=np.float32)
Muon_relIso = np.zeros(10, dtype=np.float32)
Muon_charge = np.zeros(10, dtype=np.float32)
Muon_dxy = np.zeros(10, dtype=np.float32)
Muon_dz = np.zeros(10, dtype=np.float32)
Muon_mass = np.zeros(10, dtype=np.float32)
Muon_miniPFRelIso_all = np.zeros(10, dtype=np.float32)
Muon_inTimeMuon = np.zeros(10, dtype=np.int32) # read as bool
Muon_ip3d = np.zeros(10, dtype=np.float32)
Muon_isGlobal = np.zeros(10, dtype=np.int32) # read as bool
Muon_isPFcand = np.zeros(10, dtype=np.int32) # read as bool
Muon_isTracker = np.zeros(10, dtype=np.int32) # read as bool
Muon_jetIdx = np.zeros(10, dtype=np.int32)
Muon_pdgId = np.zeros(10, dtype=np.int32)
Muon_looseId = np.zeros(10, dtype=np.int32) # read as bool
Muon_mediumId = np.zeros(10, dtype=np.int32) # read as bool
Muon_mediumPromptId = np.zeros(10, dtype=np.int32) # read as bool
Muon_mvaId = np.zeros(10, dtype=np.int32)
Muon_tightId = np.zeros(10, dtype=np.int32) # read as bool
if not isData:
    Muon_genPartFlav = np.zeros(10, dtype=np.int32)
    Muon_genPartIdx = np.zeros(10, dtype=np.int32)
Muon_mt = np.zeros(10, dtype=np.float32)
nExtraMuon = array('i',[0])
extraMuIndices = np.zeros(10, dtype=np.int32)

nElectron = array('i',[0])
Electron_pt = np.zeros(10, dtype=np.float32)
Electron_eta = np.zeros(10, dtype=np.float32)
Electron_phi = np.zeros(10, dtype=np.float32)
Electron_relIso = np.zeros(10, dtype=np.float32)
Electron_charge = np.zeros(10, dtype=np.float32)
Electron_dxy = np.zeros(10, dtype=np.float32)
Electron_dz = np.zeros(10, dtype=np.float32)
Electron_mass = np.zeros(10, dtype=np.float32)
Electron_miniPFRelIso_all = np.zeros(10, dtype=np.float32)
Electron_ip3d = np.zeros(10, dtype=np.float32)
Electron_isPFcand = np.zeros(10, dtype=np.int32) # read as bool
Electron_jetIdx = np.zeros(10, dtype=np.int32)
Electron_pdgId = np.zeros(10, dtype=np.int32)
if not isData:
    Electron_genPartFlav = np.zeros(10, dtype=np.int32)
    Electron_genPartIdx = np.zeros(10, dtype=np.int32)
Electron_mt = np.zeros(10, dtype=np.float32)
nExtraElectron = array('i',[0])
extraElIndices = np.zeros(10, dtype=np.int32)

HLT_IsoMu24 = array('i',[0]) # read as bool
HLT_Ele25_eta2p1_WPTight_Gsf = array('i',[0]) # read as bool
HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = array('i',[0]) # read as bool
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = array('i',[0]) # read as bool
HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = array('i',[0]) # read as bool
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = array('i',[0]) # read as bool
nTrigObj = array('i',[0])
TrigObj_filterBits = np.zeros(80, dtype=np.int32)
TrigObj_pt = np.zeros(80, dtype=np.float32)
TrigObj_eta = np.zeros(80, dtype=np.float32)
TrigObj_phi = np.zeros(80, dtype=np.float32)
TrigObj_id = np.zeros(80, dtype=np.int32)
Flag_goodVertices = array('i',[0]) # read as bool
Flag_HBHENoiseFilter = array('i',[0]) # read as bool 
Flag_EcalDeadCellTriggerPrimitiveFilter = array('i',[0]) # read as bool
Flag_BadPFMuonFilter = array('i',[0]) # read as bool

found3rdLept = array('i',[0])
lep1_isMu = array('i',[0]) # read as bool
lep1_index = array('i',[0])
lep2_isMu = array('i',[0]) # read as bool
lep2_index = array('i',[0])

Jet_pt = np.zeros(20, dtype=np.float32)
Jet_eta = np.zeros(20, dtype=np.float32)
Jet_phi = np.zeros(20, dtype=np.float32)
Jet_ht = array('f',[0.])
dR_lep1_jet = array('f',[0.])
dR_lep2_jet = array('f',[0.])
nJet = array('i',[0])
nbtag = array('i',[0])
nbtagLoose = array('i',[0])
nbtagTight = array('i',[0])
btag_indices = np.zeros(20, dtype=np.int32) # indices in the arr of saved good jets
btagLoose_indices = np.zeros(20, dtype=np.int32)
btagTight_indices = np.zeros(20, dtype=np.int32)

MET_pt = array('f',[0.])
MET_phi = array('f',[0.])
MET_significance = array('f',[0.])
MET_sumEt = array('f',[0.])
mt_tot = array('f',[0.])
mt_sum = array('f',[0.])
m_eff = array('f',[0.])

if not isData:
    LHE_HT = array('f',[0.])
    LHE_HTIncoming = array('f',[0.])
    LHE_Nb = array('i',[0])
    LHE_Nuds = array('i',[0])
    LHE_Nglu = array('i',[0])
    LHE_Njets = array('i',[0])
    # LHEPart_pt = array('f',[0.])
    # LHEPart_eta = array('f',[0.])
    # LHEPart_phi = array('f',[0.])
    # LHEPart_mass = array('f',[0.])
    # LHEPart_pdgId = array('i',[0])
    LHEWeight_originalXWGTUP = array('f',[0.])

PV_npvs = array('i',[0])
PV_npvsGood = array('i',[0])
luminosityBlock = array('i',[0])
run = array('i',[0])
if not isData:
    genWeight = array('f',[0.])
    puWeight = array('f',[0.])
    Pileup_nPU = array('i',[0])

#--------------------------------------------------------------------------------#
# ********************** Filling events **********************
print "Storing variables from background."

if isData: jc = jsonChecker(filein=myReferenceDataDir+"/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt")

# SET UP THE OUTPUT TREE
outDir = myDataDir+"/"+inputType+"/"+process+"/"+subprocess+"/"
if not os.path.exists(outDir): os.makedirs(outDir) 
outName = outDir+"stopCut_"
if testMode: outName += "test_"
else: outName+="all_"
outName += ntupleFileName[1+ntupleFileName.rfind("/"):ntupleFileName.rfind(".")]+\
        "_"+channelName+".root"
outFile = TFile(outName, "recreate")
outFile.cd() # cd to outFile to write to it
Events = TTree("Events", "SUSY stop cut events")
Events.Branch("nMuon", nMuon, "nMuon/I")
Events.Branch("Muon_pt", Muon_pt, "Muon_pt[10]/F")
Events.Branch("Muon_eta", Muon_eta, "Muon_eta[10]/F")
Events.Branch("Muon_phi", Muon_phi, "Muon_phi[10]/F")
Events.Branch("Muon_relIso", Muon_relIso, "Muon_relIso[10]/F")
Events.Branch("Muon_charge", Muon_charge, "Muon_charge[10]/F")
Events.Branch("Muon_dxy", Muon_dxy, "Muon_dxy[10]/F")
Events.Branch("Muon_dz", Muon_dz, "Muon_dz[10]/F")
Events.Branch("Muon_mass", Muon_mass, "Muon_mass[10]/F")
Events.Branch("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, "Muon_miniPFRelIso_all[10]/F")
Events.Branch("Muon_inTimeMuon", Muon_inTimeMuon, "Muon_inTimeMuon[10]/I")
Events.Branch("Muon_ip3d", Muon_ip3d, "Muon_ip3d[10]/F")
Events.Branch("Muon_isGlobal", Muon_isGlobal, "Muon_isGlobal[10]/I")
Events.Branch("Muon_isPFcand", Muon_isPFcand, "Muon_isPFcand[10]/I")
Events.Branch("Muon_isTracker", Muon_isTracker, "Muon_isTracker[10]/I")
Events.Branch("Muon_jetIdx", Muon_jetIdx, "Muon_jetIdx[10]/I")
Events.Branch("Muon_pdgId", Muon_pdgId, "Muon_pdgId[10]/I")
Events.Branch("Muon_looseId", Muon_looseId, "Muon_looseId[10]/I")
Events.Branch("Muon_mediumId", Muon_mediumId, "Muon_mediumId[10]/I")
Events.Branch("Muon_mediumPromptId", Muon_mediumPromptId, "Muon_mediumPromptId[10]/I")
Events.Branch("Muon_mvaId", Muon_mvaId, "Muon_mvaId[10]/I")
Events.Branch("Muon_tightId", Muon_tightId, "Muon_tightId[10]/I")
if not isData:
    Events.Branch("Muon_genPartFlav", Muon_genPartFlav, "Muon_genPartFlav[10]/I")
    Events.Branch("Muon_genPartIdx", Muon_genPartIdx, "Muon_genPartIdx[10]/I")
Events.Branch("Muon_mt", Muon_mt, "Muon_mt[10]/F")
Events.Branch("nExtraMuon", nExtraMuon, "nExtraMuon/I")
Events.Branch("extraMuIndices", extraMuIndices, "extraMuIndices[10]/I")

Events.Branch("nElectron", nElectron, "nElectron/I")
Events.Branch("Electron_pt", Electron_pt, "Electron_pt[10]/F")
Events.Branch("Electron_eta", Electron_eta, "Electron_eta[10]/F")
Events.Branch("Electron_phi", Electron_phi, "Electron_phi[10]/F")
Events.Branch("Electron_relIso", Electron_relIso, "Electron_relIso[10]/F")
Events.Branch("Electron_charge", Electron_charge, "Electron_charge[10]/F")
Events.Branch("Electron_dxy", Electron_dxy, "Electron_dxy[10]/F")
Events.Branch("Electron_dz", Electron_dz, "Electron_dz[10]/F")
Events.Branch("Electron_mass", Electron_mass, "Electron_mass[10]/F")
Events.Branch("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, "Electron_miniPFRelIso_all[10]/F")
Events.Branch("Electron_ip3d", Electron_ip3d, "Electron_ip3d[10]/F")
Events.Branch("Electron_isPFcand", Electron_isPFcand, "Electron_isPFcand[10]/I")
Events.Branch("Electron_jetIdx", Electron_jetIdx, "Electron_jetIdx[10]/I")
Events.Branch("Electron_pdgId", Electron_pdgId, "Electron_pdgId[10]/I")
if not isData:
    Events.Branch("Electron_genPartFlav", Electron_genPartFlav, "Electron_genPartFlav[10]/I")
    Events.Branch("Electron_genPartIdx", Electron_genPartIdx, "Electron_genPartIdx[10]/I")
Events.Branch("Electron_mt", Electron_mt, "Electron_mt[10]/F")
Events.Branch("nExtraElectron", nExtraElectron, "nExtraElectron/I")
Events.Branch("extraElIndices", extraElIndices, "extraElIndices[10]/I")

Events.Branch("HLT_IsoMu24", HLT_IsoMu24, "HLT_IsoMu24/I")
Events.Branch("HLT_Ele25_eta2p1_WPTight_Gsf", HLT_Ele25_eta2p1_WPTight_Gsf, "HLT_Ele25_eta2p1_WPTight_Gsf/I")
Events.Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL/I")
Events.Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/I")
Events.Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/I")
Events.Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/I")
Events.Branch("nTrigObj", nTrigObj, "nTrigObj/I")
Events.Branch("TrigObj_filterBits", TrigObj_filterBits, "TrigObj_filterBits[80]/I")
Events.Branch("TrigObj_pt", TrigObj_pt, "TrigObj_pt[80]/F")
Events.Branch("TrigObj_eta", TrigObj_eta, "TrigObj_eta[80]/F")
Events.Branch("TrigObj_phi", TrigObj_phi, "TrigObj_phi[80]/F")
Events.Branch("TrigObj_id", TrigObj_id, "TrigObj_id[80]/I")
Events.Branch("Flag_goodVertices", Flag_goodVertices, "Flag_goodVertices/I")
Events.Branch("Flag_HBHENoiseFilter", Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/I")
Events.Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/I")
Events.Branch("Flag_BadPFMuonFilter", Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/I")
Events.Branch("found3rdLept", found3rdLept, "found3rdLept/I")
Events.Branch("lep1_isMu", lep1_isMu, "lep1_isMu/I")
Events.Branch("lep1_index", lep1_index, "lep1_index/I")
Events.Branch("lep2_isMu", lep2_isMu, "lep2_isMu/I")
Events.Branch("lep2_index", lep2_index, "lep2_index/I")
Events.Branch("nJet", nJet, "nJet/I")
Events.Branch("Jet_pt", Jet_pt, "Jet_pt[20]/F")
Events.Branch("Jet_eta", Jet_eta, "Jet_eta[20]/F")
Events.Branch("Jet_phi", Jet_phi, "Jet_phi[20]/F")
Events.Branch("Jet_ht", Jet_ht, "Jet_ht/F")
Events.Branch("dR_lep1_jet", dR_lep1_jet, "dR_lep1_jet/F")
Events.Branch("dR_lep2_jet", dR_lep2_jet, "dR_lep2_jet/F")
Events.Branch("nbtag", nbtag, "nbtag/I")
Events.Branch("nbtagLoose", nbtagLoose, "nbtagLoose/I")
Events.Branch("nbtagTight", nbtagTight, "nbtagTight/I")
Events.Branch("btag_indices", btag_indices, "btag_indices[20]/I")
Events.Branch("btagLoose_indices", btagLoose_indices, "btagLoose_indices[20]/I")
Events.Branch("btagTight_indices", btagTight_indices, "btagTight_indices[20]/I")
Events.Branch("MET_pt", MET_pt, "MET_pt/F")
Events.Branch("MET_phi", MET_phi, "MET_phi/F")
Events.Branch("MET_significance", MET_significance, "MET_significance/F")
Events.Branch("MET_sumEt", MET_sumEt, "MET_sumEt/F")
Events.Branch("mt_tot", mt_tot, "mt_tot/F")
Events.Branch("mt_sum", mt_sum, "mt_sum/F")
Events.Branch("m_eff", m_eff, "m_eff/F")
if not isData:
    Events.Branch("LHE_HT", LHE_HT, "LHE_HT/F")
    Events.Branch("LHE_HTIncoming", LHE_HTIncoming, "LHE_HTIncoming/F")
    Events.Branch("LHE_Nb", LHE_Nb, "LHE_Nb/I")
    Events.Branch("LHE_Nuds", LHE_Nuds, "LHE_Nuds/I")
    Events.Branch("LHE_Nglu", LHE_Nglu, "LHE_Nglu/I")
    Events.Branch("LHE_Njets", LHE_Njets, "LHE_Njets/I")
    # Events.Branch("LHEPart_pt", LHEPart_pt, "LHEPart_pt/F")
    # Events.Branch("LHEPart_eta", LHEPart_eta, "LHEPart_eta/F")
    # Events.Branch("LHEPart_phi", LHEPart_phi, "LHEPart_phi/F")
    # Events.Branch("LHEPart_mass", LHEPart_mass, "LHEPart_mass/F")
    # Events.Branch("LHEPart_pdgId", LHEPart_pdgId, "LHEPart_pdgId/I")
    Events.Branch("LHEWeight_originalXWGTUP", LHEWeight_originalXWGTUP, "LHEWeight_originalXWGTUP/F")
Events.Branch("PV_npvs", PV_npvs, "PV_npvs/I")
Events.Branch("PV_npvsGood", PV_npvsGood, "PV_npvsGood/I")
Events.Branch("luminosityBlock", luminosityBlock, "luminosityBlock/I")
Events.Branch("run", run, "run/I")
if not isData:
    Events.Branch("genWeight", genWeight, "genWeight/F")
    Events.Branch("puWeight", puWeight, "puWeight/F")
    Events.Branch("Pileup_nPU", Pileup_nPU, "Pileup_nPU/I")

if not isData:
    if subprocess == "WJetsToLNu":
        hWxGenweightsArr = []
        for i in range(5):
            hWxGenweightsArr.append(TH1D("W"+str(i)+"genWeights",\
                    "W"+str(i)+"genWeights",1,-0.5,0.5))
    elif subprocess == "DYJetsToLL_M-50":
        hDYxGenweightsArr = []
        for i in range(5):
            hDYxGenweightsArr.append(TH1D("DY"+str(i)+"genWeights",\
                    "DY"+str(i)+"genWeights",1,-0.5,0.5))
    hGenweights = TH1D("genWeights","genWeights",1,-0.5,0.5)

try:
    if testMode: inFile = TFile.Open(ntupleFileName, "READ") # use test_WJets.root
    else: 
        inFile = TFile.Open("root://cms-xrd-global.cern.ch//"+ntupleFileName, "READ")
except: exit()


try: inTree = inFile.Get("Events")
except: exit()
nentries = inTree.GetEntries()
print "nentries =", nentries

nMax = nentries
if testMode: nMax = 10000 

# ************ BEGIN LOOPING OVER EVENTS **********
if not isData:
    # skip evts with < 0 genWeight if it's a madgraph file
    isMadgraph = False
    if "madgraph" in ntupleFileName: isMadgraph = True
for count, event in enumerate(inTree):
    if count > nMax : break
    if count == 0: start_time = time.time()
    if count % 5000 == 0: print "count =", count

    if isData:
        if not jc.checkJSON(event.luminosityBlock, event.run): continue
    else:
        if isMadgraph:
            if event.genWeight < 0: continue
        npartons = ord(event.LHE_Njets)
        if subprocess == "WJetsToLNu" and npartons <= 4:
            hWxGenweightsArr[npartons].Fill(0, event.genWeight)
        if subprocess == "DYJetsToLL_M-50" and npartons <= 4:
            hDYxGenweightsArr[npartons].Fill(0, event.genWeight)
        hGenweights.Fill(0, event.genWeight)

        # if hGenweights.GetSumOfWeights() >= hGenweights.GetEntries(): 
        #     print "WARNING: evt #", count, "genwt", event.genWeight, \
        #             "sumw", hGenweights.GetSumOfWeights(), \
        #             "nentries", hGenweights.GetEntries() 

    # ****** Loose selection for valid lep1, lep2, jets (with trigger) ******
    if findSameFlav:
        if muPref:
            lepIndices = selectMuMu(event, isData, maxOkIso=0.3)
        else: lepIndices = selectElEl(event, isData, maxOkIso=0.3)
        if lepIndices is None: continue
    else:
        lepIndices = selectMuEl(event, isData, maxOkIso=0.3)
        l1Flav = "Muon"
        l2Flav = "Electron"
        if lepIndices is None:
            lepIndices = selectElMu(event, isData, maxOkIso=0.3)
            if lepIndices is None: continue
            l1Flav = "Electron"
            l2Flav = "Muon"

    l1Index = lepIndices[0]
    l2Index = lepIndices[1]
    evt_extraMuIndices = lepIndices[2]
    evt_extraElIndices = lepIndices[3]
    print l1Index, l2Index, evt_extraMuIndices, evt_extraElIndices

    jets = findValidJets(event, l1Flav, l1Index, l2Flav, l2Index)
    numGoodJets = len(jets)

    evt_btag_indices = getBtagIndices(event, jets)
    evt_btagLoose_indices = getBtagIndices(event, jets, 0)
    evt_btagTight_indices = getBtagIndices(event, jets, 2)

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
        Muon_dxy[i] = list(event.Muon_dxy)[i]
        Muon_dz[i] = list(event.Muon_dz)[i]
        Muon_mass[i] = list(event.Muon_mass)[i]
        Muon_miniPFRelIso_all[i] = list(event.Muon_miniPFRelIso_all)[i]
        Muon_inTimeMuon[i] = list(event.Muon_inTimeMuon)[i]
        Muon_ip3d[i] = list(event.Muon_ip3d)[i]
        Muon_isGlobal[i] = list(event.Muon_isGlobal)[i]
        Muon_isPFcand[i] = list(event.Muon_isPFcand)[i]
        Muon_isTracker[i] = list(event.Muon_isTracker)[i]
        Muon_jetIdx[i] = list(event.Muon_jetIdx)[i]
        Muon_pdgId[i] = list(event.Muon_pdgId)[i]
        Muon_looseId[i] = list(event.Muon_looseId)[i]
        Muon_mediumId[i] = list(event.Muon_mediumId)[i]
        Muon_mediumPromptId[i] = list(event.Muon_mediumPromptId)[i]
        Muon_mvaId[i] = ord(list(event.Muon_mvaId)[i]) # UChar_t conversion
        Muon_tightId[i] = list(event.Muon_tightId)[i]
        if not isData:
            Muon_genPartFlav[i] = ord(list(event.Muon_genPartFlav)[i])
            Muon_genPartIdx[i] = list(event.Muon_genPartIdx)[i]
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
        Electron_dxy[i] = list(event.Electron_dxy)[i]
        Electron_dz[i] = list(event.Electron_dz)[i]
        Electron_mass[i] = list(event.Electron_mass)[i]
        Electron_miniPFRelIso_all[i] = list(event.Electron_miniPFRelIso_all)[i]
        Electron_ip3d[i] = list(event.Electron_ip3d)[i]
        Electron_isPFcand[i] = list(event.Electron_isPFcand)[i]
        Electron_jetIdx[i] = list(event.Electron_jetIdx)[i]
        Electron_pdgId[i] = list(event.Electron_pdgId)[i]
        if not isData:
            Electron_genPartFlav[i] = ord(list(event.Electron_genPartFlav)[i])
            Electron_genPartIdx[i] = list(event.Electron_genPartIdx)[i]
        Electron_mt[i] = sqrt(2 * Electron_pt[i] * event.MET_pt * \
                (1 - cos(Electron_phi[i] - event.MET_phi)))

    Flag_goodVertices[0] = event.Flag_goodVertices
    Flag_HBHENoiseFilter[0] = event.Flag_HBHENoiseFilter
    Flag_EcalDeadCellTriggerPrimitiveFilter[0] = event.Flag_EcalDeadCellTriggerPrimitiveFilter
    Flag_BadPFMuonFilter[0] = event.Flag_BadPFMuonFilter

    HLT_IsoMu24[0] = event.HLT_IsoMu24
    HLT_Ele25_eta2p1_WPTight_Gsf[0] = event.HLT_Ele25_eta2p1_WPTight_Gsf
    if isData and event.run >= 278820: # Run G or higher
        HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ[0] = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ[0] = event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
    else:
        HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL[0] = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
        HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL[0] = event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL

    nTrigObj[0] = event.nTrigObj
    for i in range(event.nTrigObj):
        if i>80: 
            sys.stderr.write("WARNING: more than 80 trig objs, skipping rest!\n")
            break
        TrigObj_filterBits[i] = list(event.TrigObj_filterBits)[i]
        TrigObj_pt[i] = list(event.TrigObj_pt)[i]
        TrigObj_eta[i] = list(event.TrigObj_eta)[i]
        TrigObj_phi[i] = list(event.TrigObj_phi)[i]
        TrigObj_id[i] = list(event.TrigObj_id)[i]

    lep1_isMu[0] = int(l1Flav == "Muon")
    lep1_index[0] = l1Index
    lep2_isMu[0] = int(l2Flav == "Muon") 
    lep2_index[0] = l2Index
    nExtraMuon[0] = len(evt_extraMuIndices)
    nExtraElectron[0] = len(evt_extraElIndices)
    for i, iExtraMu in enumerate(evt_extraMuIndices):
        extraMuIndices[i] = iExtraMu
    for i, iExtraEl in enumerate(evt_extraElIndices):
        extraElIndices[i] = iExtraEl

    # veto (3rd lepton) checks:
    found3rdLept[0] = False
    # if findSameFlav:
    #     # when looking for mumu/elel, event should not give valid muel or elmu pair
    #     if selectMuEl(event, isData) is not None: found3rdLept[0] = True
    #     if selectElMu(event, isData) is not None: found3rdLept[0] = True
    # else:
    #     # when looking for muel/elmu, event should not give valid mumu or elel pair
    #     if selectMuMu(event, isData) is not None: found3rdLept[0] = True
    #     if selectElEl(event, isData) is not None: found3rdLept[0] = True
    if nExtraMuon > 0 or nExtraElectron > 0: found3rdLept[0] = True
    
    nJet[0] = numGoodJets
    nbtag[0] = len(evt_btag_indices)
    nbtagLoose[0] = len(evt_btagLoose_indices)
    nbtagTight[0] = len(evt_btagTight_indices)

    # j = an index within the btag_indices arr; iGoodJets = an index within
    # the good jets array
    for j, iGoodJets in enumerate(evt_btag_indices):
        btag_indices[j] = iGoodJets
    for j, iGoodJets in enumerate(evt_btagLoose_indices):
        btagLoose_indices[j] = iGoodJets
    for j, iGoodJets in enumerate(evt_btagTight_indices):
        btagTight_indices[j] = iGoodJets

    if numGoodJets > 0:
        Jet_ht[0] = 0
        iMaxPtJ = jets[0] # =index within the Events jets arr
        for j in range(numGoodJets): # =index within the Events jets arr
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

    MET_pt[0] = event.MET_pt
    MET_phi[0] = event.MET_phi
    MET_significance[0] = event.MET_significance
    MET_sumEt[0] = event.MET_sumEt

    if bool(lep1_isMu[0]): 
        lep1_mt = Muon_mt[l1Index]
        lep1_pt = Muon_pt[l1Index]
    else:
        lep1_mt = Electron_mt[l1Index]
        lep1_pt = Electron_pt[l1Index]
    if bool(lep2_isMu[0]): 
        lep2_mt = Muon_mt[l2Index]
        lep2_pt = Muon_pt[l2Index]
    else:
        lep2_mt = Electron_mt[l2Index]
        lep2_pt = Electron_pt[l2Index]

    mt_tot[0] = sqrt(lep1_mt**2 + lep2_mt**2)
    mt_sum[0] = lep1_mt + lep2_mt
    m_eff[0] = MET_pt[0] + lep1_pt + lep2_pt
    if nJet[0] > 0: m_eff[0] += Jet_ht[0]

    if not isData:
        LHE_HT[0] = event.LHE_HT
        LHE_HTIncoming[0] = event.LHE_HTIncoming
        LHE_Njets[0] = ord(event.LHE_Njets)
        LHE_Nb[0] = ord(event.LHE_Nb)
        LHE_Nuds[0] = ord(event.LHE_Nuds)
        LHE_Nglu[0] = ord(event.LHE_Nglu)
        # LHEPart_pt[0] = event.LHEPart_pt
        # LHEPart_eta[0] = event.LHEPart_eta
        # LHEPart_phi[0] = event.LHEPart_phi
        # LHEPart_mass[0] = event.LHEPart_mass
        # LHEPart_pdgId[0] = event.LHEPart_pdgId
        LHEWeight_originalXWGTUP[0] = event.LHEWeight_originalXWGTUP

    PV_npvs[0] = event.PV_npvs
    PV_npvsGood[0] = event.PV_npvsGood
    luminosityBlock[0] = event.luminosityBlock
    run[0] = event.run
    if not isData:
        genWeight[0] = event.genWeight
        # calculate pileup weight
        evtPU = event.Pileup_nPU
        dataPUBin = dataPileupHist.FindBin(evtPU)
        mcPUBin = mcPileupHist.FindBin(evtPU)
        if mcPileupHist.GetBinContent(mcPUBin) == 0:
            puWeight[0] = 1
        else:
            puWeight[0] = dataPileupHist.GetBinContent(dataPUBin) / \
                    mcPileupHist.GetBinContent(mcPUBin)
        Pileup_nPU[0] = evtPU

    Events.Fill()

outFile.cd()
print "Saving", Events.GetEntries(), "events."
Events.Write()
if not isData: 
    if subprocess == "WJetsToLNu":
        for i in range(len(hWxGenweightsArr)):
            hWxGenweightsArr[i].Write()
    elif subprocess == "DYJetsToLL_M-50":
        for i in range(len(hDYxGenweightsArr)):
            hDYxGenweightsArr[i].Write()
    hGenweights.Write()
outFile.Close()
print "Finished creating", outName

#--------------------------------------------------------------------------------#
print
print int(time.time()-start_time), "secs of processing."

# f = TFile.Open(outName, "READ")
# t = f.Get("Events")
# for event in t:
#     for j in range(event.nJet):
#         print event.Jet_pt[j]
# h.Sumw2()
# h.Draw()
# raw_input()

print "Done."
