#!/usr/bin/env python

# NOTE: NEEDS 4 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# process
# Process options: 
#   data: DoubleMuon, DoubleEG, MuonEG
#   bkgd: TTBar TT+X Diboson W-Jets Drell-Yan Single-Top
#   sig: Stop-Pair
#
# Depending on the selected channel, loops over the relevant datasets (SingleMuon 
# for mumu, SingleElectron for elel) and outputs one ntuple for every dataset, 
# located in ../myData/, each containing the Events tree with all events that have
# survived loose dilepton selection cuts.
#
# Uses {data/bkgd/sig}_fileRedirector to get names of datasets and for MC, xsecs too
# Uses files in {data/bkgd/sig}NtupleLists/{process}/ dir for ntuple lists.
# Uses Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt for data 
# json checking

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

assert len(sys.argv) == 5, "need 4 command line args: testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, process"

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
redirectorFileAdr = sys.argv[2]+"_fileRedirector"

# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
if sys.argv[3] == "mumu":
    findingSameFlavor = True
    muPreference = True
    l1Flav = "Muon"
    l2Flav = "Muon"
    if isData: expectedSubfolder = "DoubleMuon"
elif sys.argv[3] == "elel":
    findingSameFlavor = True
    muPreference = False
    l1Flav = "Electron"
    l2Flav = "Electron"
    if isData: expectedSubfolder = "DoubleEG"
elif sys.argv[3] == "muel":
    findingSameFlavor = False
    muPreference = False
    l1Flav = "Muon"
    l2Flav = "Electron"
    if isData: expectedSubfolder = "MuonEG"
else: assert False, "invalid channel, need {mumu, elel, muel}"
channelName = l1Flav[:2] + l2Flav[:2]

# handle the parent folder for the dataset you want to run
if isData:
    assert sys.argv[4] == expectedSubfolder, \
            "mismatch b/w selected process and channel; for "+sys.argv[3]+\
            ", process should be "+expectedSubfolder
elif isSig:
    expectedSubfolder = "Stop-Pair"
    assert sys.argv[4] == expectedSubfolder, \
            "signal input type expects 'Stop-Pair' as the process"
else: # bkgd
    expectedSubfolder = sys.argv[4]
    processes = {"W-Jets":38, "Drell-Yan":46, "TTBar":30, "Diboson":41, \
            "Single-Top":40, "TT+X":7}
    assert expectedSubfolder in processes, "invalid process %s" % expectedSubfolder

# number of ntuples to loop on for each dataset
# note: must loop over all files and all datasets to have correct xsec
numNtuples = float("inf") 
if testMode: 
    numNtuples = 2 

outDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/Run2/"

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
nMuon = array('i',[0])
Muon_pt = np.zeros(20, dtype=np.float32)
Muon_eta = np.zeros(20, dtype=np.float32)
Muon_phi = np.zeros(20, dtype=np.float32)
Muon_relIso = np.zeros(20, dtype=np.float32)
Muon_charge = np.zeros(20, dtype=np.float32)
Muon_dxy = np.zeros(20, dtype=np.float32)
Muon_dz = np.zeros(20, dtype=np.float32)
Muon_mass = np.zeros(20, dtype=np.float32)
Muon_miniPFRelIso_all = np.zeros(20, dtype=np.float32)
Muon_inTimeMuon = np.zeros(20, dtype=np.int32) # read as bool
Muon_ip3d = np.zeros(20, dtype=np.float32)
Muon_isGlobal = np.zeros(20, dtype=np.int32) # read as bool
Muon_isPFcand = np.zeros(20, dtype=np.int32) # read as bool
Muon_isTracker = np.zeros(20, dtype=np.int32) # read as bool
Muon_jetIdx = np.zeros(20, dtype=np.int32)
Muon_pdgId = np.zeros(20, dtype=np.int32)
Muon_looseId = np.zeros(20, dtype=np.int32) # read as bool
Muon_mediumId = np.zeros(20, dtype=np.int32) # read as bool
Muon_mediumPromptId = np.zeros(20, dtype=np.int32) # read as bool
Muon_mvaId = np.zeros(20, dtype=np.int32)
Muon_tightId = np.zeros(20, dtype=np.int32) # read as bool
Muon_genPartFlav = np.zeros(20, dtype=np.int32)
Muon_genPartIdx = np.zeros(20, dtype=np.int32)
Muon_mt = np.zeros(20, dtype=np.float32)

nElectron = array('i',[0])
Electron_pt = np.zeros(20, dtype=np.float32)
Electron_eta = np.zeros(20, dtype=np.float32)
Electron_phi = np.zeros(20, dtype=np.float32)
Electron_relIso = np.zeros(20, dtype=np.float32)
Electron_charge = np.zeros(20, dtype=np.float32)
Electron_dxy = np.zeros(20, dtype=np.float32)
Electron_dz = np.zeros(20, dtype=np.float32)
Electron_mass = np.zeros(20, dtype=np.float32)
Electron_miniPFRelIso_all = np.zeros(20, dtype=np.float32)
# Electron_inTimeMuon = np.zeros(20, dtype=np.int32) # read as bool
Electron_ip3d = np.zeros(20, dtype=np.float32)
# Electron_isGlobal = np.zeros(20, dtype=np.int32) # read as bool
Electron_isPFcand = np.zeros(20, dtype=np.int32) # read as bool
# Electron_isTracker = np.zeros(20, dtype=np.int32) # read as bool
Electron_jetIdx = np.zeros(20, dtype=np.int32)
Electron_pdgId = np.zeros(20, dtype=np.int32)
# Electron_looseId = np.zeros(20, dtype=np.int32) # read as bool
# Electron_mediumId = np.zeros(20, dtype=np.int32) # read as bool
# Electron_mediumPromptId = np.zeros(20, dtype=np.int32) # read as bool
# Electron_mvaId = np.zeros(20, dtype=np.int32)
# Electron_tightId = np.zeros(20, dtype=np.int32) # read as bool
Electron_genPartFlav = np.zeros(20, dtype=np.int32)
Electron_genPartIdx = np.zeros(20, dtype=np.int32)
Electron_mt = np.zeros(20, dtype=np.float32)

found3rdLept = array('i',[0])
lep1_isMu = array('i',[0]) # read as bool
lep1_index = array('i',[0])
lep2_isMu = array('i',[0]) # read as bool
lep2_index = array('i',[0])

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

MET_pt = array('f',[0.])
MET_phi = array('f',[0.])
MET_significance = array('f',[0.])
MET_sumEt = array('f',[0.])
mt_tot = array('f',[0.])
mt_sum = array('f',[0.])
m_eff = array('f',[0.])

LHE_HT = array('f',[0.])
LHE_HTIncoming = array('f',[0.])
LHE_Njets = array('i',[0])
LHE_Nb = array('i',[0])
LHE_Nuds = array('i',[0])
LHE_Nglu = array('i',[0])
LHE_NpNLO = array('i',[0])
LHE_NpLO = array('i',[0])
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
genWeight = array('f',[0.])

#--------------------------------------------------------------------------------#
start_time = time.time()
# ********************** Filling events **********************
print "Storing variables from background."

redirectorFile = open(redirectorFileAdr)

if isData: jc = jsonChecker()

for line in redirectorFile:
    line = line.rstrip('\n').split()
    if line[0][0] == "#": continue
    datasetNtuplesList = line[0]
    subfolder = line[1]
    if not subfolder == expectedSubfolder: continue
    if not isData:
        # skip evts with < 0 genWeight if it's a madgraph file
        isMadgraph = False
        if "madgraph" in  datasetNtuplesList: isMadgraph = True
    print
    print "Filling from", datasetNtuplesList, "from", subfolder

    # assemble the outName
    outName = outDir+"stopCut_"
    if testMode: outName += "test_"
    else: outName+="all_"
    outName += inputType+"_"+datasetNtuplesList+"_"+channelName+".root"
    outFile = TFile(outName, "recreate")
    
    ntuplesListFile = open("/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1/"+inputType+"NtupleLists/"+subfolder+"/"+datasetNtuplesList)

    # SET UP THE OUTPUT TREE
    Events = TTree("Events", "SUSY stop cut events")
    Events.Branch("nMuon", nMuon, "nMuon/I")
    Events.Branch("Muon_pt", Muon_pt, "Muon_pt[20]/F")
    Events.Branch("Muon_eta", Muon_eta, "Muon_eta[20]/F")
    Events.Branch("Muon_phi", Muon_phi, "Muon_phi[20]/F")
    Events.Branch("Muon_relIso", Muon_relIso, "Muon_relIso[20]/F")
    Events.Branch("Muon_charge", Muon_charge, "Muon_charge[20]/F")
    Events.Branch("Muon_dxy", Muon_dxy, "Muon_dxy[20]/F")
    Events.Branch("Muon_dz", Muon_dz, "Muon_dz[20]/F")
    Events.Branch("Muon_mass", Muon_mass, "Muon_mass[20]/F")
    Events.Branch("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, "Muon_miniPFRelIso_all[20]/F")
    Events.Branch("Muon_inTimeMuon", Muon_inTimeMuon, "Muon_inTimeMuon[20]/I")
    Events.Branch("Muon_ip3d", Muon_ip3d, "Muon_ip3d[20]/F")
    Events.Branch("Muon_isGlobal", Muon_isGlobal, "Muon_isGlobal[20]/I")
    Events.Branch("Muon_isPFcand", Muon_isPFcand, "Muon_isPFcand[20]/I")
    Events.Branch("Muon_isTracker", Muon_isTracker, "Muon_isTracker[20]/I")
    Events.Branch("Muon_jetIdx", Muon_jetIdx, "Muon_jetIdx[20]/I")
    Events.Branch("Muon_pdgId", Muon_pdgId, "Muon_pdgId[20]/I")
    Events.Branch("Muon_looseId", Muon_looseId, "Muon_looseId[20]/I")
    Events.Branch("Muon_mediumId", Muon_mediumId, "Muon_mediumId[20]/I")
    Events.Branch("Muon_mediumPromptId", Muon_mediumPromptId, "Muon_mediumPromptId[20]/I")
    Events.Branch("Muon_mvaId", Muon_mvaId, "Muon_mvaId[20]/I")
    Events.Branch("Muon_tightId", Muon_tightId, "Muon_tightId[20]/I")
    Events.Branch("Muon_genPartFlav", Muon_genPartFlav, "Muon_genPartFlav[20]/I")
    Events.Branch("Muon_genPartFlav", Muon_genPartFlav, "Muon_genPartFlav[20]/I")
    Events.Branch("Muon_genPartIdx", Muon_genPartIdx, "Muon_genPartIdx[20]/I")
    Events.Branch("Muon_mt", Muon_mt, "Muon_mt[20]/F")
    Events.Branch("nElectron", nElectron, "nElectron/I")
    Events.Branch("Electron_pt", Electron_pt, "Electron_pt[20]/F")
    Events.Branch("Electron_eta", Electron_eta, "Electron_eta[20]/F")
    Events.Branch("Electron_phi", Electron_phi, "Electron_phi[20]/F")
    Events.Branch("Electron_relIso", Electron_relIso, "Electron_relIso[20]/F")
    Events.Branch("Electron_charge", Electron_charge, "Electron_charge[20]/F")
    Events.Branch("Electron_dxy", Electron_dxy, "Electron_dxy[20]/F")
    Events.Branch("Electron_dz", Electron_dz, "Electron_dz[20]/F")
    Events.Branch("Electron_mass", Electron_mass, "Electron_mass[20]/F")
    Events.Branch("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, "Electron_miniPFRelIso_all[20]/F")
    Events.Branch("Electron_ip3d", Electron_ip3d, "Electron_ip3d[20]/F")
    Events.Branch("Electron_isPFcand", Electron_isPFcand, "Electron_isPFcand[20]/I")
    Events.Branch("Electron_jetIdx", Electron_jetIdx, "Electron_jetIdx[20]/I")
    Events.Branch("Electron_pdgId", Electron_pdgId, "Electron_pdgId[20]/I")
    Events.Branch("Electron_genPartFlav", Electron_genPartFlav, "Electron_genPartFlav[20]/I")
    Events.Branch("Electron_genPartFlav", Electron_genPartFlav, "Electron_genPartFlav[20]/I")
    Events.Branch("Electron_genPartIdx", Electron_genPartIdx, "Electron_genPartIdx[20]/I")
    Events.Branch("Electron_mt", Electron_mt, "Electron_mt[20]/F")
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
    # Events.Branch("Jet_flavour", Jet_flavour, "Jet_flavour/F")
    Events.Branch("dR_lep1_jet", dR_lep1_jet, "dR_lep1_jet/F")
    Events.Branch("dR_lep2_jet", dR_lep2_jet, "dR_lep2_jet/F")
    Events.Branch("nbtag", nbtag, "nbtag/I")
    Events.Branch("nbtagLoose", nbtagLoose, "nbtagLoose/I")
    Events.Branch("nbtagTight", nbtagTight, "nbtagTight/I")
    Events.Branch("MET_pt", MET_pt, "MET_pt/F")
    Events.Branch("MET_phi", MET_phi, "MET_phi/F")
    Events.Branch("MET_significance", MET_significance, "MET_significance/F")
    Events.Branch("MET_sumEt", MET_sumEt, "MET_sumEt/F")
    Events.Branch("mt_tot", mt_tot, "mt_tot/F")
    Events.Branch("mt_sum", mt_sum, "mt_sum/F")
    Events.Branch("m_eff", m_eff, "m_eff/F")
    Events.Branch("LHE_HT", LHE_HT, "LHE_HT/F")
    Events.Branch("LHE_HTIncoming", LHE_HTIncoming, "LHE_HTIncoming/F")
    Events.Branch("LHE_Njets", LHE_Njets, "LHE_Njets/I")
    Events.Branch("LHE_Nb", LHE_Nb, "LHE_Nb/I")
    Events.Branch("LHE_Nuds", LHE_Nuds, "LHE_Nuds/I")
    Events.Branch("LHE_Nglu", LHE_Nglu, "LHE_Nglu/I")
    Events.Branch("LHE_NpNLO", LHE_NpNLO, "LHE_NpNLO/I")
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
    Events.Branch("genWeight", genWeight, "genWeight/F")

    if not isData:
        hGenweights = TH1D("genWeights","genWeights",1,-0.5,0.5)
    
    for fileNum, ntuplesLine in enumerate(ntuplesListFile):
        if fileNum + 1 > numNtuples: break
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

            if isData:
                if not jc.checkJSON(event.luminosityBlock, event.run): continue
            else:
                if isMadgraph:
                    if event.genWeight < 0: continue
                hGenweights.Fill(0, event.genWeight)
                # if hGenweights.GetSumOfWeights() >= hGenweights.GetEntries(): 
                #     print "WARNING: evt #", count, "genwt", event.genWeight, \
                #             "sumw", hGenweights.GetSumOfWeights(), \
                #             "nentries", hGenweights.GetEntries() 

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
                # Electron_inTimeMuon[i] = list(event.Electron_inTimeMuon)[i]
                Electron_ip3d[i] = list(event.Electron_ip3d)[i]
                # Electron_isGlobal[i] = list(event.Electron_isGlobal)[i]
                Electron_isPFcand[i] = list(event.Electron_isPFcand)[i]
                # Electron_isTracker[i] = list(event.Electron_isTracker)[i]
                Electron_jetIdx[i] = list(event.Electron_jetIdx)[i]
                Electron_pdgId[i] = list(event.Electron_pdgId)[i]
                # Electron_looseId[i] = list(event.Electron_looseId)[i]
                # Electron_mediumId[i] = list(event.Electron_mediumId)[i]
                # Electron_mediumPromptId[i] = list(event.Electron_mediumPromptId)[i]
                # Electron_mvaId[i] = ord(list(event.Electron_mvaId)[i])
                # Electron_tightId[i] = list(event.Electron_tightId)[i]
                Electron_genPartFlav[i] = ord(list(event.Electron_genPartFlav)[i])
                Electron_genPartIdx[i] = list(event.Electron_genPartIdx)[i]
                Electron_mt[i] = sqrt(2 * Electron_pt[i] * event.MET_pt * \
                        (1 - cos(Electron_phi[i] - event.MET_phi)))

            lep1_isMu[0] = int(l1Flav == "Muon")
            lep1_index[0] = l1Index
            lep2_isMu[0] = int(l2Flav == "Muon") 
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
    
            nJet[0] = numGoodJets
            nbtag[0] = numBtag
            nbtagLoose[0] = numBtagLoose
            nbtagTight[0] = numBtagTight
    
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
    
            LHE_HT[0] = event.LHE_HT
            LHE_HTIncoming[0] = event.LHE_HTIncoming
            LHE_Njets[0] = ord(event.LHE_Njets)
            LHE_Nb[0] = ord(event.LHE_Nb)
            LHE_Nuds[0] = ord(event.LHE_Nuds)
            LHE_Nglu[0] = ord(event.LHE_Nglu)
            LHE_NpNLO[0] = ord(event.LHE_NpNLO)
            LHE_NpLO[0] = ord(event.LHE_NpLO)
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
            genWeight[0] = event.genWeight

            Events.Fill()
    
    outFile.cd() # cd to outFile to write to it
    Events.Write()
    if not isData: hGenweights.Write()
    
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
