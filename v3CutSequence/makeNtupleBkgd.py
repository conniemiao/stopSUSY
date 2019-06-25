#!/usr/bin/env python

# NOTE: NEEDS 3 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, findingSameFlavor, muPreference 
# Outputs a ROOT file located in ../myData/ containing 1 tree which contains 
# the events from all files listed in bkgd_TTDiLept_file that have survived 
# loose dilepton selection cuts.
# The tree, tBkgd, has branches for the same variables as the tSig{i} outputted
# by makeNtupleSigs.py.

import sys
from ROOT import TFile, TTree, TH1F, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
import numpy as np
from math import sqrt, cos
from array import array

assert len(sys.argv) == 4, "need 3 command line args: testMode{0,1}, findingSameFlavor{0,1}, muPreference{0,1}"

# limits the number of events and files to loop over
testMode = bool(int(sys.argv[1]))
# applying cuts
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
findingSameFlavor = bool(int(sys.argv[2]))
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
muPreference = bool(int(sys.argv[3]))

print "Test mode:", testMode

if findingSameFlavor:
    if muPreference: 
        l1Flav = "muon"
        l2Flav = "muon"
        print "Selecting for pair of 2 muons."
    else: 
        l1Flav = "electron"
        l2Flav = "electron"
        print "Selecting for pair of 2 electrons."
else: 
    print "Selecting for pair of opposite flavor leptons."
    # these 2 lines just matter for creating the outFile name; actual 
    # selection of leading/trailing flavors occurs when looping over events:
    l1Flav = "muon"
    l2Flav = "electron"

# number of files to process
numBkgdFiles = 27  # need to loop over all the files in order to have correct xsec
if testMode: 
    numBkgdFiles = 2 

outDir = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/"

# assemble the outName
outName = outDir+"stopCut_"
if numBkgdFiles < 10: outName += "0"+str(numBkgdFiles)+"Bkgd_TTDiLept_"
else: outName += str(numBkgdFiles)+"Bkgd_TTDiLept_"
outName += l1Flav[:2]+l2Flav[:2]+".root"

outFile = TFile(outName, "recreate")

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
jet_pt = np.zeros(20, dtype=np.float32)
jet_eta = np.zeros(20, dtype=np.float32)
jet_phi = np.zeros(20, dtype=np.float32)
# jet_flavour = array('f',[0])
njets = array('i',[0])
nbtag = array('i',[0])
nbtagLoose = array('i',[0])
nbtagTight = array('i',[0])
met_pt = array('f',[0.])
met_phi = array('f',[0.])
# genweight = array('f',[0.])
#--------------------------------------------------------------------------------#

# ********************** Filling bkgd data summed together  **********************
print "Storing variables from background."
bkgdDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/DESY_pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/"
bkgdDataListFile = open("bkgd_TTDiLept_files")

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
tBkgd.Branch("njets", njets, "njets/i")
tBkgd.Branch("jet_pt", jet_pt, "jet_pt[20]/F")
tBkgd.Branch("jet_eta", jet_eta, "jet_eta[20]/F")
tBkgd.Branch("jet_phi", jet_phi, "jet_phi[20]/F")
# tBkgd.Branch("jet_flavour", jet_flavour, "jet_flavour/F")
tBkgd.Branch("nbtag", nbtag, "nbtag/i")
tBkgd.Branch("nbtagLoose", nbtagLoose, "nbtagLoose/i")
tBkgd.Branch("nbtagTight", nbtagTight, "nbtagTight/i")
tBkgd.Branch("met_pt", met_pt, "met_pt/F")
tBkgd.Branch("met_phi", met_phi, "met_phi/F")
# tBkgd.Branch("genweight", genweight, "genweight/F")

for fileNum, line in enumerate(bkgdDataListFile):
    if fileNum + 1 > numBkgdFiles: break
    filename = line.rstrip()
    print filename

    inFile = TFile.Open(bkgdDataDir + filename, "READ")
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

        # ****** Loose selection of events with valid lep1, lep2, jets ******
        if findingSameFlavor:
            if muPreference:
                lepIndices = selectMuMu(event, l1MinOkPt=20, maxOkIso=0.3)
            else: lepIndices = selectElEl(event, l1MinOkPt=20, maxOkIso=0.3)
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
        # Save all the leptons' and jets' info for this event if it could possibly
        # contain a good lepton pair.
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

        njets[0] = numGoodJets
        nbtag[0] = numBtag
        nbtagLoose[0] = numBtagLoose
        nbtagTight[0] = numBtagTight

        if numGoodJets > 0:
            iMaxPtJ = jets[0] 
            for j in range(numGoodJets):
                jIndex = jets[j]
                jet_pt[j] = list(event.pfjet_pt)[jIndex]
                if jet_pt[j] > list(event.pfjet_pt)[iMaxPtJ]:
                    iMaxPtJ = jIndex
                jet_eta[j] = list(event.pfjet_eta)[jIndex]
                jet_phi[j] = list(event.pfjet_phi)[jIndex]
                # jet_flavour[j] = list(event.pfjet_flavour)[jIndex]

        met_pt[0] = event.pfmet_pt
        met_phi[0] = event.pfmet_phi
        # genweight[0] = event.genweight

        tBkgd.Fill()

outFile.cd() # cd to outFile to write to it
tBkgd.Write()

#--------------------------------------------------------------------------------#

outFile.Close()

# f = TFile.Open(outName, "READ")
# t = f.Get("tBkgd")
# for event in t:
#     for j in range(event.jet_count):
#         print event.jet_pt[j]
# h = f.Get("bkgd_cutflow")
# h.Sumw2()
# h.Draw()
# raw_input()

print "Finished creating", outName
print "Done."
