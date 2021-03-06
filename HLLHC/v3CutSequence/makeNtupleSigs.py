#!/usr/bin/env python

# NOTE: NEEDS 3 CMD LINE ARGS with values {0 (false) or 1 (true)}: 
# testMode, findingSameFlavor, muPreference 
# Outputs a ROOT file located in ../myData/ containing numSigFiles trees 
# Each tree, tSig{i}, corresponds to a file listed in sig_SingleStop_files,
# contains events that have survived loose dilepton selection cuts, and
# has branches for the same variables as tBkgd outputted by makeNtupleBkgd.py.

print "Importing modules."
import sys
from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
from stopSelection import deltaR,  getNumBtag, findValidJets
from stopSelection import selectMuMu, selectElEl, selectMuEl, selectElMu
import numpy as np
from math import sqrt, cos
from array import array
import time
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

# number of files to process
numSigFiles = 3 # max 25

outDir = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/"

# assemble the outName
outName = outDir+"stopCut_"
if numSigFiles < 10: outName += "0"+str(numSigFiles)
else: outName += str(numSigFiles)
outName += "Sig_"+channelName+".root"

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
jet_pt = np.zeros(20, dtype=np.float32)
jet_eta = np.zeros(20, dtype=np.float32)
jet_phi = np.zeros(20, dtype=np.float32)
jet_ht = array('f',[0.])
# jet_flavour = array('f',[0])
dR_lep1_jet = array('f',[0.])
dR_lep2_jet = array('f',[0.])
njets = array('i',[0])
nbtag = array('i',[0])
nbtagLoose = array('i',[0])
nbtagTight = array('i',[0])
met_pt = array('f',[0.])
met_phi = array('f',[0.])
mt_tot = array('f',[0.])
mt_sum = array('f',[0.])
m_eff = array('f',[0.])
genweight = array('f',[0.])

#--------------------------------------------------------------------------------#
start_time = time.time()
# *************** Filling each signal data in a separate tree  **************
print "Storing variables from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

hGenweightsArr = []

for fileNum, line in enumerate(sigDataListFile):
    if fileNum + 1 > numSigFiles: break

    # SET UP THE OUTPUT TREE
    tSig = TTree("tSig"+str(fileNum), "SUSY stop cut events")
    tSig.Branch("muon_count", muon_count, "muon_count/i")
    tSig.Branch("muon_pt", muon_pt, "muon_pt[20]/F")
    tSig.Branch("muon_eta", muon_eta, "muon_eta[20]/F")
    tSig.Branch("muon_phi", muon_phi, "muon_phi[20]/F")
    tSig.Branch("muon_relIso", muon_relIso, "muon_relIso[20]/F")
    tSig.Branch("muon_charge", muon_charge, "muon_charge[20]/F")
    tSig.Branch("muon_mt", muon_mt, "muon_mt[20]/F")
    tSig.Branch("electron_count", electron_count, "electron_count/i")
    tSig.Branch("electron_pt", electron_pt, "electron_pt[20]/F")
    tSig.Branch("electron_eta", electron_eta, "electron_eta[20]/F")
    tSig.Branch("electron_phi", electron_phi, "electron_phi[20]/F")
    tSig.Branch("electron_relIso", electron_relIso, "electron_relIso[20]/F")
    tSig.Branch("electron_charge", electron_charge, "electron_charge[20]/F")
    tSig.Branch("electron_mt", electron_mt, "electron_mt[20]/F")
    tSig.Branch("found3rdLept", found3rdLept, "found3rdLept/i")
    tSig.Branch("lep1_isMu", lep1_isMu, "lep1_isMu/i")
    tSig.Branch("lep1_index", lep1_index, "lep1_index/i")
    tSig.Branch("lep1_pt", lep1_pt, "lep1_pt/F")
    tSig.Branch("lep1_eta", lep1_eta, "lep1_eta/F")
    tSig.Branch("lep1_phi", lep1_phi, "lep1_phi/F")
    tSig.Branch("lep1_relIso", lep1_relIso, "lep1_relIso/F")
    tSig.Branch("lep1_charge", lep1_charge, "lep1_charge/F")
    tSig.Branch("lep1_mt", lep1_mt, "lep1_mt/F")
    tSig.Branch("lep2_isMu", lep2_isMu, "lep2_isMu/i")
    tSig.Branch("lep2_index", lep2_index, "lep2_index/i")
    tSig.Branch("lep2_pt", lep2_pt, "lep2_pt/F")
    tSig.Branch("lep2_eta", lep2_eta, "lep2_eta/F")
    tSig.Branch("lep2_phi", lep2_phi, "lep2_phi/F")
    tSig.Branch("lep2_relIso", lep2_relIso, "lep2_relIso/F")
    tSig.Branch("lep2_charge", lep2_charge, "lep2_charge/F")
    tSig.Branch("lep2_mt", lep2_mt, "lep2_mt/F")
    tSig.Branch("njets", njets, "njets/i")
    tSig.Branch("jet_pt", jet_pt, "jet_pt[20]/F")
    tSig.Branch("jet_eta", jet_eta, "jet_eta[20]/F")
    tSig.Branch("jet_phi", jet_phi, "jet_phi[20]/F")
    tSig.Branch("jet_ht", jet_ht, "jet_ht/F")
    # tSig.Branch("jet_flavour", jet_flavour, "jet_flavour/F")
    tSig.Branch("dR_lep1_jet", dR_lep1_jet, "dR_lep1_jet/F")
    tSig.Branch("dR_lep2_jet", dR_lep2_jet, "dR_lep2_jet/F")
    tSig.Branch("nbtag", nbtag, "nbtag/i")
    tSig.Branch("nbtagLoose", nbtagLoose, "nbtagLoose/i")
    tSig.Branch("nbtagTight", nbtagTight, "nbtagTight/i")
    tSig.Branch("met_pt", met_pt, "met_pt/F")
    tSig.Branch("met_phi", met_phi, "met_phi/F")
    tSig.Branch("mt_tot", mt_tot, "mt_tot/F")
    tSig.Branch("mt_sum", mt_sum, "mt_sum/F")
    tSig.Branch("m_eff", m_eff, "m_eff/F")
    tSig.Branch("genweight", genweight, "genweight/F")

    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename
    inFile = TFile.Open(sigDataDir + filename + ".root", "READ")

    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    nMax = nentries
    if testMode: nMax = 5000 

    hGenweightsArr.append(TH1D("genweights"+str(fileNum),\
            "genweights"+str(fileNum),1,-0.5,0.5))

    # ***** EVERYTHING BELOW THIS LINE MUST MATCH makeNtupleBkgd.py *****
    # ************ BEGIN LOOPING OVER EVENTS **********
    for count, event in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        hGenweightsArr[fileNum].Fill(0,event.genweight)
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
            muon_pt[i] = list(event.muon_pt)[i]
            muon_eta[i] = list(event.muon_eta)[i]
            muon_phi[i] = list(event.muon_phi)[i]
            muon_relIso[i] = list(event.muon_relIso)[i]
            muon_charge[i] = list(event.muon_charge)[i]
            muon_mt[i] = sqrt(2 * muon_pt[i] * event.pfmet_pt * \
                    (1 - cos(muon_phi[i] - event.pfmet_phi)))
    
        assert l2Index > -1
        electron_count[0] = event.electron_count
        for i in range(event.electron_count):
            electron_pt[i] = list(event.electron_pt)[i]
            electron_eta[i] = list(event.electron_eta)[i]
            electron_phi[i] = list(event.electron_phi)[i]
            electron_relIso[i] = list(event.electron_relIso)[i]
            electron_charge[i] = list(event.electron_charge)[i]
            electron_mt[i] = sqrt(2 * electron_pt[i] * event.pfmet_pt * \
                    (1 - cos(electron_phi[i] - event.pfmet_phi)))
    
        lep1_isMu[0] = int(l1Flav == "muon")
        lep1_index[0] = l1Index
        lep1_pt[0] = list(getattr(event, l1Flav+"_pt"))[l1Index]
        lep1_eta[0] = list(getattr(event, l1Flav+"_eta"))[l1Index]
        lep1_phi[0] = list(getattr(event, l1Flav+"_phi"))[l1Index]
        lep1_relIso[0] = list(getattr(event, l1Flav+"_relIso"))[l1Index]
        lep1_charge[0] = list(getattr(event, l1Flav+"_charge"))[l1Index]
        lep1_mt[0] = sqrt(2 * lep1_pt[0] * event.pfmet_pt * \
                    (1 - cos(lep1_phi[0] - event.pfmet_phi)))
        lep2_isMu[0] = int(l2Flav == "muon") 
        lep2_index[0] = l2Index
        lep2_pt[0] = list(getattr(event, l2Flav+"_pt"))[l2Index]
        lep2_eta[0] = list(getattr(event, l2Flav+"_eta"))[l2Index]
        lep2_phi[0] = list(getattr(event, l2Flav+"_phi"))[l2Index]
        lep2_relIso[0] = list(getattr(event, l2Flav+"_relIso"))[l2Index]
        lep2_charge[0] = list(getattr(event, l2Flav+"_charge"))[l2Index]
        lep2_mt[0] = sqrt(2 * lep2_pt[0] * event.pfmet_pt * \
                    (1 - cos(lep2_phi[0] - event.pfmet_phi)))
    
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
            iMaxPtJ = jets[0] # =index within the AC1B jets arr
            for j in range(numGoodJets): # =index within the tBkgd jets arr
                jIndex = jets[j] #  =index within the AC1B jets arr
                jet_pt[j] = list(event.pfjet_pt)[jIndex]
                if jet_pt[j] > list(event.pfjet_pt)[iMaxPtJ]:
                    iMaxPtJ = jIndex
                jet_eta[j] = list(event.pfjet_eta)[jIndex]
                jet_phi[j] = list(event.pfjet_phi)[jIndex]
                jet_ht[0] += jet_pt[j]
                # jet_flavour[j] = list(event.pfjet_flavour)[jIndex]

            dR_lep1_jet[0] = deltaR(event, l1Flav, l1Index, "pfjet", iMaxPtJ)
            dR_lep2_jet[0] = deltaR(event, l2Flav, l2Index, "pfjet", iMaxPtJ)
    
        met_pt[0] = event.pfmet_pt
        met_phi[0] = event.pfmet_phi
        mt_tot[0] = sqrt(lep1_mt[0]**2+ lep2_mt[0]**2)
        mt_sum[0] = lep1_mt[0] + lep2_mt[0]
        genweight[0] = event.genweight
        m_eff[0] = met_pt[0] + lep1_pt[0]+ lep2_pt[0]
        if njets[0] > 0: m_eff[0] += jet_ht[0]

        tSig.Fill()

    outFile.cd() # cd to outfile to write to it
    tSig.Write()
    hGenweightsArr[fileNum].Write()

#--------------------------------------------------------------------------------#
print
print int(time.time()-start_time), "secs of processing."

outFile.Close()

# f = TFile.Open(outName, "READ")
# t = f.Get("tSig")
# for event in t:
#     for j in range(event.jet_count):
#         print event.jet_pt[j]
# h = f.Get("bkgd_cutflow")
# h.Sumw2()
# h.Draw()
# raw_input()

print "Finished creating", outName
print "Done."

