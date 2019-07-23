#!/usr/bin/env python

# Outputs a ROOT file located in ../myData/ containing 1+numSigFiles trees and
# 1+numSigFiles cutflow hists.
# Each tree has one branch for each relevant variable in the input ROOT files
# after performing additional SUSY cuts.
# The bkgd output tree tBkgd is the sum of all files listed in bkgd_TTDiLept_files.
# Each signal tree tSig{i} corresponds to a file listed in sig_SingleStop_files.

from ROOT import TFile, TTree, TH1F, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
from stopSelection import deltaR, selectLepts, getNumBtag, findValidJets
import numpy as np
from math import sqrt, cos
from array import array

testMode = True # limits the number of events and files to loop over 
cutMode = True # applying cuts
print "Test mode: ", testMode
print "Cut mode: ", cutMode

findingSameFlavor = False 
# selecting for either mu-mu or el-el (as opposed to mu-el or el-mu)
muPreference = True 
# only applies if findingSameFlav; selects for mu-mu as opposed to el-el
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

# max number of files to process
numBkgdFiles = 27  # max 27
numSigFiles = 25  # max 25
if testMode: 
    numBkgdFiles = 2 
    numSigFiles = 2

outDir = "~/private/CMSSW_9_4_9/s2019_SUSY/myData/"

# assemble the outName
outName = outDir+"stopCut_"
if numBkgdFiles < 10: outName += "0"+str(numBkgdFiles)+"Bkgd_TTDiLept_"
else: outName += str(numBkgdFiles)+"Bkgd_TTDiLept_"
if numSigFiles < 10: outName += "0"+str(numSigFiles)+"Sig_"+l1Flav[:2]+l2Flav[:2]
else: outName += str(numSigFiles)+"Sig_"+l1Flav[:2]+l2Flav[:2]
if not cutMode: outName += "_baseline.root"
else: outName += ".root"

outFile = TFile(outName, "recreate")

# number of events surviving after each cut.
cuts = ["no cut","dilepton","no 3rd lepton","deltaR(ll)>0.3","njets<5","nbtag<2",\
        "MET>80"]

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
lep1_pt = array('f',[0.])
lep1_eta = array('f',[0.])
lep1_phi = array('f',[0.])
lep1_relIso = array('f', [0.])
lep2_pt = array('f',[0.])
lep2_eta = array('f',[0.])
lep2_phi = array('f',[0.])
lep2_relIso = array('f', [0.])
njets = array('i',[0])
jet_pt = np.zeros(20, dtype=np.float32)
jet_eta = np.zeros(20, dtype=np.float32)
jet_phi = np.zeros(20, dtype=np.float32)
# jet_flavour = array('f',[0])
nbtag = array('i',[0])
deltaR_lep1_jet = array('f',[0.]) # deltaR(lep1, jet with max pt)
deltaR_lep2_jet = array('f',[0.])
met_pt = array('f',[0.])
met_phi = array('f',[0.])
# genweight = array('f',[0.])
mtlep2 = array('f',[0.])
mtlep1 = array('f',[0.])
#--------------------------------------------------------------------------------#

# ********************** Filling bkgd data summed together  **********************
print "Storing variables from background."
bkgdDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/DESY_pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/"
bkgdDataListFile = open("bkgd_TTDiLept_files")

# SET UP THE OUTPUT TREE
tBkgd = TTree("tBkgd", "SUSY stop cut events")
tBkgd.Branch("lep1_pt", lep1_pt, "lep1_pt/F")
tBkgd.Branch("lep1_eta", lep1_eta, "lep1_eta/F")
tBkgd.Branch("lep1_phi", lep1_phi, "lep1_phi/F")
tBkgd.Branch("lep1_relIso", lep1_relIso, "lep1_relIso/F")
tBkgd.Branch("lep2_pt", lep2_pt, "lep2_pt/F")
tBkgd.Branch("lep2_eta", lep2_eta, "lep2_eta/F")
tBkgd.Branch("lep2_phi", lep2_phi, "lep2_phi/F")
tBkgd.Branch("lep2_relIso", lep2_relIso, "lep2_relIso/F")
tBkgd.Branch("njets", njets, "njets/i")
tBkgd.Branch("jet_pt", jet_pt, "jet_pt[20]/F")
tBkgd.Branch("jet_eta", jet_eta, "jet_eta[20]/F")
tBkgd.Branch("jet_phi", jet_phi, "jet_phi[20]/F")
# tBkgd.Branch("jet_flavour", jet_flavour, "jet_flavour/F")
tBkgd.Branch("nbtag", nbtag, "nbtag/i")
tBkgd.Branch("deltaR_lep1_jet", deltaR_lep1_jet, "deltaR_lep1_jet/F")
tBkgd.Branch("deltaR_lep2_jet", deltaR_lep2_jet, "deltaR_lep2_jet/F")
tBkgd.Branch("met_pt", met_pt, "met_pt/F")
tBkgd.Branch("met_phi", met_phi, "met_phi/F")
# tBkgd.Branch("genweight", genweight, "genweight/F")
tBkgd.Branch("mtlep1", mtlep1, "mtlep1/F")
tBkgd.Branch("mtlep2", mtlep2, "mtlep2/F")

bkgdCutflowHist = TH1F("bkgd_cutflow","bkgd_cutflow", len(cuts), 0, len(cuts))
for i in range(len(cuts)):
    bkgdCutflowHist.GetXaxis().SetBinLabel(i+1, cuts[i])

for fileNum, line in enumerate(bkgdDataListFile):
    if fileNum + 1 > numBkgdFiles: break
    filename = line.rstrip()
    print filename

    inFile = TFile.Open(bkgdDataDir + filename, "READ")
    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))

    nMax = nentries
    if testMode: nMax = 100000

    # ************ BEGIN LOOPING OVER EVENTS **********
    for count, event in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        bkgdCutflowHist.Fill(0) # no cut

        # *** Selecting lep1, lep2, jet, must be same for sig and bkgd. *** 
        if findingSameFlavor:
            lepIndices = selectLepts(event, True, muPreference)
            if lepIndices is None: continue
            bkgdCutflowHist.Fill(1) # dilepton

            # veto checks:
            # event should not give valid muel or elmu pair
            if selectLepts(event, False, True) is not None: continue
            if selectLepts(event, False, False) is not None: continue
            bkgdCutflowHist.Fill(2) # no 3rd lepton 

        else:
            lepIndices = selectLepts(event, False, True)
            l1Flav = "muon"
            l2Flav = "electron"
            if lepIndices is None:
                lepIndices = selectLepts(event, False, False)
                if lepIndices is None: continue
                l1Flav = "electron"
                l2Flav = "muon"
            bkgdCutflowHist.Fill(1) # dilepton

            # veto check: event should not give valid mumu or elel pair
            if selectLepts(event, True, True) is not None: continue
            if selectLepts(event, True, False) is not None: continue
            bkgdCutflowHist.Fill(2) # no 3rd lepton

        l1Index = lepIndices[0]
        l2Index = lepIndices[1]

        if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        bkgdCutflowHist.Fill(3) # deltaR(ll) > 0.3

        jets = findValidJets(event)
        numJets = len(jets)
        
        # ************* CUTS: must be same as for sig and bkgd ************
        if cutMode:
            if numJets > 4: continue
            bkgdCutflowHist.Fill(4) # njets < 5
            if getNumBtag(event, jets) > 1: continue
            bkgdCutflowHist.Fill(5) # nbtag < 2
            if event.pfmet_pt < 80: continue
            bkgdCutflowHist.Fill(6) # MET > 80

        # *********** STORE THE DATA *************
        # only events that pass all cuts will be stored
        assert l1Index > -1
        lep1_pt[0] = list(getattr(event, l1Flav+"_pt"))[l1Index]
        lep1_eta[0] = list(getattr(event, l1Flav+"_eta"))[l1Index]
        lep1_phi[0] = list(getattr(event, l1Flav+"_phi"))[l1Index]
        lep1_relIso[0] = list(getattr(event, l1Flav+"_relIso"))[l1Index]
        mtlep1[0] = sqrt(2 * lep1_pt[0] * event.pfmet_pt * \
                (1 - cos(lep1_phi[0] - event.pfmet_phi)))

        assert l2Index > -1
        lep2_pt[0] = list(getattr(event, l2Flav+"_pt"))[l2Index]
        lep2_eta[0] = list(getattr(event, l2Flav+"_eta"))[l2Index]
        lep2_phi[0] = list(getattr(event, l2Flav+"_phi"))[l2Index]
        lep2_relIso[0] = list(getattr(event, l2Flav+"_relIso"))[l2Index]
        mtlep2[0] = sqrt(2 * lep2_pt[0] * event.pfmet_pt * \
                (1 - cos(lep2_phi[0] - event.pfmet_phi)))

        jet_pt.fill(0)
        jet_eta.fill(0)
        jet_phi.fill(0)
        njets[0] = numJets

        iMaxPtJ = 0
        for j in range(numJets):
            jIndex = jets[j]
            jet_pt[j] = list(event.pfjet_pt)[jIndex]
            if jet_pt[j] > list(event.pfjet_pt)[iMaxPtJ]:
                maxPtJIndex = j
            jet_eta[j] = list(event.pfjet_eta)[jIndex]
            jet_phi[j] = list(event.pfjet_phi)[jIndex]
            # jet_flavour[j] = list(event.pfjet_flavour)[jIndex]
        if numJets > 0:
            deltaR_lep1_jet[0] = deltaR(event, l1Flav, l1Index, "pfjet", iMaxPtJ)
            deltaR_lep2_jet[0] = deltaR(event, l2Flav, l2Index, "pfjet", iMaxPtJ)

        met_pt[0] = event.pfmet_pt
        met_phi[0] = event.pfmet_phi
        # genweight[0] = event.genweight

        tBkgd.Fill()

outFile.cd() # cd to outFile to write to it
tBkgd.Write()
bkgdCutflowHist.Write()

#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate tree  ************
print "Storing variables from signal."
# sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

sigCutflowHists = []
for fileNum, line in enumerate(sigDataListFile):
    if fileNum + 1 > numSigFiles: break

    # SET UP THE OUTPUT TREE
    tSig = TTree("tSig"+str(fileNum), "SUSY stop cut events")
    # tSig.Branch("lep1_px", lep1_px, "lep1_px/F")
    # tSig.Branch("lep1_py", lep1_py, "lep1_py/F")
    # tSig.Branch("lep1_pz", lep1_pz, "lep1_pz/F")
    tSig.Branch("lep1_pt", lep1_pt, "lep1_pt/F")
    tSig.Branch("lep1_eta", lep1_eta, "lep1_eta/F")
    tSig.Branch("lep1_phi", lep1_phi, "lep1_phi/F")
    tSig.Branch("lep1_relIso", lep1_relIso, "lep1_relIso/F")
    tSig.Branch("lep2_pt", lep2_pt, "lep2_pt/F")
    tSig.Branch("lep2_eta", lep2_eta, "lep2_eta/F")
    tSig.Branch("lep2_phi", lep2_phi, "lep2_phi/F")
    tSig.Branch("lep2_relIso", lep2_relIso, "lep2_relIso/F")
    tSig.Branch("njets", njets, "njets/i")
    tSig.Branch("jet_pt", jet_pt, "jet_pt[20]/F")
    tSig.Branch("jet_eta", jet_eta, "jet_eta[20]/F")
    tSig.Branch("jet_phi", jet_phi, "jet_phi[20]/F")
    # tSig.Branch("jet_flavour", jet_flavour, "jet_flavour[20]/F")
    tSig.Branch("nbtag", nbtag, "nbtag/i")
    tSig.Branch("deltaR_lep1_jet", deltaR_lep1_jet, "deltaR_lep1_jet/F")
    tSig.Branch("deltaR_lep2_jet", deltaR_lep2_jet, "deltaR_lep2_jet/F")
    tSig.Branch("met_pt", met_pt, "met_pt/F")
    tSig.Branch("met_phi", met_phi, "met_phi/F")
    # tSig.Branch("genweight", genweight, "genweight/F")
    tSig.Branch("mtlep1", mtlep1, "mtlep1/F")
    tSig.Branch("mtlep2", mtlep2, "mtlep2/F")

    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename
    inFile = TFile.Open(sigDataDir + filename, "READ")

    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    nMax = nentries
    if testMode: nMax = 100000

    sigCutflowHists.append(TH1F("sig_"+filename[19:24]+"_cutflow",\
            "sig_"+filename[19:24]+"_cutflow", len(cuts), 0, len(cuts)))
    for i in range(len(cuts)):
        sigCutflowHists[fileNum].GetXaxis().SetBinLabel(i+1, cuts[i])

    # ************ BEGIN LOOPING OVER EVENTS **********
    for count, event in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        sigCutflowHists[fileNum].Fill(0) # no cuts

        # *** Selecting lep1, lep2, jet, must be same for sig and bkgd. *** 
        if findingSameFlavor:
            lepIndices = selectLepts(event, True, muPreference)
            if lepIndices is None: continue
            sigCutflowHists[fileNum].Fill(1) # dilepton

            # veto checks:
            # event should not give valid muel or elmu pair
            if selectLepts(event, False, True) is not None: continue
            if selectLepts(event, False, False) is not None: continue
            sigCutflowHists[fileNum].Fill(2) # no 3rd lepton

        else:
            lepIndices = selectLepts(event, False, True)
            l1Flav = "muon"
            l2Flav = "electron"
            if lepIndices is None:
                lepIndices = selectLepts(event, False, False)
                if lepIndices is None: continue
                l1Flav = "electron"
                l2Flav = "muon"
            sigCutflowHists[fileNum].Fill(1) # dilepton

            # veto check: event should not give valid mumu or elel pair
            if selectLepts(event, True, True) is not None: continue
            if selectLepts(event, True, False) is not None: continue
            sigCutflowHists[fileNum].Fill(2) # no 3rd lepton

        l1Index = lepIndices[0]
        l2Index = lepIndices[1]

        if deltaR(event, l1Flav, l1Index, l2Flav, l2Index) < 0.3: continue
        sigCutflowHists[fileNum].Fill(3) # deltaR > 0.3

        jets = findValidJets(event)
        numJets = len(jets)
        
        # ************* CUTS: must be same as for sig and bkgd ************
        if cutMode:
            if numJets > 4: continue
            sigCutflowHists[fileNum].Fill(4) # njets < 5
            if getNumBtag(event, jets) > 1: continue
            sigCutflowHists[fileNum].Fill(5) # nbtag < 2
            if event.pfmet_pt < 80: continue
            sigCutflowHists[fileNum].Fill(6) # MET > 80

        # *********** STORE THE DATA *************
        # only events that pass all cuts will be stored
        assert l1Index > -1
        lep1_pt[0] = list(getattr(event, l1Flav+"_pt"))[l1Index]
        lep1_eta[0] = list(getattr(event, l1Flav+"_eta"))[l1Index]
        lep1_phi[0] = list(getattr(event, l1Flav+"_phi"))[l1Index]
        lep1_relIso[0] = list(getattr(event, l1Flav+"_relIso"))[l1Index]
        mtlep1[0] = sqrt(2 * lep1_pt[0] * event.pfmet_pt * \
                (1 - cos(lep1_phi[0] - event.pfmet_phi)))

        assert l2Index > -1
        lep2_pt[0] = list(getattr(event, l2Flav+"_pt"))[l2Index]
        lep2_eta[0] = list(getattr(event, l2Flav+"_eta"))[l2Index]
        lep2_phi[0] = list(getattr(event, l2Flav+"_phi"))[l2Index]
        lep2_relIso[0] = list(getattr(event, l2Flav+"_relIso"))[l2Index]
        mtlep2[0] = sqrt(2 * lep2_pt[0] * event.pfmet_pt * \
                (1 - cos(lep2_phi[0] - event.pfmet_phi)))

        jet_pt.fill(0)
        jet_eta.fill(0)
        jet_phi.fill(0)
        njets[0] = numJets

        iMaxPtJ = 0
        for j in range(numJets):
            jIndex = jets[j]
            jet_pt[j] = list(event.pfjet_pt)[jIndex]
            if jet_pt[j] > list(event.pfjet_pt)[iMaxPtJ]:
                maxPtJIndex = j
            jet_eta[j] = list(event.pfjet_eta)[jIndex]
            jet_phi[j] = list(event.pfjet_phi)[jIndex]
            # jet_flavour[j] = list(event.pfjet_flavour)[jIndex]
        if numJets > 0:
            deltaR_lep1_jet[0] = deltaR(event, l1Flav, l1Index, "pfjet", iMaxPtJ)
            deltaR_lep2_jet[0] = deltaR(event, l2Flav, l2Index, "pfjet", iMaxPtJ)

        met_pt[0] = event.pfmet_pt
        met_phi[0] = event.pfmet_phi
        # genweight[0] = event.genweight

        tSig.Fill()

    outFile.cd() # cd to outfile to write to it
    tSig.Write()
    sigCutflowHists[fileNum].Write()

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

print "Finished creating " + outName + "\n"
print "Done."

