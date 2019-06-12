#!/usr/bin/env python

# Outputs a ROOT file located in ../myData/ containing 1+numSigFiles trees. 
# Each tree has one branch for each relevant variable in the input ROOT files
# after performing additional SUSY cuts.
# The bkgd output tree tBkgd is the sum of all files listed in bkgd_TTDiLept_files.
# Each signal tree tSig{i} corresponds to a file listed in sig_SingleStop_files.

from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
from stopSelection.py import passesCut
import numpy as np
from math import sqrt, cos
from array import array

testMode = True # limits the number of events and files to loop over 
cutMode = False # applying cuts
print "Test mode: " + str(testMode)
print "Cut mode: " + str(cutMode)

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
if numSigFiles < 10: outName += "0"+str(numSigFiles)+"Sig"
else: outName += str(numSigFiles)+"Sig"
if not cutMode: outName += "_baseline.root"
else: outName += ".root"

outFile = TFile(outName, "recreate")

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
muon_count = array('f',[0.])
# muon_px = array('f',[0.])
# muon_py = array('f',[0.])
# muon_pz = array('f',[0.])
muon_pt = array('f',[0.])
muon_eta = array('f',[0.])
muon_phi = array('f',[0.])
pfjet_count = array('i',[0])
# pfjet_px = array('f',[0])
# pfjet_py = array('f',[0])
# pfjet_pz = array('f',[0])
pfjet_pt = array('f',[0])
pfjet_eta = array('f',[0])
pfjet_phi = array('f',[0])
# pfjet_flavour = array('f',[0])
# pfjet_btag = array('f',[0])
electron_count = array('i',[0])
# electron_px = array('f',[0.])
# electron_py = array('f',[0.])
# electron_pz = array('f',[0.])
electron_pt = array('f',[0.])
electron_eta = array('f',[0.])
electron_phi = array('f',[0.])
# pfmet_ex = array('f',[0.])
# pfmet_ey = array('f',[0.])
# pfmet_ez = array('f',[0.])
pfmet_pt = array('f',[0.])
pfmet_phi = array('f',[0.])
# genweight = array('f',[0.])
mtel = array('f',[0.])
mtmu = array('f',[0.])
#--------------------------------------------------------------------------------#

# ********************** Filling bkgd data summed together  **********************
print "Storing variables from background."
bkgdDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/DESY_pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/"
bkgdDataListFile = open("bkgd_TTDiLept_files")

# SET UP THE OUTPUT TREE
tBkgd = TTree("tBkgd", "SUSY stop cut events")
tBkgd.Branch("muon_count", muon_count, "muon_count/F")
# tBkgd.Branch("muon_px", muon_px, "muon_px/F")
# tBkgd.Branch("muon_py", muon_py, "muon_py/F")
# tBkgd.Branch("muon_pz", muon_pz, "muon_pz/F")
tBkgd.Branch("muon_pt", muon_pt, "muon_pt/F")
tBkgd.Branch("muon_eta", muon_eta, "muon_eta/F")
tBkgd.Branch("muon_phi", muon_phi, "muon_phi/F")

tBkgd.Branch("pfjet_count", pfjet_count, "pfjet_count/i")
# tBkgd.Branch("pfjet_px", pfjet_px, "pfjet_px/F")
# tBkgd.Branch("pfjet_py", pfjet_py, "pfjet_py/F")
# tBkgd.Branch("pfjet_pz", pfjet_pz, "pfjet_pz/F")
tBkgd.Branch("pfjet_pt", pfjet_pt, "pfjet_pt/F")
tBkgd.Branch("pfjet_eta", pfjet_eta, "pfjet_eta/F")
tBkgd.Branch("pfjet_phi", pfjet_phi, "pfjet_phi/F")
# tBkgd.Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour/F")
# tBkgd.Branch("pfjet_btag", pfjet_btag, "pfjet_btag/F")

tBkgd.Branch("electron_count", electron_count, "electron_count/i")
# tBkgd.Branch("electron_px", electron_px, "electron_px/F")
# tBkgd.Branch("electron_py", electron_py, "electron_py/F")
# tBkgd.Branch("electron_pz", electron_pz, "electron_pz/F")
tBkgd.Branch("electron_pt", electron_pt, "electron_pt/F")
tBkgd.Branch("electron_eta", electron_eta, "electron_eta/F")
tBkgd.Branch("electron_phi", electron_phi, "electron_phi/F")

# tBkgd.Branch("pfmet_ex", pfmet_ex, "pfmet_ex/F")
# tBkgd.Branch("pfmet_ey", pfmet_ey, "pfmet_ey/F")
# tBkgd.Branch("pfmet_ez", pfmet_ez, "pfmet_ez/F")
tBkgd.Branch("pfmet_pt", pfmet_pt, "pfmet_pt/F")
tBkgd.Branch("pfmet_phi", pfmet_phi, "pfmet_phi/F")

# tBkgd.Branch("genweight", genweight, "genweight/F")

tBkgd.Branch("mtmu", mtmu, "mtmu/F")
tBkgd.Branch("mtel", mtel, "mtel/F")

for fileNum, line in enumerate(bkgdDataListFile):
    if fileNum + 1 > numSigFiles: break
    filename = line.rstrip()
    print filename

    inFile = TFile.Open(bkgdDataDir + filename, "READ")
    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))

    nMax = nentries
    if testMode: nMax = 5000
    for count, entry in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        # *** Selecting muon, electron, jet, must be same for sig and bkgd. *** 
        mIndex = -1
        maxPt = 0
        for im in range(entry.muon_count):
            pt = list(entry.muon_pt)[im]
            if pt > 25 and pt > maxPt and abs(list(entry.muon_eta)[im]) < 2.4: 
                maxPt = pt
                mIndex = im
        if mIndex == -1: continue

        eIndex = -1
        maxPt = 0
        for ie in range(entry.electron_count):
            pt = list(entry.electron_pt)[ie]
            if pt > 25 and pt > maxPt and abs(list(entry.electron_eta)[ie]) < 2.4:
                maxPt = pt
                eIndex = ie
        if eIndex == -1: continue

        jIndex = -1
        maxPt = 0
        for ij in range(entry.pfjet_count):
            pt = list(entry.pfjet_pt)[ij]
            if pt > 25 and pt > maxPt and abs(list(entry.pfjet_eta)[ij]) < 2.4:
                maxPt = pt
                jIndex = ij
        if jIndex == -1: continue

        # ************* CUTS: must be same as for sig and bkgd ************
        if cutMode:
            if not passesCut(entry): continue

        # *********** STORE THE DATA *************
        # only events that pass all cuts will be stored
        muon_count[0] = entry.muon_count
        assert mIndex > -1
        # muon_px[0] = list(entry.muon_px)[mIndex]
        # muon_py[0] = list(entry.muon_py)[mIndex]
        # muon_pz[0] = list(entry.muon_pz)[mIndex]
        muon_pt[0] = list(entry.muon_pt)[mIndex]
        muon_eta[0] = list(entry.muon_eta)[mIndex]
        muon_phi[0] = list(entry.muon_phi)[mIndex]
        mtmu[0] = sqrt(2 * muon_pt[0] * entry.pfmet_pt * \
                (1 - cos(muon_phi[0] - entry.pfmet_phi)))

        electron_count[0] = entry.electron_count
        assert eIndex > -1
        # electron_px[0] = list(entry.electron_px)[eIndex]
        # electron_py[0] = list(entry.electron_py)[eIndex]
        # electron_pz[0] = list(entry.electron_pz)[eIndex]
        electron_pt[0] = list(entry.electron_pt)[eIndex]
        electron_eta[0] = list(entry.electron_eta)[eIndex]
        electron_phi[0] = list(entry.electron_phi)[eIndex]
        mtel[0] = sqrt(2 * electron_pt[0] * entry.pfmet_pt * \
                (1 - cos(electron_phi[0] - entry.pfmet_phi)))

        pfjet_count[0] = entry.pfjet_count
        assert jIndex > -1
        # pfjet_px[0] = list(entry.pfjet_px)[jIndex]
        # pfjet_py[0] = list(entry.pfjet_py)[jIndex]
        # pfjet_pz[0] = list(entry.pfjet_pz)[jIndex]
        pfjet_pt[0] = list(entry.pfjet_pt)[jIndex]
        pfjet_eta[0] = list(entry.pfjet_eta)[jIndex]
        pfjet_phi[0] = list(entry.pfjet_phi)[jIndex]
        # pfjet_flavour[0] = list(entry.pfjet_flavour)[jIndex]

        # pfmet_ex[0] = entry.pfmet_ex
        # pfmet_ey[0] = entry.pfmet_ey
        # pfmet_ez[0] = entry.pfmet_ez
        pfmet_pt[0] = entry.pfmet_pt
        pfmet_phi[0] = entry.pfmet_phi
        # genweight[0] = entry.genweight

        tBkgd.Fill()

outFile.cd() # cd to outFile to write to it
tBkgd.Write()


#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate root file  ************
print "Storing variables from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

for fileNum, line in enumerate(sigDataListFile):
    if fileNum + 1 > numSigFiles: break

    # SET UP THE OUTPUT TREE
    tSig = TTree("tSig"+str(fileNum), "SUSY stop cut events")
    tSig.Branch("muon_count", muon_count, "muon_count/F")
    # tSig.Branch("muon_px", muon_px, "muon_px/F")
    # tSig.Branch("muon_py", muon_py, "muon_py/F")
    # tSig.Branch("muon_pz", muon_pz, "muon_pz/F")
    tSig.Branch("muon_pt", muon_pt, "muon_pt/F")
    tSig.Branch("muon_eta", muon_eta, "muon_eta/F")
    tSig.Branch("muon_phi", muon_phi, "muon_phi/F")
    
    tSig.Branch("pfjet_count", pfjet_count, "pfjet_count/i")
    # tSig.Branch("pfjet_px", pfjet_px, "pfjet_px/F")
    # tSig.Branch("pfjet_py", pfjet_py, "pfjet_py/F")
    # tSig.Branch("pfjet_pz", pfjet_pz, "pfjet_pz/F")
    tSig.Branch("pfjet_pt", pfjet_pt, "pfjet_pt/F")
    tSig.Branch("pfjet_eta", pfjet_eta, "pfjet_eta/F")
    tSig.Branch("pfjet_phi", pfjet_phi, "pfjet_phi/F")
    # tSig.Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour/F")
    # tSig.Branch("pfjet_btag", pfjet_btag, "pfjet_btag/F")
    
    tSig.Branch("electron_count", electron_count, "electron_count/i")
    # tSig.Branch("electron_px", electron_px, "electron_px/F")
    # tSig.Branch("electron_py", electron_py, "electron_py/F")
    # tSig.Branch("electron_pz", electron_pz, "electron_pz/F")
    tSig.Branch("electron_pt", electron_pt, "electron_pt/F")
    tSig.Branch("electron_eta", electron_eta, "electron_eta/F")
    tSig.Branch("electron_phi", electron_phi, "electron_phi/F")
    
    # tSig.Branch("pfmet_ex", pfmet_ex, "pfmet_ex/F")
    # tSig.Branch("pfmet_ey", pfmet_ey, "pfmet_ey/F")
    # tSig.Branch("pfmet_ez", pfmet_ez, "pfmet_ez/F")
    tSig.Branch("pfmet_pt", pfmet_pt, "pfmet_pt/F")
    tSig.Branch("pfmet_phi", pfmet_phi, "pfmet_phi/F")

    # tSig.Branch("genweight", genweight, "genweight/F")

    tSig.Branch("mtmu", mtmu, "mtmu/F")
    tSig.Branch("mtel", mtel, "mtel/F")

    # ************ BEGIN LOOPING OVER EVENTS **********
    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    xsec = float(xsec)
    print filename
    inFile = TFile.Open(sigDataDir + filename, "READ")

    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))
    nMax = nentries
    if testMode: nMax = 5000
    for count, entry in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        # *** Selecting muon, electron, jet, must be same for sig and bkgd. *** 
        mIndex = -1
        maxPt = 0
        for im in range(entry.muon_count):
            pt = list(entry.muon_pt)[im]
            if pt > 25 and pt > maxPt and abs(list(entry.muon_eta)[im]) < 2.4: 
                maxPt = pt
                mIndex = im
        if mIndex == -1: continue

        eIndex = -1
        maxPt = 0
        for ie in range(entry.electron_count):
            pt = list(entry.electron_pt)[ie]
            if pt > 25 and pt > maxPt and abs(list(entry.electron_eta)[ie]) < 2.4:
                maxPt = pt
                eIndex = ie
        if eIndex == -1: continue

        jIndex = -1
        maxPt = 0
        for ij in range(entry.pfjet_count):
            pt = list(entry.pfjet_pt)[ij]
            if pt > 25 and pt > maxPt and abs(list(entry.pfjet_eta)[ij]) < 2.4:
                maxPt = pt
                jIndex = ij
        if jIndex == -1: continue

        # ************* CUTS: must be same as for sig and bkgd ************
        if cutMode:
            if not passesCut(entry): continue

        # *********** STORE THE DATA *************
        # only events that pass all cuts will be stored
        muon_count[0] = entry.muon_count
        assert mIndex > -1
        # muon_px[0] = list(entry.muon_px)[mIndex]
        # muon_py[0] = list(entry.muon_py)[mIndex]
        # muon_pz[0] = list(entry.muon_pz)[mIndex]
        muon_pt[0] = list(entry.muon_pt)[mIndex]
        muon_eta[0] = list(entry.muon_eta)[mIndex]
        muon_phi[0] = list(entry.muon_phi)[mIndex]
        mtmu[0] = sqrt(2 * muon_pt[0] * entry.pfmet_pt * \
                (1 - cos(muon_phi[0] - entry.pfmet_phi)))

        electron_count[0] = entry.electron_count
        assert eIndex > -1
        # electron_px[0] = list(entry.electron_px)[eIndex]
        # electron_py[0] = list(entry.electron_py)[eIndex]
        # electron_pz[0] = list(entry.electron_pz)[eIndex]
        electron_pt[0] = list(entry.electron_pt)[eIndex]
        electron_eta[0] = list(entry.electron_eta)[eIndex]
        electron_phi[0] = list(entry.electron_phi)[eIndex]
        mtel[0] = sqrt(2 * electron_pt[0] * entry.pfmet_pt * \
                (1 - cos(electron_phi[0] - entry.pfmet_phi)))

        pfjet_count[0] = entry.pfjet_count
        assert jIndex > -1
        # pfjet_px[0] = list(entry.pfjet_px)[jIndex]
        # pfjet_py[0] = list(entry.pfjet_py)[jIndex]
        # pfjet_pz[0] = list(entry.pfjet_pz)[jIndex]
        pfjet_pt[0] = list(entry.pfjet_pt)[jIndex]
        pfjet_eta[0] = list(entry.pfjet_eta)[jIndex]
        pfjet_phi[0] = list(entry.pfjet_phi)[jIndex]
        # pfjet_flavour[0] = list(entry.pfjet_flavour)[jIndex]

        # pfmet_ex[0] = entry.pfmet_ex
        # pfmet_ey[0] = entry.pfmet_ey
        # pfmet_ez[0] = entry.pfmet_ez
        pfmet_pt[0] = entry.pfmet_pt
        pfmet_phi[0] = entry.pfmet_phi
        # genweight[0] = entry.genweight

        tSig.Fill()

    outFile.cd() # cd to outfile to write to it
    tSig.Write()

#--------------------------------------------------------------------------------#

outFile.Close()

# f = TFile.Open(outName, "READ")
# t = f.Get("tBkgd")
# for entry in t:
#     print entry.muon_pt

print "Finished creating " + outName + "\n"
print "Done."

