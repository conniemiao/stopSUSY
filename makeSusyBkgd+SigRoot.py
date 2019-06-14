#!/usr/bin/env python

# Outputs a ROOT file located in ../myData/ containing 1+numSigFiles trees. 
# Each tree has one branch for each relevant variable in the input ROOT files
# after performing additional SUSY cuts.
# The bkgd output tree tBkgd is the sum of all files listed in bkgd_TTDiLept_files.
# Each signal tree tSig{i} corresponds to a file listed in sig_SingleStop_files.

from ROOT import TFile, TTree, TH1F, TCanvas, TLorentzVector, TImage, TLegend
from ROOT import gSystem, gStyle
from stopSelection import selectLepts, passesBTags, findValidJets
import numpy as np
from math import sqrt, cos
from array import array

testMode = True # limits the number of events and files to loop over 
cutMode = True # applying cuts
print "Test mode: " + str(testMode)
print "Cut mode: " + str(cutMode)

findingSameFlavor = True 
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
cuts = ["no cut", "dilepton", "no 3rd lepton", "njets<4", "nbtag<2"]

#--------------------------------------------------------------------------------#
# ************* Make all the arrays. *************
# lep1_px = array('f',[0.])
# lep1_py = array('f',[0.])
# lep1_pz = array('f',[0.])
lep1_pt = array('f',[0.])
lep1_eta = array('f',[0.])
lep1_phi = array('f',[0.])
pfjet_count = array('i',[0])
# pfjet_px = np.zeros(20, dtype=np.float32)
# pfjet_py = np.zeros(20, dtype=np.float32)
# pfjet_pz = np.zeros(20, dtype=np.float32)
pfjet_pt = np.zeros(20, dtype=np.float32)
pfjet_eta = np.zeros(20, dtype=np.float32)
pfjet_phi = np.zeros(20, dtype=np.float32)
# pfjet_flavour = array('f',[0])
# pfjet_btag = array('f',[0])
# lep2_px = array('f',[0.])
# lep2_py = array('f',[0.])
# lep2_pz = array('f',[0.])
lep2_pt = array('f',[0.])
lep2_eta = array('f',[0.])
lep2_phi = array('f',[0.])
# pfmet_ex = array('f',[0.])
# pfmet_ey = array('f',[0.])
# pfmet_ez = array('f',[0.])
pfmet_pt = array('f',[0.])
pfmet_phi = array('f',[0.])
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
# tBkgd.Branch("lep1_px", lep1_px, "lep1_px/F")
# tBkgd.Branch("lep1_py", lep1_py, "lep1_py/F")
# tBkgd.Branch("lep1_pz", lep1_pz, "lep1_pz/F")
tBkgd.Branch("lep1_pt", lep1_pt, "lep1_pt/F")
tBkgd.Branch("lep1_eta", lep1_eta, "lep1_eta/F")
tBkgd.Branch("lep1_phi", lep1_phi, "lep1_phi/F")

tBkgd.Branch("pfjet_count", pfjet_count, "pfjet_count/i")
# tBkgd.Branch("pfjet_px", pfjet_px, "pfjet_px[20]/F")
# tBkgd.Branch("pfjet_py", pfjet_py, "pfjet_py[20]/F")
# tBkgd.Branch("pfjet_pz", pfjet_pz, "pfjet_pz[20]/F")
tBkgd.Branch("pfjet_pt", pfjet_pt, "pfjet_pt[20]/F")
tBkgd.Branch("pfjet_eta", pfjet_eta, "pfjet_eta[20]/F")
tBkgd.Branch("pfjet_phi", pfjet_phi, "pfjet_phi[20]/F")
# tBkgd.Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour/F")
# tBkgd.Branch("pfjet_btag", pfjet_btag, "pfjet_btag/F")

# tBkgd.Branch("lep2_px", lep2_px, "lep2_px/F")
# tBkgd.Branch("lep2_py", lep2_py, "lep2_py/F")
# tBkgd.Branch("lep2_pz", lep2_pz, "lep2_pz/F")
tBkgd.Branch("lep2_pt", lep2_pt, "lep2_pt/F")
tBkgd.Branch("lep2_eta", lep2_eta, "lep2_eta/F")
tBkgd.Branch("lep2_phi", lep2_phi, "lep2_phi/F")

# tBkgd.Branch("pfmet_ex", pfmet_ex, "pfmet_ex/F")
# tBkgd.Branch("pfmet_ey", pfmet_ey, "pfmet_ey/F")
# tBkgd.Branch("pfmet_ez", pfmet_ez, "pfmet_ez/F")
tBkgd.Branch("pfmet_pt", pfmet_pt, "pfmet_pt/F")
tBkgd.Branch("pfmet_phi", pfmet_phi, "pfmet_phi/F")

# tBkgd.Branch("genweight", genweight, "genweight/F")

tBkgd.Branch("mtlep1", mtlep1, "mtlep1/F")
tBkgd.Branch("mtlep2", mtlep2, "mtlep2/F")

bkgd_cutflow_hist = TH1F("bkgd_cutflow","bkgd_cutflow", len(cuts), 0, len(cuts))
for i in range(len(cuts)):
    bkgd_cutflow_hist.GetXaxis().SetBinLabel(i+1, cuts[i])

for fileNum, line in enumerate(bkgdDataListFile):
    if fileNum + 1 > numSigFiles: break
    filename = line.rstrip()
    print filename

    inFile = TFile.Open(bkgdDataDir + filename, "READ")
    inTree = inFile.Get("AC1B")
    nentries = inTree.GetEntries()
    print("nentries={0:d}".format(nentries))

    nMax = nentries
    if testMode: nMax = 100000
    for count, event in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        bkgd_cutflow_hist.Fill(0) # no cut

        # *** Selecting lep1, lep2, jet, must be same for sig and bkgd. *** 
        if findingSameFlavor:
            lepIndices = selectLepts(event, True, muPreference)
            if lepIndices is None: continue
            bkgd_cutflow_hist.Fill(1) # dilepton

            # veto checks:
            if muPreference: # mumu
                # event should not give valid lead mu, trail el pair
                if not selectLepts(event, False, True) is None: continue
            else: # elel
                # event should not give valid lead el, trail mu pair
                if not selectLepts(event, False, False) is None: continue
            bkgd_cutflow_hist.Fill(2) # no 3rd lepton 

        else:
            lepIndices = selectLepts(event, False, True)
            l1Flav = "muon"
            l2Flav = "electron"
            if lepIndices is None:
                lepIndices = selectLepts(event, False, False)
                if lepIndices is None: continue
                l1Flav = "electron"
                l2Flav = "muon"
            bkgd_cutflow_hist.Fill(1) # dilepton

            # veto check: event should not give valid mumu or elel pair
            if not selectLepts(event, True, True) is None: continue
            if not selectLepts(event, True, False) is None: continue
            bkgd_cutflow_hist.Fill(2) # no 3rd lepton

        l1Index = lepIndices[0]
        l2Index = lepIndices[1]

        jets = findValidJets(event)
        
        # ************* CUTS: must be same as for sig and bkgd ************
        if cutMode:
            if event.pfjet_count >= 4: continue
            bkgd_cutflow_hist.Fill(3) # njets < 4
            if not passesBTags(event, jets): continue
            bkgd_cutflow_hist.Fill(4) # nbtag < 2

        # *********** STORE THE DATA *************
        # only events that pass all cuts will be stored
        assert l1Index > -1
        # lep1_px[0] = list(getattr(event, l1Flav+"_px")[l1Index]
        # lep1_py[0] = list(getattr(event, l1Flav+"_py")[l1Index]
        # lep1_pz[0] = list(getattr(event, l1Flav+"_pz")[l1Index]
        lep1_pt[0] = list(getattr(event, l1Flav+"_pt"))[l1Index]
        lep1_eta[0] = list(getattr(event, l1Flav+"_eta"))[l1Index]
        lep1_phi[0] = list(getattr(event, l1Flav+"_phi"))[l1Index]
        mtlep1[0] = sqrt(2 * lep1_pt[0] * event.pfmet_pt * \
                (1 - cos(lep1_phi[0] - event.pfmet_phi)))

        assert l2Index > -1
        # lep2_px[0] = list(getattr(event, l2Flav+"_px"))[l2Index]
        # lep2_py[0] = list(getattr(event, l2Flav+"_py"))[l2Index]
        # lep2_pz[0] = list(getattr(event, l2Flav+"_pz"))[l2Index]
        lep2_pt[0] = list(getattr(event, l2Flav+"_pt"))[l2Index]
        lep2_eta[0] = list(getattr(event, l2Flav+"_eta"))[l2Index]
        lep2_phi[0] = list(getattr(event, l2Flav+"_phi"))[l2Index]
        mtlep2[0] = sqrt(2 * lep2_pt[0] * event.pfmet_pt * \
                (1 - cos(lep2_phi[0] - event.pfmet_phi)))

        # pfjet_px.fill(0)
        # pfjet_py.fill(0)
        # pfjet_pz.fill(0)
        pfjet_pt.fill(0)
        pfjet_eta.fill(0)
        pfjet_phi.fill(0)
        pfjet_count[0] = len(jets)
        for j in range(pfjet_count[0]):
            jIndex = jets[j]
            # pfjet_px[j] = list(event.pfjet_px)[jIndex]
            # pfjet_py[j] = list(event.pfjet_py)[jIndex]
            # pfjet_pz[j] = list(event.pfjet_pz)[jIndex]
            pfjet_pt[j] = list(event.pfjet_pt)[jIndex]
            pfjet_eta[j] = list(event.pfjet_eta)[jIndex]
            pfjet_phi[j] = list(event.pfjet_phi)[jIndex]
            # pfjet_flavour[j] = list(event.pfjet_flavour)[jIndex]

        # pfmet_ex[0] = event.pfmet_ex
        # pfmet_ey[0] = event.pfmet_ey
        # pfmet_ez[0] = event.pfmet_ez
        pfmet_pt[0] = event.pfmet_pt
        pfmet_phi[0] = event.pfmet_phi
        # genweight[0] = event.genweight

        tBkgd.Fill()

outFile.cd() # cd to outFile to write to it
tBkgd.Write()
bkgd_cutflow_hist.Write()


#--------------------------------------------------------------------------------#
# *************** Filling each signal data in a separate root file  ************
print "Storing variables from signal."
sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v2/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")

sig_cutflow_hists = []
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
    
    tSig.Branch("pfjet_count", pfjet_count, "pfjet_count/i")
    # tSig.Branch("pfjet_px", pfjet_px, "pfjet_px[20]/F")
    # tSig.Branch("pfjet_py", pfjet_py, "pfjet_py[20]/F")
    # tSig.Branch("pfjet_pz", pfjet_pz, "pfjet_pz[20]/F")
    tSig.Branch("pfjet_pt", pfjet_pt, "pfjet_pt[20]/F")
    tSig.Branch("pfjet_eta", pfjet_eta, "pfjet_eta[20]/F")
    tSig.Branch("pfjet_phi", pfjet_phi, "pfjet_phi[20]/F")
    # tSig.Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour[20]/F")
    # tSig.Branch("pfjet_btag", pfjet_btag, "pfjet_btag[20]/F")
    
    # tSig.Branch("lep2_px", lep2_px, "lep2_px/F")
    # tSig.Branch("lep2_py", lep2_py, "lep2_py/F")
    # tSig.Branch("lep2_pz", lep2_pz, "lep2_pz/F")
    tSig.Branch("lep2_pt", lep2_pt, "lep2_pt/F")
    tSig.Branch("lep2_eta", lep2_eta, "lep2_eta/F")
    tSig.Branch("lep2_phi", lep2_phi, "lep2_phi/F")
    
    # tSig.Branch("pfmet_ex", pfmet_ex, "pfmet_ex/F")
    # tSig.Branch("pfmet_ey", pfmet_ey, "pfmet_ey/F")
    # tSig.Branch("pfmet_ez", pfmet_ez, "pfmet_ez/F")
    tSig.Branch("pfmet_pt", pfmet_pt, "pfmet_pt/F")
    tSig.Branch("pfmet_phi", pfmet_phi, "pfmet_phi/F")

    # tSig.Branch("genweight", genweight, "genweight/F")

    tSig.Branch("mtlep1", mtlep1, "mtlep1/F")
    tSig.Branch("mtlep2", mtlep2, "mtlep2/F")


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
    if testMode: nMax = 100000

    sig_cutflow_hists.append(TH1F("sig_"+filename[19:24]+"_cutflow",\
            "sig_"+filename[19:24]+"_cutflow", len(cuts), 0, len(cuts)))
    for i in range(len(cuts)):
        sig_cutflow_hists[fileNum].GetXaxis().SetBinLabel(i+1, cuts[i])

    for count, event in enumerate(inTree):
        if count > nMax : break
        if count % 500000 == 0: print("count={0:d}".format(count))

        sig_cutflow_hists[fileNum].Fill(0) # no cuts

        # *** Selecting lep1, lep2, jet, must be same for sig and bkgd. *** 
        if findingSameFlavor:
            lepIndices = selectLepts(event, True, muPreference)
            if lepIndices is None: continue
            sig_cutflow_hists[fileNum].Fill(1) # dilepton

            # veto checks:
            if muPreference: # mumu
                # event should not give valid lead mu, trail el pair
                if not selectLepts(event, False, True) is None: continue
            else: # elel
                # event should not give valid lead el, trail mu pair
                if not selectLepts(event, False, False) is None: continue
            sig_cutflow_hists[fileNum].Fill(2) # no 3rd lepton

        else:
            lepIndices = selectLepts(event, False, True)
            l1Flav = "muon"
            l2Flav = "electron"
            if lepIndices is None:
                lepIndices = selectLepts(event, False, False)
                if lepIndices is None: continue
                l1Flav = "electron"
                l2Flav = "muon"
            sig_cutflow_hists[fileNum].Fill(1) # dilepton

            # veto check: event should not give valid mumu or elel pair
            if not selectLepts(event, True, True) is None: continue
            if not selectLepts(event, True, False) is None: continue
            sig_cutflow_hists[fileNum].Fill(2) # no 3rd lepton

        l1Index = lepIndices[0]
        l2Index = lepIndices[1]

        jets = findValidJets(event)
        
        # ************* CUTS: must be same as for sig and bkgd ************
        if cutMode:
            if event.pfjet_count >= 4: continue
            sig_cutflow_hists[fileNum].Fill(3) # njets < 4
            if not passesBTags(event, jets): continue
            sig_cutflow_hists[fileNum].Fill(4) # nbtag < 2

        # *********** STORE THE DATA *************
        # only events that pass all cuts will be stored
        assert l1Index > -1
        # lep1_px[0] = list(getattr(event, l1Flav+"_px")[l1Index]
        # lep1_py[0] = list(getattr(event, l1Flav+"_py")[l1Index]
        # lep1_pz[0] = list(getattr(event, l1Flav+"_pz")[l1Index]
        lep1_pt[0] = list(getattr(event, l1Flav+"_pt"))[l1Index]
        lep1_eta[0] = list(getattr(event, l1Flav+"_eta"))[l1Index]
        lep1_phi[0] = list(getattr(event, l1Flav+"_phi"))[l1Index]
        mtlep1[0] = sqrt(2 * lep1_pt[0] * event.pfmet_pt * \
                (1 - cos(lep1_phi[0] - event.pfmet_phi)))

        assert l2Index > -1
        # lep2_px[0] = list(getattr(event, l2Flav+"_px"))[l2Index]
        # lep2_py[0] = list(getattr(event, l2Flav+"_py"))[l2Index]
        # lep2_pz[0] = list(getattr(event, l2Flav+"_pz"))[l2Index]
        lep2_pt[0] = list(getattr(event, l2Flav+"_pt"))[l2Index]
        lep2_eta[0] = list(getattr(event, l2Flav+"_eta"))[l2Index]
        lep2_phi[0] = list(getattr(event, l2Flav+"_phi"))[l2Index]
        mtlep2[0] = sqrt(2 * lep2_pt[0] * event.pfmet_pt * \
                (1 - cos(lep2_phi[0] - event.pfmet_phi)))

        # pfjet_px.fill(0)
        # pfjet_py.fill(0)
        # pfjet_pz.fill(0)
        pfjet_pt.fill(0)
        pfjet_eta.fill(0)
        pfjet_phi.fill(0)
        pfjet_count[0] = len(jets)
        for j in range(pfjet_count[0]):
            jIndex = jets[j]
            # pfjet_px[j] = list(event.pfjet_px)[jIndex]
            # pfjet_py[j] = list(event.pfjet_py)[jIndex]
            # pfjet_pz[j] = list(event.pfjet_pz)[jIndex]
            pfjet_pt[j] = list(event.pfjet_pt)[jIndex]
            pfjet_eta[j] = list(event.pfjet_eta)[jIndex]
            pfjet_phi[j] = list(event.pfjet_phi)[jIndex]
            # pfjet_flavour[j] = list(event.pfjet_flavour)[jIndex]

        # pfmet_ex[0] = event.pfmet_ex
        # pfmet_ey[0] = event.pfmet_ey
        # pfmet_ez[0] = event.pfmet_ez
        pfmet_pt[0] = event.pfmet_pt
        pfmet_phi[0] = event.pfmet_phi
        # genweight[0] = event.genweight

        tSig.Fill()

    outFile.cd() # cd to outfile to write to it
    tSig.Write()
    sig_cutflow_hists[fileNum].Write()

#--------------------------------------------------------------------------------#

outFile.Close()

# f = TFile.Open(outName, "READ")
# t = f.Get("tBkgd")
# for event in t:
#     for j in range(event.pfjet_count):
#         print event.pfjet_pt[j]
# h = f.Get("bkgd_cutflow")
# h.Sumw2()
# h.Draw()
# raw_input()

print "Finished creating " + outName + "\n"
print "Done."

