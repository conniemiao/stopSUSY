from os import getenv
from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector
# from ROOT import gInterpreter, gSystem # allows C++ include statements
# import ROOT
import numpy as np
from math import sqrt, cos
from array import array
from looseJetiD import looseJetiD

cmsswBase = getenv("CMSSW_BASE")
# ROOT.gInterpreter.ProcessLine('#include "' + cmsswBase + '/src/DesyTauAnalyses/NTupleMaker/interface/Jets.h"')

verbose = 0
printEventStats = False
createRoot = True

filename = "../myData/Stop_175_LSP1_small"
inFile = TFile.Open(filename + ".root")
print filename

inTree = inFile.Get("AC1B")
nentries = inTree.GetEntries()
print("nentries={0:d}".format(nentries))
nMax = nentries
# nMax = 5000
count = -1 # = event currently on

if createRoot:
    outName = "selectedMuEl_" + filename + ".root"
    outFile = TFile(outName, "recreate")
    t = TTree("t", "Information on selected mu, el, jets")
    muon_pt = array('f',[0.])
    t.Branch("muon_pt", muon_pt, "muon_pt/F")
    muon_phi = array('f',[0.])
    t.Branch("muon_phi", muon_phi, "muon_phi/F")
    muon_eta = array('f',[0.])
    t.Branch("muon_eta", muon_eta, "muon_eta/F")
    muon_px = array('f',[0.])
    t.Branch("muon_px", muon_px, "muon_px/F")
    muon_py = array('f',[0.])
    t.Branch("muon_py", muon_py, "muon_py/F")
    muon_pz = array('f',[0.])
    t.Branch("muon_pz", muon_pz, "muon_pz/F")
    muon_charge = array('f',[0.])
    t.Branch("muon_charge", muon_charge, "muon_charge/F")
    muon_relIso = array('f',[0.])
    t.Branch("muon_relIso", muon_relIso, "muon_relIso/F")

    electron_pt = array('f',[0.])
    t.Branch("electron_pt", electron_pt, "electron_pt/F")
    electron_phi = array('f',[0.])
    t.Branch("electron_phi", electron_phi, "electron_phi/F")
    electron_eta = array('f',[0.])
    t.Branch("electron_eta", electron_eta, "electron_eta/F")
    electron_px = array('f',[0.])
    t.Branch("electron_px", electron_px, "electron_px/F")
    electron_py = array('f',[0.])
    t.Branch("electron_py", electron_py, "electron_py/F")
    electron_pz = array('f',[0.])
    t.Branch("electron_pz", electron_pz, "electron_pz/F")
    electron_charge = array('f',[0.])
    t.Branch("electron_charge", electron_charge, "electron_charge/F")
    electron_relIso = array('f',[0.])
    t.Branch("electron_relIso", electron_relIso, "electron_relIso/F")

    njets = array('i',[0])
    t.Branch("njets", njets, "njets/I")
    nbtag = array('i',[0])
    t.Branch("nbtag", nbtag, "nbtag/I")
    nbtagTight = array('i',[0])
    t.Branch("nbtagTight", nbtagTight, "nbtagTight/I")

    pfmet_corr_x = array('f',[0.])
    t.Branch("pfmet_corr_x", pfmet_corr_x, "pfmet_corr_x/F")
    pfmet_corr_y = array('f',[0.])
    t.Branch("pfmet_corr_y", pfmet_corr_y, "pfmet_corr_y/F")

    primvertex_count = array('i',[0])
    t.Branch("primvertex_count", primvertex_count, "primvertex_count/I")
    genweight = array('f',[0.])
    t.Branch("genweight", genweight, "genweight/F")

    mtmu = array('f',[0.])
    t.Branch("mtmu", mtmu, "mtmu/F")
    mtel = array('f',[0.])
    t.Branch("mtel", mtel, "mtel/F")

########## CUT VALUES ##########  
zVertexCut = 25
ndofVertexCut = 4
dVertexCut = 2
dRleptonsCut = 0.3

ptElectronLowCut = 13
ptElectronHighCut = 24
etaElectronCut = 2.5
dxyElectronCut = 0.045
dzElectronCut = 0.2        

ptMuonLowCut = 10
ptMuonHighCut = 24
etaMuonCut = 2.4
dxyMuonCut = 0.045
dzMuonCut = 0.2 

isoMuMin = 10**(10)
isoEleMin = 10**(10)

ptVetoElectronCut = 10
etaVetoElectronCut = 2.5
dxyVetoElectronCut = 0.045
dzVetoElectronCut = 0.2
isoVetoElectronCut = 0.3

ptVetoMuonCut = 10
etaVetoMuonCut = 2.4
dxyVetoMuonCut = 0.045
dzVetoMuonCut = 0.2
isoVetoMuonCut = 0.3

etaJetCut = 2.4
ptJetCut = 20
bJetEtaCut = 2.4
bTag = 0.8484
looseBTagCut = 0.2217
medBTagCut = 0.6321
tightBTagCut = 0.8953

drMax = 0.5

def deltaR(eta1, phi1, eta2, phi2):
    return sqrt((eta1-eta2)**2 + (phi1-phi2)**2)

########## B TAGGING SET UP ########## 
# fileTagging = TFile(cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/data/"+\
#         "tagging_efficiencies_ichep2016.root")
# tagEff_B = fileTagging.Get("btag_eff_b")
# tagEff_C = fileTagging.Get("btag_eff_c")
# tagEff_Light = fileTagging.Get("btag_eff_oth")


########## BEGIN LOOPING OVER EVENTS ########## 
pairCount = 0
madeItCount = 0
for entry in inTree:
    if count > nMax: break
    if not printEventStats and (verbose > 1 or (count % 5000 == 0)): 
        print "count={0:d}".format(count)
    count += 1

    # require a good primary vertex
    if abs(entry.primvertex_z) > zVertexCut: continue
    if entry.primvertex_ndof < ndofVertexCut: continue
    dVertex = entry.primvertex_x * entry.primvertex_x + \
            entry.primvertex_y * entry.primvertex_y
    if dVertex > dVertexCut: continue
    if entry.primvertex_count < 2: continue

    # require a collection of good electrons
    electrons = []
    for ie in range(entry.electron_count):
        if entry.electron_pt[ie] < ptElectronLowCut: continue
        if abs(entry.electron_eta[ie]) > etaElectronCut: continue
        if abs(entry.electron_dxy[ie]) > dxyElectronCut: continue
        if abs(entry.electron_dz[ie]) > dzElectronCut: continue
        if not entry.electron_mva_wp80_general_Spring16_v1[ie]: continue
        if not entry.electron_pass_conversion[ie]: continue
        if ord(entry.electron_nmissinginnerhits[ie]) > 1 : continue
        madeItCount += 1 
        if abs(entry.electron_charge[ie]) != 1: continue
        electrons.append(ie)
    if len(electrons) == 0: continue

    # require a collection of good muons
    muons = []
    for im in range(entry.muon_count):
        # if entry.muon_isDuplicate[im]: continue
        # if entry.muon_isBad[im]: continue
        if entry.muon_pt[im] < ptMuonLowCut: continue
        if abs(entry.muon_eta[im]) > etaMuonCut: continue
        if abs(entry.muon_dxy[im]) > dxyMuonCut: continue
        if abs(entry.muon_dz[im]) > dzMuonCut: continue
        if not entry.muon_isMedium[im]: continue
        if abs(entry.muon_charge[im]) != 1: continue
        muons.append(im)
    if len(muons) == 0: continue

    mu_index, el_index = -1, -1 
    isoMuMin, isoEleMin = 10**(10), 10**(10)
    for im in range(len(muons)):
        mIndex = muons[im]

        neutralHadIsoMu = entry.muon_r04_sumNeutralHadronEt[mIndex]
        photonIsoMu = entry.muon_r04_sumPhotonEt[mIndex]
        chargedHadIsoMu = entry.muon_r04_sumChargedHadronPt[mIndex]
        puIsoMu = entry.muon_r04_sumPUPt[mIndex]
        neutralIsoMu = max(0, neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu)
        relIsoMu = (chargedHadIsoMu + neutralIsoMu)/entry.muon_pt[mIndex]

        # apply trigger muon cuts
        isMu23 = False
        isMu8 = False
        if entry.muon_pt[mIndex] > ptMuonHighCut: isMu23 = True
        if entry.muon_pt[mIndex] > ptMuonLowCut: isMu8 = True
        if (not isMu23) and (not isMu8): continue

        for ie in range(len(electrons)):
            eIndex = electrons[ie]

            neutralHadIsoEle = entry.electron_r03_sumNeutralHadronEt[eIndex]
            photonIsoEle = entry.electron_r03_sumPhotonEt[eIndex]
            chargedHadIsoEle = entry.electron_r03_sumChargedHadronPt[eIndex]
            puIsoEle = entry.electron_r03_sumPUPt[eIndex]
            neutralIsoEle = max(0, neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle)
            relIsoEle = (chargedHadIsoEle + neutralIsoEle)/\
                    entry.electron_pt[eIndex]

            # require leptons are not overlapping
            dR = deltaR(entry.electron_eta[eIndex], entry.electron_phi[eIndex], \
                    entry.muon_eta[mIndex], entry.muon_phi[mIndex])
            if dR < dRleptonsCut: continue

            # apply trigger electron cuts
            isEle23 = False
            isEle12 = False
            if entry.electron_pt[eIndex] > ptElectronHighCut: isEle23 = True
            if entry.electron_pt[eIndex] > ptElectronLowCut: isEle12 = True
            trigMatchA = isMu23 and isEle12
            trigMatchB = isMu8 and isEle23
            if (not isEle23) and (not isEle12): continue
            if (not trigMatchA) and (not trigMatchB): continue
        
            # require leptons to have minimum isolation value to indicate signal 
            # (find leptons with minimum isolation). lines 1304-1307 etc. not
            # implemented
            if mIndex != mu_index:
                if relIsoMu < isoMuMin or (relIsoMu == isoMuMin and \
                        entry.muon_pt[mIndex] > entry.muon_pt[mu_index]):
                    isoMuMin = relIsoMu
                    mu_index = mIndex
                    isoEleMin = relIsoEle
                    el_index = eIndex
            elif relIsoEle < isoEleMin or (relIsoEle == isoEleMin and \
                    entry.electron_pt[eIndex] > entry.electron_pt[el_index]):
                isoEleMin = relIsoEle
                el_index = eIndex

    if el_index < 0 or mu_index < 0: continue
    if isoEleMin > 0.3 or isoMuMin > 0.3: continue
    # q = entry.electron_charge[el_index] * entry.muon_charge[mu_index] 

    # look for an extra electron
    foundExtraElectron = False
    for ie in range(entry.electron_count):
        if ie == el_index: continue
        if entry.electron_pt[ie] < ptVetoElectronCut: continue
        if abs(entry.electron_eta[ie]) > etaVetoElectronCut: continue
        if abs(entry.electron_dxy[ie]) > dxyVetoElectronCut: continue
        if abs(entry.electron_dz[ie]) > dzVetoElectronCut: continue
        if not entry.electron_mva_wp90_general_Spring16_v1[ie]: continue
        if not entry.electron_pass_conversion[ie]: continue
        if ord(entry.electron_nmissinginnerhits[ie]) > 1: continue
        neutralHadIsoEle = entry.electron_r03_sumNeutralHadronEt[ie]
        photonIsoEle = entry.electron_r03_sumPhotonEt[ie]
        chargedHadIsoEle = entry.electron_r03_sumChargedHadronPt[ie]
        puIsoEle = entry.electron_r03_sumPUPt[ie]
        neutralIsoEle = max(0, neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle)
        relIsoEle = (chargedHadIsoEle + neutralIsoEle)/entry.electron_pt[ie]
        if relIsoEle > isoVetoElectronCut: continue
        foundExtraElectron = True

    # look for an extra muon
    foundExtraMuon = False
    for im in range(entry.muon_count):
        if im == mu_index: continue
        # if entry.muon_isDuplicate[im]: continue
        # if entry.muon_isBad[im]: continue
        if entry.muon_pt[im] < ptVetoMuonCut: continue
        if abs(entry.muon_eta[im]) > etaVetoMuonCut: continue
        if abs(entry.muon_dxy[im]) > dxyVetoMuonCut: continue
        if abs(entry.muon_dz[im]) > dzVetoMuonCut: continue
        if not entry.muon_isMedium[im]: continue
        neutralHadIsoMu = entry.muon_r04_sumNeutralHadronEt[im]
        photonIsoMu = entry.muon_r04_sumPhotonEt[im]
        chargedHadIsoMu = entry.muon_r04_sumChargedHadronPt[im]
        puIsoMu = entry.muon_r04_sumPUPt[im]
        neutralIsoMu = max(0, neutralHadIsoMu + photonIsoMu - 0.5 * puIsoMu)
        relIsoMu = (chargedHadIsoMu + neutralIsoMu)/entry.muon_pt[im]
        if relIsoMu > isoVetoMuonCut: continue
        foundExtraMuon = True
 
    if foundExtraElectron or foundExtraMuon: continue
    # this continue actually seems to be commented out in line 1435?

    ############### JETS ###############
    # cleans jets
    bjets = []
    bjetsTight = []
    jets = []
    for jet in range(entry.pfjet_count):
        if abs(entry.pfjet_pt[jet]) < ptJetCut: continue
        absJetEta = abs(entry.pfjet_eta[jet])
        if absJetEta > etaJetCut: continue
        
        isPfJetId = looseJetiD(entry, jet)
        if not isPfJetId: continue
        cleanedJet = True
 
        drMu = deltaR(entry.muon_eta[mu_index], entry.muon_phi[mu_index],\
                entry.pfjet_eta[jet], entry.pfjet_phi[jet])
        if drMu < drMax: continue 
        drE = deltaR(entry.electron_eta[el_index], entry.electron_phi[el_index],\
                entry.pfjet_eta[jet], entry.pfjet_phi[jet])
        if drE < drMax: continue 

        # SKIPPING UPGRADING/DOWNGRADING

        if absJetEta<bJetEtaCut: # jet within b-tagging acceptance
            # print(type(entry.pfjet_btag))
            # need to reshape PyRootBuffer object pfjet_btag into a count x 2 array
            entry_pfjet_btag = np.ndarray((entry.pfjet_count, 2), 'f', entry.pfjet_btag)
            # print count
            # print entry_pfjet_btag.item((jet,0)), entry_pfjet_btag.item((jet,1))
            btaggedTight = entry_pfjet_btag.item((jet, 0))  > tightBTagCut
            btagged = entry_pfjet_btag.item((jet, 0)) > bTag
            if btagged:
                bjets.append(jet)
            if btaggedTight:
                bjetsTight.append(jet)
        jets.append(jet)

    ########## CREATE ROOT ##########
    if createRoot:
        muon_pt[0] = entry.muon_pt[mu_index]
        muon_phi[0] = entry.muon_phi[mu_index]
        muon_eta[0] = entry.muon_eta[mu_index]
        muon_px[0] = entry.muon_px[mu_index]
        muon_py[0] = entry.muon_px[mu_index]
        muon_pz[0] = entry.muon_pz[mu_index]
        muon_charge[0] = entry.muon_charge[mu_index]
        muon_relIso[0] = isoMuMin

        electron_pt[0] = entry.electron_pt[el_index]
        electron_phi[0] = entry.electron_phi[el_index]
        electron_eta[0] = entry.electron_eta[el_index]
        electron_px[0] = entry.electron_px[el_index]
        electron_py[0] = entry.electron_px[el_index]
        electron_pz[0] = entry.electron_pz[el_index]
        electron_charge[0] = entry.electron_charge[el_index]
        electron_relIso[0] = isoEleMin

        njets[0] = len(jets)
        nbtag[0] = len(bjets)
        nbtagTight[0] = len(bjetsTight)

        pfmet_corr_x[0] = entry.pfmetcorr_ex
        pfmet_corr_y[0] = entry.pfmetcorr_ey
        primvertex_count[0] = entry.primvertex_count
        genweight[0] = entry.genweight

        mtmu[0] = sqrt(2 * muon_pt[0] * entry.pfmet_pt * \
                (1 - cos(muon_phi[0] - entry.pfmet_phi)))
        mtel[0] = sqrt(2 * electron_pt[0] * entry.pfmet_pt * \
                (1 - cos(electron_phi[0] - entry.pfmet_phi)))

        t.Fill()

    pairCount += 1
    if printEventStats:
        print "Entry %d" % count
        print("  muon_pt: {0}   mu_relIso: {1}   el_pt: {2}   el_relIso: {3}"\
                .format(entry.muon_pt[mu_index], isoMuMin, \
                entry.electron_pt[el_index], isoEleMin))
        print

# print "made it: %d" % madeItCount
print "Number of events with a valid l/l pair found: %d" % pairCount
if createRoot:
    outFile.Write()
    outFile.Close()

print "Done."
