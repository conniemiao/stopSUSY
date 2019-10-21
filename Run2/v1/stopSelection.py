# Defines functions to perform object and event selection.

import numpy as np
from math import sqrt, cos

#--------------------------------------------------------------------------------#

# returns delta R for an event between 2 leptons or a lepton and a jet, given a 
# 2 types (type1, type2) and the indices of the 2 objects within the list of all
# objects of that type (index1, index2). types allowed are "Muon", "Electron", or
# "Jet"
def deltaR(event, type1, index1, type2, index2):
    eta1 = list(getattr(event, type1+"_eta"))[index1]
    eta2 = list(getattr(event, type2+"_eta"))[index2]
    phi1 = list(getattr(event, type1+"_phi"))[index1]
    phi2 = list(getattr(event, type2+"_phi"))[index2]
    
    return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2))
    
#--------------------------------------------------------------------------------#

# The following methods all call one of the channels of __selectLepts with default
# max eta, min pt, and max iso values.

def selectMuMu(event, isData, maxL1OkEta=2.4, maxL2OkEta=2.4, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, isData, True, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElEl(event, isData, maxL1OkEta=1.6, maxL2OkEta=1.6, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, isData, True, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectMuEl(event, isData, maxL1OkEta=2.4, maxL2OkEta=1.6, l1MinOkPt=12, \
        l2MinOkPt=15, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, isData, False, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElMu(event, isData, maxL1OkEta=1.6, maxL2OkEta=2.4, l1MinOkPt=25, \
        l2MinOkPt=5, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, isData, False, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

#--------------------------------------------------------------------------------#

# Selects a pair of leptons from an event (general private helper method). 
# If findingSameFlav is True: selects for pair of muons if muPreference is True, 
# or for pair of electrons if False. 
# If findingSameFlav is False: selects for leading mu and trailing el if
# muPreference is True, or for leading el and trailing mu if False.
# isData is used to determine proper trigger checks.
# Returns an array of 2 entries where a[0] = l1Index, a[1] = l2Index, where 
# l1 is always the leading lepton (satisfies the higher pt cut), l2 is trailing.
# All the max/min parameters should be specified by selectMuMu/ElEl/MuEl/ElMu.
def __selectLepts(event, isData, findingSameFlav, muPreference, maxL1OkEta, \
        maxL2OkEta, l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz):
    if findingSameFlav:
        if muPreference: # mumu
            l1Flav = "Muon"
            l2Flav = "Muon"
        else: # elel
            l1Flav = "Electron"
            l2Flav = "Electron"
    else:
        if muPreference: # muel
            l1Flav = "Muon"
            l2Flav = "Electron"
        else: # elmu
            l1Flav = "Electron"
            l2Flav = "Muon"

    l1Index = -1
    l2Index = -1
    minSumIso = 1

    # select leading lepton
    l1count = getattr(event, "n"+l1Flav)
    if findingSameFlav and l1count<2: return None
    for i1 in range(l1count):
        pt1 = list(getattr(event, l1Flav+"_pt"))[i1]
        absDxy1 = abs(list(getattr(event, l1Flav+"_dxy"))[i1])
        absDz1 = abs(list(getattr(event, l1Flav+"_dz"))[i1])
        absEta1 = abs(list(getattr(event, l1Flav+"_eta"))[i1])
        if l1Flav[0] == "M": # Muon:
            iso1 = list(getattr(event, "Muon_pfRelIso04_all"))[i1]
        else: # Electron 
            iso1 = list(getattr(event, "Electron_pfRelIso03_all"))[i1]

        if not(pt1 > l1MinOkPt and iso1 < maxOkIso and absEta1 < maxL1OkEta\
                and absDxy1 < maxOkDxy and absDz1 < maxOkDz): continue

        # select trailing lepton
        l2count = getattr(event, "n"+l2Flav)
        if findingSameFlav: startl2Search = i1+1 # don't want l1 = l2
        else: startl2Search = 0
        for i2 in range(startl2Search, l2count):
            pt2 = list(getattr(event, l2Flav+"_pt"))[i2]
            absDxy2 = abs(list(getattr(event, l2Flav+"_dxy"))[i2])
            absDz2 = abs(list(getattr(event, l2Flav+"_dz"))[i2])
            absEta2 = abs(list(getattr(event, l2Flav+"_eta"))[i2])
            if l2Flav[0] == "M": # Muon:
                iso2 = list(getattr(event, "Muon_pfRelIso04_all"))[i2]
            else: # Electron 
                iso2 = list(getattr(event, "Electron_pfRelIso03_all"))[i2]

            if not(pt2 > l2MinOkPt and iso2 < maxOkIso and absEta2 < maxL2OkEta \
                    and absDxy2 < maxOkDxy and absDz2 < maxOkDz): continue

            # total iso check
            if (iso1+iso2) >= minSumIso: continue
            minSumIso = iso1+iso2

            # trigger checks
            if not passTrigger(event, findingSameFlav, muPreference, i1, isData, \
                    l1Flav): continue
            # only do trailing lepton trigger check for muel/elmu
            if not findingSameFlav:
                if not passTrigger(event, findingSameFlav, muPreference, i2, isData, \
                        l2Flav): continue

            # don't check for opposite charge since needed for regions

            # found a pair!
            l1Index = i1
            l2Index = i2

    if l1Index == -1: return None
    if l2Index == -1: return None

    # print l1Index, l2Index, "n1", l1count, "n2", l2count

    # print
    # print "l1 pt: " + str(list(getattr(event, l1Flav+"_pt"))[l1Index])
    # print "l1 eta: " + str(list(getattr(event, l1Flav+"_eta"))[l1Index])
    # print "l1 relIso: " + str(list(getattr(event, l1Flav+"_relIso"))[l1Index])

    # print "l2 pt: " + str(list(getattr(event, l2Flav+"_pt"))[l2Index])
    # print "l2 eta: " + str(list(getattr(event, l2Flav+"_eta"))[l2Index])
    # print "l2 relIso: " + str(list(getattr(event, l2Flav+"_relIso"))[l2Index])

    return [l1Index, l2Index]

#--------------------------------------------------------------------------------#

# Handles redirecting the pass trigger check to the correct function. flav only
# needed for muel/elmu, where it refers to the actual lepton
# being checked at that point (since both mu and el need to be checked and the
# checks are different depending on if they are leading or trailing).
def passTrigger(event, findingSameFlav, muPreference, index, isData, flav):
    if findingSameFlav:
        if muPreference: return passMu_mumuTrigger(event, index) # mumu
        else: return passEl_elelTrigger(event, index) # elel
    else:
        if muPreference: # leading mu, trailing el
            if flav[0] == "M": return passMu_muelTrigger(event, index, isData)
            else: return passEl_muelTrigger(event, index, isData)
        else: # leading el, trailing mu
            if flav[0] == "M": return passMu_elmuTrigger(event, index, isData)
            else: return passEl_elmuTrigger(event, index, isData)

# Returns true if the muon at muIndex passes the mumu trigger check.
def passMu_mumuTrigger(event, muIndex):
    if not event.HLT_IsoMu24: return False
    if list(event.Muon_pt)[muIndex] < 26: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Muon", muIndex, "TrigObj", iObj)
        # filter bit 2: iso
        if event.TrigObj_filterBits[iObj] & 2 and dR < 0.5: return True
    return False

# Returns true if the electron at elIndex passes the elel trigger check.
def passEl_elelTrigger(event, elIndex):
    if not event.HLT_Ele25_eta2p1_WPTight_Gsf: return False
    if list(event.Electron_pt)[elIndex] < 27: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Electron", elIndex, "TrigObj", iObj)
        # filter bit 2: el working point tight
        if event.TrigObj_filterBits[iObj] & 2 and dR < 0.5: return True
    return False

# Returns true if the mu at muIndex passes the muon trigger check for leading mu,
# trailing el. 
def passMu_muelTrigger(event, muIndex, isData):
    if isData and event.run >= 278820: # Run G or higher
        if not event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ:
            return False
    else:
        if not event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL: return False
    if list(event.Muon_pt)[muIndex] < 25: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Muon", muIndex, "TrigObj", iObj)
        # filter bit 32: HLT trigger bits for 1 mu 1 el
        if event.TrigObj_filterBits[iObj] & 32 and dR < 0.5: return True
    return False

# Returns true if the el at elIndex passes the muon trigger check for leading mu,
# trailing el. 
def passEl_muelTrigger(event, elIndex, isData):
    if isData and event.run >= 278820: # Run G or higher
        if not event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ:
            return False
    else:
        if not event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL: return False
    if list(event.Electron_pt)[elIndex] < 14: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Electron", elIndex, "TrigObj", iObj)
        # filter bit 32: HLT trigger bits for 1 mu 1 el
        if event.TrigObj_filterBits[iObj] & 32 and dR < 0.5: return True
    return False

# Returns true if the mu at muIndex passes the muon trigger check for leading el,
# trailing mu. 
def passMu_elmuTrigger(event, muIndex, isData):
    if isData and event.run >= 278820: # Run G or higher
        if not event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ:
            return False
    else:
        if not event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL: return False
    if list(event.Muon_pt)[muIndex] < 10: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Muon", muIndex, "TrigObj", iObj)
        # filter bit 32: HLT trigger bits for 1 mu 1 el
        if event.TrigObj_filterBits[iObj] & 32 and dR < 0.5: return True
    return False

# Returns true if the el at elIndex passes the muon trigger check for leading el,
# trailing mu. 
def passEl_elmuTrigger(event, elIndex, isData):
    if isData and event.run >= 278820: # Run G or higher
        if not event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ:
            return False
    else:
        if not event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL: return False
    if list(event.Electron_pt)[elIndex] < 25: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Electron", elIndex, "TrigObj", iObj)
        # filter bit 32: HLT trigger bits for 1 mu 1 el
        if event.TrigObj_filterBits[iObj] & 32 and dR < 0.5: return True
    return False

#--------------------------------------------------------------------------------#

# loops over all jets for this event and returns an array of all the indices
# that contain valid jets, given the flavs and indices of the 2 selected
# leptons (to clean the jets)
def findValidJets(event, l1Flav, l1Index, l2Flav, l2Index):
    jets = []
    numJets = event.nJet
    for j in range(numJets):
        pt = list(event.Jet_pt)[j]
        # pt0 = list(event.Jet_pt)[0]
        # if pt != pt0: 
        #     print "pt",pt,"pt0",pt0,
        # print "pt", pt, "eta", list(event.Jet_eta)[i],
        dRl1 = deltaR(event, l1Flav, l1Index, "Jet", j)
        dRl2 = deltaR(event, l2Flav, l2Index, "Jet", j)
        if pt > 20 and abs(list(event.Jet_eta)[j]) < 2.4\
                and dRl1 > 0.5 and dRl2 > 0.5:
            jets.append(j)
    return jets

#--------------------------------------------------------------------------------#
    
# Given the jets array containing the indices (within all the jets) of the accepted
# jets, returns an array containing the indices (within the good jets) of btagged jets
# fulfilling some strictness (0=loose, 1=medium, 2=tight).
# threshold = [0.5803, 0.8838, 0.9693] # Jet_btagCSSV2
# threshold = [0.1522, 0.4941, 0.8001] # Jet_btagDeepB 
threshold = [0.0521, 0.3033, 0.7489] # Jet_btagDeepFlavB 
def getBtagIndices(event, jets, strictness=1):
    # bool Jet_btag[jet][0, 1, 2]: 0, 1, 2 = passed loose, medium, tight cuts
    # stored as float so > 0.5 = True
    indices = []
    Jet_btag = list(event.Jet_btagDeepFlavB)
    for iGoodJets, iAllJets in enumerate(jets):
        if Jet_btag[iAllJets] > threshold[strictness]: indices.append(iGoodJets)
    return indices 

#--------------------------------------------------------------------------------#

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findingSameFlavor
# parameters falls in region A (same sign, nominal rel iso)
def isRegionA(l1Charge, l2Charge, l1RelIso, l2RelIso, findingSameFlavor):
    if findingSameFlavor: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge > 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findingSameFlavor
# parameters falls in region B (opposite sign, nominal rel iso - signal region)
def isRegionB(l1Charge, l2Charge, l1RelIso, l2RelIso, findingSameFlavor):
    if findingSameFlavor: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge < 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findingSameFlavor
# parameters falls in region C (opposite sign, inverted rel iso)
def isRegionC(l1Charge, l2Charge, l1RelIso, l2RelIso, findingSameFlavor):
    if findingSameFlavor: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge < 0 and l1RelIso > maxRelIso and l2RelIso > maxRelIso \
            and l1RelIso < 2*maxRelIso and l2RelIso < 2*maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findingSameFlavor
# parameters falls in region D (same sign, inverted rel iso)
def isRegionD(l1Charge, l2Charge, l1RelIso, l2RelIso, findingSameFlavor):
    if findingSameFlavor: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge > 0 and l1RelIso > maxRelIso and l2RelIso > maxRelIso \
            and l1RelIso < 2*maxRelIso and l2RelIso < 2*maxRelIso:
        return True
    return False

#--------------------------------------------------------------------------------#
