# Defines functions to perform object and event selection.

import numpy as np
from math import sqrt, cos

#--------------------------------------------------------------------------------#

# returns delta R for an event between 2 leptons or a lepton and a jet, given a 
# 2 types (type1, type2) and the indices of the 2 objects within the list of all
# objects of that type (index1, index2). types allowed are "Muon", "Electron", or
# "Jet"
def deltaR(event, type1, index1, type2, index2):
    assert type1=="Muon" or type1=="Electron" or type1=="Jet"
    assert type2=="Muon" or type2=="Electron" or type2=="Jet"
    
    eta1 = list(getattr(event, type1+"_eta"))[index1]
    eta2 = list(getattr(event, type2+"_eta"))[index2]
    phi1 = list(getattr(event, type1+"_phi"))[index1]
    phi2 = list(getattr(event, type2+"_phi"))[index2]
    
    return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2))
    
#--------------------------------------------------------------------------------#

# Selects a pair of leptons from an event (general private helper method). 
# If findingSameFlav is True: selects for pair of muons if muPreference is True, 
# or for pair of electrons if False. 
# If findingSameFlav is False: selects for leading mu and trailing el if
# muPreference is True, or for leading el and trailing mu if False.
# Returns an array of 2 entries where a[0] = l1Index, a[1] = l2Index, where 
# l1 is always the leading lepton (satisfies the higher pt cut), l2 is trailing.
# All the max/min parameters should be specified by selectMuMu/ElEl/MuEl/ElMu.
def __selectLepts(event, findingSameFlav, muPreference, maxL1OkEta, maxL2OkEta, \
        l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz):
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

    # select leading lepton
    l1Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l1count = getattr(event, "n"+l1Flav)
    for i1 in range(l1count):
        pt = list(getattr(event, l1Flav+"_pt"))[i1]
        absDxy = abs(list(getattr(event, l1Flav+"_dxy"))[i1])
        absDz = abs(list(getattr(event, l1Flav+"_dz"))[i1])
        if l1Flav[0] == "M": # Muon:
            iso = list(getattr(event, "Muon_pfRelIso04_all"))[i1]
        else: # Electron 
            iso = list(getattr(event, "Electron_pfRelIso03_all"))[i1]
        absEta = abs(list(getattr(event, l1Flav+"_eta"))[i1])
        if pt > l1MinOkPt and iso < maxOkIso and absEta < maxL1OkEta \
                and ((iso < minFoundIso) or (iso == minFoundIso \
                and pt > maxFoundPt)) and absDxy < maxOkDxy and absDz < maxOkDz:
            minFoundIso = iso
            maxFoundPt = pt
            l1Index = i1
    if l1Index == -1: return None

    # select trailing lepton
    l2Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l2count = getattr(event, "n"+l2Flav)
    for i2 in range(l2count):
        absDxy = abs(list(getattr(event, l2Flav+"_dxy"))[i2])
        absDz = abs(list(getattr(event, l2Flav+"_dz"))[i2])
        if l1Flav[0] == "M": # Muon:
            iso = list(getattr(event, "Muon_pfRelIso04_all"))[i1]
        else: # Electron 
            iso = list(getattr(event, "Electron_pfRelIso03_all"))[i1]
        absEta = abs(list(getattr(event, l2Flav+"_eta"))[i2])
        if pt > l2MinOkPt and absEta < maxL2OkEta \
                and iso < maxOkIso and absDxy < maxOkDxy and absDz < maxOkDz:
            pt = list(getattr(event, l2Flav+"_pt"))[i2]
            if (iso < minFoundIso) or (iso == minFoundIso and pt > maxFoundPt):
                minFoundIso = iso
                maxFoundPt = pt
                l2Index = i2
    if l2Index == -1: return None

    # print
    # print "l1 pt: " + str(list(getattr(event, l1Flav+"_pt"))[l1Index])
    # print "l1 eta: " + str(list(getattr(event, l1Flav+"_eta"))[l1Index])
    # print "l1 relIso: " + str(list(getattr(event, l1Flav+"_relIso"))[l1Index])

    # print "l2 pt: " + str(list(getattr(event, l2Flav+"_pt"))[l2Index])
    # print "l2 eta: " + str(list(getattr(event, l2Flav+"_eta"))[l2Index])
    # print "l2 relIso: " + str(list(getattr(event, l2Flav+"_relIso"))[l2Index])

    return [l1Index, l2Index]

#--------------------------------------------------------------------------------#
# The following methods all call one of the channels of __selectLepts with default
# max eta, min pt, and max iso values.

def selectMuMu(event, maxL1OkEta=2.4, maxL2OkEta=2.4, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, True, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElEl(event, maxL1OkEta=1.6, maxL2OkEta=1.6, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, True, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectMuEl(event, maxL1OkEta=2.4, maxL2OkEta=1.6, l1MinOkPt=12, \
        l2MinOkPt=15, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, False, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElMu(event, maxL1OkEta=1.6, maxL2OkEta=2.4, l1MinOkPt=25, \
        l2MinOkPt=5, maxOkIso=1, maxOkDxy=1, maxOkDz=1):
    return __selectLepts(event, False, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

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
