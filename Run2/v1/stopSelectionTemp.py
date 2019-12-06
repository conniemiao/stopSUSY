
# Defines functions to perform object and event selection.

import numpy as np
from math import sqrt, cos, pi, sin, acos

#--------------------------------------------------------------------------------#

# returns delta R for an event between 2 leptons or a lepton and a jet, given a 
# 2 types (type1, type2) and the indices of the 2 objects within the list of all
# objects of that type (index1, index2). types allowed are "Muon", "Electron", or
# "Jet"
def deltaR(event, type1, index1, type2, index2):
    eta1 = list(getattr(event, type1+"_eta"))[index1]
    phi1 = list(getattr(event, type1+"_phi"))[index1]
    eta2 = list(getattr(event, type2+"_eta"))[index2]
    phi2 = list(getattr(event, type2+"_phi"))[index2]
  
    dPhi = min(abs(phi2-phi1), 2.0*pi-abs(phi2-phi1))
    return sqrt(dPhi**2 + (eta2-eta1)**2)
    
# general version of deltaR with values of eta, phi given
def dR(eta1, phi1, eta2, phi2):
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    return sqrt(dPhi**2 + (eta2-eta1)**2)

# another way of calculating dR
def dPhiFrom2P(Px1, Py1, Px2, Py2) :
    prod = Px1*Px2 + Py1*Py2
    mod1 = sqrt(Px1*Px1+Py1*Py1)
    mod2 = sqrt(Px2*Px2+Py2*Py2)
    cosDPhi = prod/(mod1*mod2)
    return acos(cosDPhi)
def dR2(eta1, phi1, eta2, phi2) : 
    Px1 = cos(phi1)
    Py1 = sin(phi1)
    Px2 = cos(phi2)
    Py2 = sin(phi2)
    dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2)
    dEta = eta1 - eta2
    dR = sqrt(dPhi*dPhi+dEta*dEta)
    return dR

#--------------------------------------------------------------------------------#

# The following methods all call one of the channels of __selectLepts with default
# max eta, min pt, and max iso values.

def selectMuMu(event, isData, maxL1OkEta=2.4, maxL2OkEta=2.4, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=0.5, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, True, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElEl(event, isData, maxL1OkEta=2.1, maxL2OkEta=2.1, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=0.5, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, True, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectMuEl(event, isData, maxL1OkEta=2.4, maxL2OkEta=2.1, l1MinOkPt=12, \
        l2MinOkPt=15, maxOkIso=0.5, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, False, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElMu(event, isData, maxL1OkEta=2.1, maxL2OkEta=2.4, l1MinOkPt=25, \
        l2MinOkPt=5, maxOkIso=0.5, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, False, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

#--------------------------------------------------------------------------------#

# Generalized method to determine whether some event's lepton of flavor flav and
# at index i passes cuts determined by the eta, pt, iso, dxy, and dz cuts and also
# passes the appropriate reconstruction flags.
def passCuts(event, flav, i, maxOkEta, minOkPt, maxOkIso, maxOkDxy, maxOkDz):
    pt = list(getattr(event, flav+"_pt"))[i]
    phi = list(getattr(event, flav+"_phi"))[i]
    absDxy = abs(list(getattr(event, flav+"_dxy"))[i])
    absDz = abs(list(getattr(event, flav+"_dz"))[i])
    absEta = abs(list(getattr(event, flav+"_eta"))[i])
    eta = list(getattr(event, flav+"_eta"))[i]
    if flav[0] == "M": # Muon:
        iso = list(getattr(event, "Muon_pfRelIso04_all"))[i]
        isGlobal = list(getattr(event, "Muon_isGlobal"))[i]
        isTracker = list(getattr(event, "Muon_isTracker"))[i]
        #isTight = list(getattr(event, "tightId"))[i]
        isMedium = list(getattr(event, "Muon_mediumId"))[i]
    else: # Electron 
        iso = list(getattr(event, "Electron_pfRelIso03_all"))[i]
        hits = ord(list(getattr(event, "Electron_lostHits"))[i])
        convVeto = list(getattr(event, "Electron_convVeto"))[i]
        mvaId = list(getattr(event, "Electron_mvaFall17V2noIso_WP90"))[i]

    if pt > minOkPt and iso < maxOkIso and absEta < maxOkEta and \
            absDxy < maxOkDxy and absDz < maxOkDz: return True
    return False

#--------------------------------------------------------------------------------#

# Selects a pair of leptons from an event (general private helper method). 
# If findSameFlav is True: selects for pair of muons if muPref is True, 
# or for pair of electrons if False. 
# If findSameFlav is False: selects for leading mu and trailing el if
# muPref is True, or for leading el and trailing mu if False.
# isData is used to determine proper trigger checks.
# Returns an array of length 4, where a[0] = l1Index, a[1] = l2Index, a[2] = an array
# of the indices of all extra 3rd leptons that are muons, a[3] = an array of the
# indices of all extra 3rd leptons that are electrons. a[2] and a[3] will be of
# length 0 either if it's not valid to have extra leptons for that flavor (e.g. 
# don't look for extra electrons if it's mumu), or if no extra leptons were found.
# l1 is always the leading lepton (satisfies the higher pt cut), l2 is trailing.
# All the max/min parameters should be specified by selectMuMu/ElEl/MuEl/ElMu.
def __selectLepts(event, isData, findSameFlav, muPref, maxL1OkEta, \
        maxL2OkEta, l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz):
    if findSameFlav:
        if muPref: # mumu
            l1Flav = "Muon"
            l2Flav = "Muon"
        else: # elel
            l1Flav = "Electron"
            l2Flav = "Electron"
    else:
        if muPref: # muel
            l1Flav = "Muon"
            l2Flav = "Electron"
        else: # elmu
            l1Flav = "Electron"
            l2Flav = "Muon"

    l1Index = -1
    l2Index = -1
    minSumIso = 9999
    # select leading lepton
    l1count = getattr(event, "n"+l1Flav)
    if findSameFlav and l1count<2: return None
    for i1 in range(l1count):
        if not passCuts(event, l1Flav, i1, maxL1OkEta, l1MinOkPt, maxOkIso, maxOkDxy,\
                maxOkDz): continue

        # select trailing lepton
        l2count = getattr(event, "n"+l2Flav)
        if findSameFlav: startl2Search = i1+1 # don't want l1 = l2
        else: startl2Search = 0
        for i2 in range(startl2Search, l2count):
            if not passCuts(event, l2Flav, i2, maxL2OkEta, l2MinOkPt, maxOkIso, \
                    maxOkDxy, maxOkDz): continue

            dr = deltaR(event, l1Flav, i1, l2Flav, i2)
            if dr < 0.3 : continue

            # total iso check; want to keep pair with min total iso
            if l1Flav[0] == "M":
                iso1 = list(getattr(event, "Muon_pfRelIso04_all"))[i1]
            else: iso1 = list(getattr(event, "Electron_pfRelIso03_all"))[i1]
            if l2Flav[0] == "M":
                iso2 = list(getattr(event, "Muon_pfRelIso04_all"))[i2]
            else: iso2 = list(getattr(event, "Electron_pfRelIso03_all"))[i2]
            if (iso1+iso2) >= minSumIso: continue
            minSumIso = iso1+iso2

            # trigger checks: either leading or trailing lep passes
            if passSingleLeptTrig(event, i1, l1Flav) or \
                    passSingleLeptTrig(event, i2, l2Flav) or \
                    passCrossTrig(event, findSameFlav, muPref, i1, isData, l1Flav) \
                    or passCrossTrig(event, findSameFlav, muPref, i2, isData, l2Flav):
                # found a pair!
                l1Index = i1
                l2Index = i2

            # don't check for opposite charge since needed for QCD regions
    if l1Index == -1: return None
    if l2Index == -1: return None

    extraMuIndices = []
    extraElIndices = []
    if findSameFlav:
        if muPref: # mumu
            extraMuIndices = findExtraLeptSameFlav(event, l1Flav, l1Index, l2Index,\
                    maxL1OkEta, maxOkDxy, maxOkDz)
        else: # elel
            extraElIndices = findExtraLeptSameFlav(event, l1Flav, l1Index, l2Index,\
                    maxL1OkEta, maxOkDxy, maxOkDz)
    else:
        if muPref: # muel
            extraMuIndices = findExtraLeptOppFlav(event, l1Flav, l1Index,\
                    maxL1OkEta, maxOkDxy, maxOkDz)
            extraElIndices = findExtraLeptOppFlav(event, l2Flav, l2Index,\
                    maxL2OkEta, maxOkDxy, maxOkDz)
        else: # elmu
            extraElIndices = findExtraLeptOppFlav(event, l1Flav, l1Index,\
                    maxL1OkEta, maxOkDxy, maxOkDz)
            extraMuIndices = findExtraLeptOppFlav(event, l2Flav, l2Index,\
                    maxL2OkEta, maxOkDxy, maxOkDz)

    return [l1Index, l2Index, extraMuIndices, extraElIndices]

#--------------------------------------------------------------------------------#

# Given the indices of a pair of same flav leptons already found, returns an array of
# the indices of all additional leptons of the same flav that satisfy the same eta,
# dxy, and dz cuts, plus a very loose minimum pt cut and no iso cut.
def findExtraLeptSameFlav(event, flav, l1Index, l2Index, maxOkEta, maxOkDxy, maxOkDz):
    nLeptons = getattr(event, "n"+flav)
    extraLeptsIndices = []
    for i3 in range(nLeptons):
        if i3 == l1Index or i3 == l2Index: continue
        if passCuts(event, flav, i3, maxOkEta, 5, 9999, maxOkDxy, maxOkDz):
            extraLeptsIndices.append(i3)
    return extraLeptsIndices

# Given that a pair of opp flav leptons was already found, and the index of one of
# those leptons (which is of type flav), returns an array of the indices of all
# additional leptons of type flav that satisfy the same eta, dxy, and dz cuts,
# plus a very loose minimum pt cut and no iso cut.
def findExtraLeptOppFlav(event, flav, index, maxOkEta, maxOkDxy, maxOkDz):
    nLeptons = getattr(event, "n"+flav)
    extraLeptsIndices = []
    for i3 in range(nLeptons):
        if i3 == index: continue
        if passCuts(event, flav, i3, maxOkEta, 5, 9999, maxOkDxy, maxOkDz):
            extraLeptsIndices.append(i3)
    return extraLeptsIndices

#--------------------------------------------------------------------------------#

# Handles redirecting the pass single lepton trigger check to the correct function. 
def passSingleLeptTrig(event, index, flav):
    if flav[0] == "M": return passMuTrigger(event, index)
    else: return passElTrigger(event, index)

# Handles redirecting the pass cross trigger check to the correct function. flav only
# needed for muel/elmu, where it refers to the actual lepton
# being checked at that point (since checks are different depending on if they are 
# leading or trailing).
def passCrossTrig(event, findSameFlav, muPref, index, isData, flav):
    if findSameFlav:
        if muPref: return passMu_mumuTrigger(event, index) # mumu
        else: return passEl_elelTrigger(event, index) # elel
    else:
        if muPref: # leading mu, trailing el
            if flav[0] == "M": return passMu_muelTrigger(event, index, isData)
            else: return passEl_muelTrigger(event, index, isData)
        else: # leading el, trailing mu
            if flav[0] == "M": return passMu_elmuTrigger(event, index, isData)
            else: return passEl_elmuTrigger(event, index, isData)

# Returns true if the muon at muIndex passes the single muon trigger check.
def passMuTrigger(event, muIndex):
    if (not (event.HLT_IsoMu27 and list(event.Muon_pt)[muIndex] > 29)) and \
            (not (event.HLT_IsoMu24 and list(event.Muon_pt)[muIndex] > 26)):
                return False

    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Muon", muIndex, "TrigObj", iObj)
        # filter bit 2: iso
        if event.TrigObj_filterBits[iObj] & 2 and dR < 0.5: return True
    return False

# Returns true if the electron at elIndex passes the single electron trigger check.
def passElTrigger(event, elIndex):
    if (not (event.HLT_Ele27_eta2p1_WPTight_Gsf and \
            list(event.Electron_pt)[elIndex] > 29)) and \
            (not (event.HLT_Ele25_eta2p1_WPTight_Gsf and \
            list(event.Electron_pt)[elIndex] > 27)):
                return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Electron", elIndex, "TrigObj", iObj)
        # filter bit 2: el working point tight
        if event.TrigObj_filterBits[iObj] & 2 and dR < 0.5: return True
    return False

# Returns true if the el at elIndex passes the double electron trigger check.
def passEl_elelTrigger(event, elIndex):
    if not event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ: return False
    if list(event.Electron_pt)[elIndex] < 14: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Electron", elIndex, "TrigObj", iObj)
        # filter bit 16: HLT trigger bits for 2 el
        if event.TrigObj_filterBits[iObj] & 16 and dR < 0.5: return True
    return False

# Returns true if the mu at muIndex passes the double muon trigger check.
def passMu_mumuTrigger(event, muIndex):
    if not event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ and \
            not event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ: return False
    if list(event.Muon_pt)[muIndex] < 10: return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Muon", muIndex, "TrigObj", iObj)
        # filter bit 16: HLT trigger bits for 2 mu 
        if event.TrigObj_filterBits[iObj] & 16 and dR < 0.5: return True
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
# threshold = [0.0521, 0.3033, 0.7489] # Jet_btagDeepFlavB 
threshold = [0.2217, 0.6321, 0.8953] # Jet_btagDeepCSVB 
def getBtagIndices(event, jets, strictness=1):
    # bool Jet_btag[jet][0, 1, 2]: 0, 1, 2 = passed loose, medium, tight cuts
    # stored as float so > 0.5 = True
    indices = []
    #Jet_btag = list(event.Jet_btagDeepFlavB)
    Jet_btag = list(event.Jet_btagDeepB)
    for iGoodJets, iAllJets in enumerate(jets):
        if Jet_btag[iAllJets] > threshold[strictness]: indices.append(iGoodJets)
    return indices 

#--------------------------------------------------------------------------------#
# Regions for the ABCD method for QCD estimation

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region A (same sign, nominal rel iso)
def isRegionA(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.15
    else: maxRelIso = 0.25
    if l1Charge*l2Charge > 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region B (opposite sign, nominal rel iso; signal region)
def isRegionB(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.15
    else: maxRelIso = 0.25
    if l1Charge*l2Charge < 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region C (opposite sign, inverted rel iso)
def isRegionC(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.15
    else: maxRelIso = 0.25
    if l1Charge*l2Charge < 0 and l1RelIso > maxRelIso and l2RelIso > maxRelIso \
            and l1RelIso < 2*maxRelIso and l2RelIso < 2*maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region D (same sign, inverted rel iso)
def isRegionD(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.15
    else: maxRelIso = 0.25
    if l1Charge*l2Charge > 0 and l1RelIso > maxRelIso and l2RelIso > maxRelIso \
            and l1RelIso < 2*maxRelIso and l2RelIso < 2*maxRelIso:
        return True
    return False

#--------------------------------------------------------------------------------#

# Regions for fake estimation

# sr: Given the flav/index of 1st and 2nd leptons in pair from the baseline
# selection, returns True if the event has nbtag < 2, nJet < 4, no 3rd lepton was
# found, meets signal region iso and charge reqs, and |mll-mZ|>15. Otherwise returns
# False.
def isSR(event, l1Flav, l1Index, l2Flav, l2Index):
    if event.found3rdLept: return False 
    if event.nbtag > 1: return False 
    if event.nJet >= 4: return False

    if l1Flav[0] == l2Flav[0]: maxRelIso = 0.15
    else: maxRelIso = 0.25
    l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
    l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
    l2Charge = list(getattr(event, l2Flav+"_charge"))[l2Index]
    l2RelIso = list(getattr(event, l2Flav+"_relIso"))[l2Index]
    if not (l1Charge*l2Charge < 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso):
        return False
    mll = mass(l1, l2)
    if abs(mll-80) < 15: return False
    return True 

# cr1a: Given the flav/index of 1st lepton in pair and the flavor of the 2nd as 
# found in the baseline selection, returns the index of a 2nd lepton for the event,
# which is the highest pt 3rd lepton found (no iso cuts), if the event has nbtag < 2,
# nJet < 4, meets signal region iso and charge reqs, and |mll-mZ|>15. If the event 
# didn't pass the CR1a selection, returns -1.
def getCR1al2Index(event, l1Flav, l1Index, l2Flav):
    if getattr(event, "nExtra"+l2Flav) == 0: return -1
    if event.nbtag > 1: return -1 
    if event.nJet >= 4: return -1

    if l1Flav[0] == l2Flav[0]: maxRelIso = 0.15
    else: maxRelIso = 0.25
    l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
    l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
    if l1RelIso >= maxRelIso: return -1

    for extraLeptInd in list(getattr(event, "extra"+l2Flav[:2]+"Indices")):
        l2Charge = list(getattr(event, l2Flav+"_charge"))[extraLeptInd]
        mll = mass(l1, l2)
        if l1Charge*l2Charge < 0 and abs(mll-80) > 15:
            return extraLeptInd
    return -1

# cr1b: Given the flav/index of 1st lepton in pair and the flavor of the 2nd as 
# found in the baseline selection, returns the index of a 2nd lepton for the event,
# which is the highest pt 3rd lepton found with inverted iso, if the event has 
# nbtag < 2, nJet < 4, meets signal region iso and charge reqs, and |mll-mZ|>15. 
# If the event didn't pass the CR1b selection, returns -1.
def getCR1bl2Index(event, l1Flav, l1Index, l2Flav):
    if getattr(event, "nExtra"+l2Flav) == 0: return -1
    if event.nbtag > 1: return -1 
    if event.nJet >= 4: return -1

    if l1Flav[0] == l2Flav[0]: maxRelIso = 0.15
    else: maxRelIso = 0.25
    l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
    l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
    if l1RelIso >= maxRelIso: return -1

    for extraLeptInd in list(getattr(event, "extra"+l2Flav[:2]+"Indices")):
        l2Charge = list(getattr(event, l2Flav+"_charge"))[extraLeptInd]
        l2RelIso = list(getattr(event, l1Flav+"_relIso"))[extraLeptInd]
        mll = mass(l1, l2)
        if l1Charge*l2Charge < 0 and \
                (l2RelIso > maxRelIso and l2RelIso < 2 * maxRelIso) and \
                abs(mll-80) > 15:
            return extraLeptInd
    return -1

# cr3: Given the flav/index of 1st and 2nd leptons in pair from the baseline
# selection, returns True if the event has nbtag = 0, nJet < 4, a 3rd lepton was
# found, meets signal region iso and charge reqs, and |mll-mZ|<15. Otherwise returns
# False.
def isCR3(event, l1Flav, l1Index, l2Flav):
    if not event.found3rdLept: return False 
    if event.nbtag > 0: return False 
    if event.nJet >= 4: return False

    if l1Flav[0] == l2Flav[0]: maxRelIso = 0.15
    else: maxRelIso = 0.25
    l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]
    l1RelIso = list(getattr(event, l1Flav+"_relIso"))[l1Index]
    l2Charge = list(getattr(event, l2Flav+"_charge"))[l2Index]
    l2RelIso = list(getattr(event, l2Flav+"_relIso"))[l2Index]
    if not (l1Charge*l2Charge < 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso):
        return False
    mll = mass(l1, l2)
    if abs(mll-80) > 15: return False
    return True 
