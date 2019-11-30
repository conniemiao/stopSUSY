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
  
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    #return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2))
    return sqrt(dPhi**2 + (eta2-eta1)**2)
    
def DR(eta1,phi1,eta2,phi2):
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    #return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2))
    return sqrt(dPhi**2 + (eta2-eta1)**2)

def dPhiFrom2P(Px1, Py1, Px2, Py2) :


    prod = Px1*Px2 + Py1*Py2
    mod1 = sqrt(Px1*Px1+Py1*Py1)
    mod2 = sqrt(Px2*Px2+Py2*Py2)
  
    cosDPhi = prod/(mod1*mod2)
  
    return acos(cosDPhi)


def DR2(eta1, phi1, eta2, phi2) : 

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
        l2MinOkPt=-0.01, maxOkIso=1, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, True, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElEl(event, isData, maxL1OkEta=2.1, maxL2OkEta=2.1, l1MinOkPt=20, \
        l2MinOkPt=-0.01, maxOkIso=1, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, True, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectMuEl(event, isData, maxL1OkEta=2.4, maxL2OkEta=2.1, l1MinOkPt=12, \
        l2MinOkPt=15, maxOkIso=1, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, False, True, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

def selectElMu(event, isData, maxL1OkEta=2.1, maxL2OkEta=2.4, l1MinOkPt=25, \
        l2MinOkPt=5, maxOkIso=1, maxOkDxy=0.045, maxOkDz=0.2):
    return __selectLepts(event, isData, False, False, maxL1OkEta, maxL2OkEta, \
            l1MinOkPt, l2MinOkPt, maxOkIso, maxOkDxy, maxOkDz)

#--------------------------------------------------------------------------------#

# Selects a pair of leptons from an event (general private helper method). 
# If findSameFlav is True: selects for pair of muons if muPref is True, 
# or for pair of electrons if False. 
# If findSameFlav is False: selects for leading mu and trailing el if
# muPref is True, or for leading el and trailing mu if False.
# isData is used to determine proper trigger checks.
# Returns an array of 2 entries where a[0] = l1Index, a[1] = l2Index, where 
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
    minSumIso = 999
    isoMin = 999
    conveto1 = False
    conveto2 = False
    mvaID1 = False
    mvaID2 = False
    hits1 = 2
    hits2 = 2
    isGlobal1 = False
    isTracker1 = False
    isGlobal2 = False
    isTracker2 = False
    isMedium1 = False
    isMedium2 = False
    iso1=999
    iso2=999
    pt2=-1
    pt1=-1
    eta1=999
    eta2=999
    sumIso=999
    # select leading lepton
    l1count = getattr(event, "n"+l1Flav)
    if findSameFlav and l1count<2: return None
    for i1 in range(l1count):
        pt1 = list(getattr(event, l1Flav+"_pt"))[i1]
        phi1 = list(getattr(event, l1Flav+"_phi"))[i1]
        absDxy1 = abs(list(getattr(event, l1Flav+"_dxy"))[i1])
        absDz1 = abs(list(getattr(event, l1Flav+"_dz"))[i1])
        absEta1 = abs(list(getattr(event, l1Flav+"_eta"))[i1])
        eta1 = list(getattr(event, l1Flav+"_eta"))[i1]
        if l1Flav[0] == "M": # Muon:
            iso1 = list(getattr(event, "Muon_pfRelIso04_all"))[i1]
            isGlobal1 = list(getattr(event, "Muon_isGlobal"))[i1]
            isTracker1 = list(getattr(event, "Muon_isTracker"))[i1]
            #isTight = list(getattr(event, "tightId"))[i1]
            isMedium1 = list(getattr(event, "Muon_mediumId"))[i1]
        else: # Electron 
            iso1 = list(getattr(event, "Electron_pfRelIso03_all"))[i1]
            hits1 = ord(list(getattr(event, "Electron_lostHits"))[i1])
            conveto1 = list(getattr(event, "Electron_convVeto"))[i1]
            mvaID1 = list(getattr(event, "Electron_mvaFall17V2noIso_WP90"))[i1]

        if not(pt1 > l1MinOkPt and iso1 < maxOkIso and absEta1 < maxL1OkEta\
                and absDxy1 < maxOkDxy and absDz1 < maxOkDz): continue
        #if l1Flav[0] == "M" and not isGlobal1 and not isTracker1 : continue
        if l1Flav[0] == "M" and not isMedium1 : continue
        if l1Flav[0] == "E" :
            if hits1 > 1 or not conveto1 or not mvaID1 : 
                continue

        # select trailing lepton
        l2count = getattr(event, "n"+l2Flav)
        if findSameFlav: startl2Search = i1+1 # don't want l1 = l2
        else: startl2Search = 0
        for i2 in range(startl2Search, l2count):
            if findSameFlav and i2 == i1 : continue
            pt2 = list(getattr(event, l2Flav+"_pt"))[i2]
            phi2 = list(getattr(event, l2Flav+"_phi"))[i2]
            absDxy2 = abs(list(getattr(event, l2Flav+"_dxy"))[i2])
            absDz2 = abs(list(getattr(event, l2Flav+"_dz"))[i2])
            absEta2 = abs(list(getattr(event, l2Flav+"_eta"))[i2])
            eta2 = list(getattr(event, l2Flav+"_eta"))[i2]
            if l2Flav[0] == "M": # Muon:
                iso2 = list(getattr(event, "Muon_pfRelIso04_all"))[i2]
                isGlobal2 = list(getattr(event, "Muon_isGlobal"))[i2]
                isTracker2 = list(getattr(event, "Muon_isTracker"))[i2]
                #isTight = list(getattr(event, "tightId"))[i2]
                isMedium2 = list(getattr(event, "Muon_mediumId"))[i2]

            else: # Electron 
                iso2 = list(getattr(event, "Electron_pfRelIso03_all"))[i2]
                hits2 = ord(list(getattr(event, "Electron_lostHits"))[i2])
                conveto2 = list(getattr(event, "Electron_convVeto"))[i2]
                mvaID2 = list(getattr(event, "Electron_mvaFall17V2noIso_WP90"))[i2]

            if not(pt2 > l2MinOkPt and iso2 < maxOkIso and absEta2 < maxL2OkEta \
                    and absDxy2 < maxOkDxy and absDz2 < maxOkDz): continue
            #if l2Flav[0] == "M" and not isGlobal2 and not isTracker2 : continue
            if l2Flav[0] == "M" and not isMedium2 : continue
            if l2Flav[0] == "E" :
                if hits2 > 1 or not conveto2 or not mvaID2 : 
                    continue

            dr = DR(eta1,phi1,eta2,phi2)
            dr2 = DR2(eta1,phi1,eta2,phi2)
            dr3 = deltaR(event, l1Flav, i1, l2Flav, i2)
            #if dr != dr3 : print 'iso ', iso1, iso2, 'pt1==========', pt1,pt2, 'dr========', dr, 'dr2========', dr2, 'dr3=============', dr3, 'maxOkDz======', maxOkDz, 'maxOkIso', maxOkIso
            #print 'iso ', eta1, eta2, 'pt1==========', pt1,pt2, 'dr========', dr, 'dr2========', dr2, 'dr3=============', dr3, l1Flav, list(getattr(event, l1Flav+"_eta"))[i1]
            
            if dr < 0.3 : continue

            # total iso check and keep the pair with the min iso sum
            sumIso  = iso1+iso2
            if sumIso > isoMin: continue
           
            if sumIso < isoMin :
	        isoMin = sumIso


            # trigger checks: either leading or trailing lep passes
            if passSingleLeptTrig(event, i1, l1Flav) or \
                    passSingleLeptTrig(event, i2, l2Flav) or \
                    passCrossTrig(event, findSameFlav, muPref, i1, isData, l1Flav) \
                    or passCrossTrig(event, findSameFlav, muPref, i2, isData, l2Flav):
                # found a pair!
                l1Index = i1
                l2Index = i2

            # don't check for opposite charge since needed for regions
    
    #if l1Index>-1 and l2Index > -1 : print '======================', iso1, iso2, pt1,pt2, dr
    if l1Index == -1: return None
    if l2Index == -1: return None


    l3Index = -1
    l4Index = -1
    # search for extra lepton
    l1extra = getattr(event, "nMuon")
    for i3 in range(l1extra):
        if l1Flav == 'Muon' and l2Flav == 'Muon' and (i3 == l1Index or i3==l2Index) : continue
        if l1Flav == 'Muon' and l2Flav == 'Electron' and i3 == l1Index : continue
        if l1Flav == 'Electron' and l2Flav == 'Muon' and i3==l2Index : continue
        if l1Flav == 'Electron' and l2Flav == 'Electron' : continue #this is make sure that we get an extra lepton of the same flavor like the main ones
        if l3Index > -1 : continue 

        pt3 = list(getattr(event, "Muon_pt"))[i3]
        phi3 = list(getattr(event, "Muon_phi"))[i3]
        absDxy3 = abs(list(getattr(event, "Muon_dxy"))[i3])
        absDz3 = abs(list(getattr(event, "Muon_dz"))[i3])
        absEta3 = abs(list(getattr(event, "Muon_eta"))[i3])
        iso3 = list(getattr(event, "Muon_pfRelIso04_all"))[i3]
        isGlobal3 = list(getattr(event, "Muon_isGlobal"))[i3]
        isTracker3 = list(getattr(event, "Muon_isTracker"))[i3]
        isMedium3 = list(getattr(event, "Muon_mediumId"))[i3]

        if  pt3 < 5 or iso3 < maxOkIso or absEta3 > 2.4\
                or absDxy3 > maxOkDxy or absDz3 > maxOkDz : continue
        if  not isMedium3 : continue
        l3Index = i3
        #print 'for ', l1Flav, l2Flav, pt3, 'iso3=========', iso3, 'dxy========', absDxy3, 'dz=========', absDz3, 'l3Ind=======', l3Index, isMedium3

    l2extra = getattr(event, "nElectron")
    for i4 in range(l2extra):
        if l1Flav == 'Electron' and l2Flav == 'Electron' and (i4 == l1Index or i4==l2Index) : continue
        if l1Flav == 'Muon' and l2Flav == 'Electron' and i4 == l2Index : continue
        if l1Flav == 'Electron' and l2Flav == 'Muon' and i4==l1Index : continue
        if l1Flav == 'Muon' and l2Flav == 'Muon' :  continue
        if l4Index > -1 : continue

        pt4 = list(getattr(event, "Electron_pt"))[i4]
        phi4 = list(getattr(event, "Electron_phi"))[i4]
        absDxy4 = abs(list(getattr(event, "Electron_dxy"))[i4])
        absDz4 = abs(list(getattr(event, "Electron_dz"))[i4])
        absEta4 = abs(list(getattr(event, "Electron_eta"))[i4])
        iso4 = list(getattr(event, "Electron_pfRelIso03_all"))[i4]
        hits4 = ord(list(getattr(event, "Electron_lostHits"))[i4])
        conveto4 = list(getattr(event, "Electron_convVeto"))[i4]
        mvaID4 = list(getattr(event, "Electron_mvaFall17V2noIso_WP90"))[i4]

        if  pt4 < 5 or iso4 < maxOkIso  or absEta4 > 2.1\
                or absDxy4 > maxOkDxy or absDz4 > maxOkDz : continue

        if hits4 > 1 or not conveto4 or not mvaID4 : continue
        l4Index = i4
        #print 'for ', l1Flav, l2Flav, pt4, 'iso4=========', iso4, 'dxy========', absDxy4, 'dz=========', absDz4, 'l4Ind=======', l4Index, mvaID4, conveto4
        
    # print l1Index, l2Index, "n1", l1count, "n2", l2count

    # print
    # print "l1 pt: " + str(list(getattr(event, l1Flav+"_pt"))[l1Index])
    # print "l1 eta: " + str(list(getattr(event, l1Flav+"_eta"))[l1Index])
    # print "l1 relIso: " + str(list(getattr(event, l1Flav+"_relIso"))[l1Index])

    # print "l2 pt: " + str(list(getattr(event, l2Flav+"_pt"))[l2Index])
    # print "l2 eta: " + str(list(getattr(event, l2Flav+"_eta"))[l2Index])
    # print "l2 relIso: " + str(list(getattr(event, l2Flav+"_relIso"))[l2Index])
    #print 'now.....', l1Flav, l2Flav, l1Index, l2Index, l3Index, l4Index
    return [l1Index, l2Index, l3Index, l4Index]

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
    if not event.HLT_IsoMu24 and not event.HLT_IsoMu27: return False
    if list(event.Muon_pt)[muIndex] < 26 and event.HLT_IsoMu24 and not event.HLT_IsoMu27 : return False
    if list(event.Muon_pt)[muIndex] < 29 and not event.HLT_IsoMu24 and  event.HLT_IsoMu27 : return False
    for iObj in range(event.nTrigObj):
        dR = deltaR(event, "Muon", muIndex, "TrigObj", iObj)
        # filter bit 2: iso
        if event.TrigObj_filterBits[iObj] & 2 and dR < 0.5: return True
    return False

# Returns true if the electron at elIndex passes the single electron trigger check.
def passElTrigger(event, elIndex):
    if not event.HLT_Ele25_eta2p1_WPTight_Gsf and \
            not event.HLT_Ele27_eta2p1_WPTight_Gsf: return False
    if list(event.Electron_pt)[elIndex] < 27 and event.HLT_Ele25_eta2p1_WPTight_Gsf and not event.HLT_Ele27_eta2p1_WPTight_Gsf: return False
    if list(event.Electron_pt)[elIndex] < 29 and not event.HLT_Ele25_eta2p1_WPTight_Gsf and event.HLT_Ele27_eta2p1_WPTight_Gsf: return False
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

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region A (same sign, nominal rel iso)
def isRegionA(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge > 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region B (opposite sign, nominal rel iso - signal region)
def isRegionB(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge < 0 and l1RelIso < maxRelIso and l2RelIso < maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region C (opposite sign, inverted rel iso)
def isRegionC(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge < 0 and l1RelIso > maxRelIso and l2RelIso > maxRelIso \
            and l1RelIso < 2*maxRelIso and l2RelIso < 2*maxRelIso:
        return True
    return False

# Returns true if an event with l1/l2Charge, l1/l2RelIso, and findSameFlav
# parameters falls in region D (same sign, inverted rel iso)
def isRegionD(l1Charge, l2Charge, l1RelIso, l2RelIso, findSameFlav):
    if findSameFlav: maxRelIso = 0.1
    else: maxRelIso = 0.2
    if l1Charge*l2Charge > 0 and l1RelIso > maxRelIso and l2RelIso > maxRelIso \
            and l1RelIso < 2*maxRelIso and l2RelIso < 2*maxRelIso:
        return True
    return False

#--------------------------------------------------------------------------------#
