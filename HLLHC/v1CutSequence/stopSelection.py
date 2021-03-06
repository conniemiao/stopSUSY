# Defines functions to perform object and event selection.

import numpy as np
from math import sqrt, cos

#--------------------------------------------------------------------------------#

# returns delta R for an event between 2 leptons or a lepton and a jet, given a 
# 2 types (type1, type2) and the indices of the 2 objects within the list of all
# objets of that type (index1, index2). types allowed are "muon", "electron", or
# "pfjet"
def deltaR(event, type1, index1, type2, index2):
    assert type1 == "muon" or type1 == "electron" or type1 == "pfjet"
    assert type2 == "muon" or type2 == "electron" or type2 == "pfjet"
    
    eta1 = list(getattr(event, type1+"_eta"))[index1]
    eta2 = list(getattr(event, type2+"_eta"))[index2]
    phi1 = list(getattr(event, type1+"_phi"))[index1]
    phi2 = list(getattr(event, type2+"_phi"))[index2]
    
    return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2))
    
#--------------------------------------------------------------------------------#


# selects a pair of leptons from an event. 
# if findingSameFlav is True: selects for pair of muons if muPreference is True, 
# or for pair of electrons if False. 
# if findingSameFlav is False: selects for leading mu and trailing el if
# muPreference is True, or for leading el and trailing mu if False.
# returns an array of 2 entries where a[0] = l1Index, a[1] = l2Index, where 
# l1 is always the leading lepton (satisfies the higher pt cut), l2 is trailing.
def selectLepts(event, findingSameFlav, muPreference):
    if findingSameFlav:
        if muPreference: # mumu
            l1Flav = "muon"
            l2Flav = "muon"
            maxL1OkEta = 2.4
            maxL2OkEta = 2.4
        else: # elel
            l1Flav = "electron"
            l2Flav = "electron"
            maxL1OkEta = 1.6
            maxL2OkEta = 1.6
        l1MinOkPt = 30
        l2MinOkPt = -0.01 
        maxOkIso = 0.1
    else:
        if muPreference: # muel
            l1Flav = "muon"
            l2Flav = "electron"
            l1MinOkPt = 12
            l2MinOkPt = 15
            maxL1OkEta = 2.4
            maxL2OkEta = 1.6 
        else: # elmu
            l1Flav = "electron"
            l2Flav = "muon"
            l1MinOkPt = 25
            l2MinOkPt = 5
            maxL1OkEta = 1.6
            maxL2OkEta = 2.4
        maxOkIso = 0.2

    # select leading lepton
    l1Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l1count = getattr(event, l1Flav+"_count")
    for i1 in range(l1count):
        pt = list(getattr(event, l1Flav+"_pt"))[i1]
        iso = list(getattr(event, l1Flav+"_relIso"))[i1]
        absEta = abs(list(getattr(event, l1Flav+"_eta"))[i1])
        if pt > l1MinOkPt and iso < maxOkIso and absEta < maxL1OkEta \
                and ((iso < minFoundIso) or (iso == minFoundIso \
                and pt > maxFoundPt)):
            minFoundIso = iso
            maxFoundPt = pt
            l1Index = i1
    if l1Index == -1: return None

    l1Charge = list(getattr(event, l1Flav+"_charge"))[l1Index]

    # select trailing lepton
    l2Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l2count = getattr(event, l2Flav+"_count")
    for i2 in range(l2count):
        iso = list(getattr(event, l2Flav+"_relIso"))[i2]
        absEta = abs(list(getattr(event, l2Flav+"_eta"))[i2])
        l2Charge = list(getattr(event, l2Flav+"_charge"))[i2]
        if l1Charge*l2Charge < 0 and pt > l2MinOkPt and absEta < maxL2OkEta \
                and iso < maxOkIso:
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

# loops over all jets for this event and returns an array of all the indices
# that contain valid jets
def findValidJets(event, l1Flav, l1Index, l2Flav, l2Index):
    jets = []
    numJets = event.pfjet_count
    for j in range(numJets):
        pt = list(event.pfjet_pt)[j]
        # pt0 = list(event.pfjet_pt)[0]
        # if pt != pt0: 
        #     print "pt",pt,"pt0",pt0,
        # print "pt", pt, "eta", list(event.pfjet_eta)[i],
        dRl1 = deltaR(event, l1Flav, l1Index, "pfjet", j)
        dRl2 = deltaR(event, l2Flav, l2Index, "pfjet", j)
        if pt > 20 and abs(list(event.pfjet_eta)[j]) < 2.4\
                and dRl1 > 0.5 and dRl2 > 0.5:
            jets.append(j)
    return jets

#--------------------------------------------------------------------------------#
    
# performs btag cuts on an event, given the jets array containing the indices of
# the accepted jets. returns number of medium cut btags.
def getNumBtag(event, jets):
    # bool pfjet_btag[jet][0, 1, 2]: 0, 1, 2 = passed loose, medium, tight cuts
    # stored as float so > 0.5 = True
    pfjet_btag = np.reshape(event.pfjet_btag, (event.pfjet_count,10))
    numBTagLoose = 0
    numBTag = 0
    numBTagTight = 0
    for jIndex in jets:
        if pfjet_btag[jIndex, 0] > 0.5:
            numBTagLoose += 1
        if pfjet_btag[jIndex, 1] > 0.5:
            numBTag += 1
        if pfjet_btag[jIndex, 2] > 0.5:
            numBTagTight += 1
    return numBTag

#--------------------------------------------------------------------------------#

