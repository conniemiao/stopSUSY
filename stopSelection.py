# defines functions to perform object and event selection.

import numpy as np

# selects a pair of same flavor leptons from an event. selects for muons if 
# findingMuons is True, or for electrons if False. returns an array of 
# 2 entries where a[0] = l1Index, a[1] = l2Index
def selectSameFlavLepts(event, findingMuons):
    if findingMuons:
        flav = "muon"
        maxOkEta = 2.4
    else:
        flav = "electron"
        maxOkEta = 1.6
    l1MinOkPt = 30
    maxOkIso = 0.1

    l1Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l1count = getattr(event, flav+"_count")
    for i1 in range(l1count):
        pt = list(getattr(event, flav+"_pt"))[i1]
        iso = list(getattr(event, flav+"_relIso"))[i1]
        absEta = abs(list(getattr(event, flav+"_eta"))[i1])
        if pt > l1MinOkPt and iso < maxOkIso and absEta < maxOkEta \
        and ((iso < minFoundIso) or (iso == minFoundIso and pt > maxFoundPt)):
            minFoundIso = iso
            maxFoundPt = pt
            l1Index = i1
    if l1Index == -1: return None

    l1Charge = list(getattr(event, flav+"_charge"))[l1Index]

    l2Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l2count = getattr(event, flav+"_count")
    for i2 in range(l2count):
        iso = list(getattr(event, flav+"_relIso"))[i2]
        absEta = abs(list(getattr(event, flav+"_eta"))[i2])
        l2Charge = list(getattr(event, flav+"_charge"))[i2]
        if l1Charge*l2Charge < 0 and absEta < maxOkEta and iso < maxOkIso:
            pt = list(getattr(event, flav+"_pt"))[i2]
            if (iso < minFoundIso) or (iso == minFoundIso and pt > maxFoundPt):
                minFoundIso = iso
                maxFoundPt = pt
                l2Index = i2
    if l2Index == -1: return None

    return [l1Index, l2Index]

# loops over all jets for this event and returns an array of all the indices
# that contain valid jets
def findValidJets(event):
    jets = []
    numJets = event.pfjet_count
    for i in range(numJets):
        pt = list(event.pfjet_pt)[i]
        if pt > 20 and abs(list(event.pfjet_eta)[i]) < 2.4:
            jets.append(i)
    return jets

    
# performs btag cuts on an event. returns True if passes, False if not.
def passesCut(event):
    # if event.pfjet_count > 4: return False

    # bool pfjet_btag[jet][0, 1, 2]: 0, 1, 2 = passed loose, medium, tight cuts
    # stored as float so > 0.5 = True
    pfjet_btag = np.reshape(event.pfjet_btag, (event.pfjet_count,10))
    numBTagLoose = 0
    numBTag = 0
    numBTagTight = 0
    for jet in range(event.pfjet_count):
        if pfjet_btag[jet, 0] > 0.5:
            numBTagLoose += 1
        if pfjet_btag[jet, 1] > 0.5:
            numBTag += 1
        if pfjet_btag[jet, 2] > 0.5:
            numBTagTight += 1
    if numBTag > 1: return False
    return True

