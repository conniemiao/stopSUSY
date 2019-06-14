# defines functions to perform object and event selection.

import numpy as np

# selects a pair of leptons from an event. 
# if findingSameFlav is True: selects for pair of muons if muPreference is True, 
# or for pair of electrons if False. 
# if findingSameFlav is False: selects for leading mu and trailing el if
# muPreference is True, or for leading el and trailing mu if False.
# returns an array of 2 entries where a[0] = l1Index, a[1] = l2Index, where in the
# case of finding opposite flavors, l1 = leading, l2 = trailing.
def selectLepts(event, findingSameFlav, muPreference):
    if findingSameFlav:
        if muPreference: # mu-mu
            flav1 = "muon"
            flav2 = "muon"
            maxL1OkEta = 2.4
            maxL2OkEta = 2.4
        else: # el-el
            flav2 = "electron"
            flav2 = "electron"
            maxL2OkEta = 1.6
            maxL2OkEta = 1.6
        l1MinOkPt = 30
        l2MinOkPt = 0
        maxOkIso = 0.1
    else:
        if muPreference: # mu-el
            flav1 = "muon"
            flav2 = "electron"
            maxL1OkEta = 2.4
            maxL2OkEta = 1.6 
        else: # el-mu
            flav1 = "electron"
            flav2 = "muon"
            maxL1OkEta = 1.6
            maxL2OkEta = 2.4
        l1MinOkPt = 25 
        l2MinOkPt = 15
        maxOkIso = 0.2


    l1Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l1count = getattr(event, flav1+"_count")
    for i1 in range(l1count):
        pt = list(getattr(event, flav1+"_pt"))[i1]
        iso = list(getattr(event, flav1+"_relIso"))[i1]
        absEta = abs(list(getattr(event, flav1+"_eta"))[i1])
        if pt > l1MinOkPt and iso < maxOkIso and absEta < maxL1OkEta \
                and ((iso < minFoundIso) or (iso == minFoundIso \
                and pt > maxFoundPt)):
            minFoundIso = iso
            maxFoundPt = pt
            l1Index = i1
    if l1Index == -1: return None

    l1Charge = list(getattr(event, flav1+"_charge"))[l1Index]

    l2Index = -1
    minFoundIso = 1
    maxFoundPt = 0
    l2count = getattr(event, flav2+"_count")
    for i2 in range(l2count):
        iso = list(getattr(event, flav2+"_relIso"))[i2]
        absEta = abs(list(getattr(event, flav2+"_eta"))[i2])
        l2Charge = list(getattr(event, flav2+"_charge"))[i2]
        if l1Charge*l2Charge < 0 and pt > l2MinOkPt and absEta < maxL2OkEta \
                and iso < maxOkIso:
            pt = list(getattr(event, flav2+"_pt"))[i2]
            if (iso < minFoundIso) or (iso == minFoundIso and pt > maxFoundPt):
                minFoundIso = iso
                maxFoundPt = pt
                l2Index = i2
    if l2Index == -1: return None

    print
    print "l1 pt: " + str(list(getattr(event, flav1+"_pt"))[l1Index])
    print "l1 eta: " + str(list(getattr(event, flav1+"_eta"))[l1Index])
    print "l1 relIso: " + str(list(getattr(event, flav1+"_relIso"))[l1Index])

    print "l2 pt: " + str(list(getattr(event, flav2+"_pt"))[l2Index])
    print "l2 eta: " + str(list(getattr(event, flav2+"_eta"))[l2Index])
    print "l2 relIso: " + str(list(getattr(event, flav2+"_relIso"))[l2Index])
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

