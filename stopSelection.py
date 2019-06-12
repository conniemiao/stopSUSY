# defines functions to perform object and event selection.

# selects a pair of same-flavor leptons from an entry. selects for muons if
# findingMuons is True, or for electrons if False. returns an array of 3 entries
# where a[0] = l
def selectSameFlaveLepts(entry, findingMuons)

# performs btag cuts on an event. returns True if passes, False if not.
def passesCut(entry):
    # if entry.pfjet_count > 4: return False

    # bool pfjet_btag[jet][0, 1, 2]: 0, 1, 2 = passed loose, medium, tight cuts
    # stored as float so > 0.5 = True
    pfjet_btag = np.ndarray((entry.pfjet_count, 2), 'f', entry.pfjet_btag)
    numBTagLoose = 0
    numBTag = 0
    numBTagTight = 0
    for jet in range(entry.pfjet_count):
        if pfjet_btag.item((jet, 0)) > 0.5:
            numBTagLoose += 1
        if pfjet_btag.item((jet, 1)) > 0.5:
            numBTag += 1
        # if pfjet_btag.item((jet, 2)) > 0.5:
        #     numBTagTight += 1
    if numBTag > 1: return False
    return True
