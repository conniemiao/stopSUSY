# see DesyTauAnalyses/NTupleMaker/interface/Jets.h

def looseJetiD(entry, jet):
    looseJetID = False
    energy = entry.pfjet_e[jet]
    energy *= entry.pfjet_energycorr[jet] # uncorrected energy must be used
    eta = entry.pfjet_eta[jet]
    chf = entry.pfjet_chargedhadronicenergy[jet]/energy
    nhf = entry.pfjet_neutralhadronicenergy[jet]/energy
    nem = entry.pfjet_neutralemenergy[jet]/energy
    elf = entry.pfjet_chargedemenergy[jet]/energy
    muf = entry.pfjet_muonenergy[jet]/energy
    chm = entry.pfjet_chargedmulti[jet] 
    nm  = entry.pfjet_neutralmulti[jet]
    npr = entry.pfjet_chargedmulti[jet] + entry.pfjet_neutralmulti[jet]

    if abs(eta) <= 2.7: looseJetID = (nhf < 0.99 and nem < 0.99 and npr > 1) and \
            (abs(eta) > 2.4 or (chf > 0 and chm > 0 and elf < 0.99))
    elif abs(eta) <= 3.0: looseJetID = nem < 0.9 and nm > 2
    else: looseJetID = nem < 0.9 and nm > 10
    return looseJetID
