import ROOT, array

def mt2davis(pt1, eta1, phi1, pt2, eta2, phi2, met, metphi):
    ## NOTE THAT THIS FUNCTION ASSUMES MASSLESS OBJECTS. NOT ADVISED TO USE WITH HEMISPHERES ETC.
    p1 = ROOT.TLorentzVector()
    p2 = ROOT.TLorentzVector()
    mv = ROOT.TLorentzVector()
    p1.SetPtEtaPhiM(pt1,eta1,phi1,0.);
    p2.SetPtEtaPhiM(pt2,eta2,phi2,0.);
    mv.SetPtEtaPhiM(met,0.,metphi,0.);
    a = array.array('d', [p1.M(), p1.Px(), p1.Py()])
    b = array.array('d', [p2.M(), p2.Px(), p2.Py()])
    c = array.array('d', [mv.M(), mv.Px(), mv.Py()])
                                                                                                                                   
    mt2obj = ROOT.Davismt2()
    mt2obj.set_momenta( a, b, c );
    mt2obj.set_mn( 0. );

    result = mt2obj.get_mt2();
    return result

def mt2(event, l1Flav, l1Index, l2Flav, l2Index):
    pt1 = list(getattr(event, l1Flav + "_pt"))[l1Index]
    pt2 = list(getattr(event, l2Flav + "_pt"))[l2Index]
    eta1 = list(getattr(event, l1Flav + "_eta"))[l1Index]
    eta2 = list(getattr(event, l2Flav + "_eta"))[l2Index]
    phi1 = list(getattr(event, l1Flav + "_phi"))[l1Index]
    phi2 = list(getattr(event, l2Flav + "_phi"))[l2Index]
    met = event.MET_pt
    metphi = event.MET_phi
    return mt2davis(pt1, eta1, phi1, pt2, eta2, phi2, met, metphi)

if __name__ == "__main__":
    ROOT.gROOT.ProcessLine(".include mt2/")
    ROOT.gROOT.ProcessLine(".L Davismt2.cc+" )
    print mt2davis(35., -1.2, 1.4, 43.2, 0.8, 0.8, 65., 0.3) 
