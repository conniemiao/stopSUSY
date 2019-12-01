import ROOT, array

ROOT.gROOT.ProcessLine(".include ./")
ROOT.gROOT.ProcessLine(".L Davismt2.cc+" )


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


print mt2davis(35., -1.2, 1.4, 43.2, 0.8, 0.8, 65., 0.3) 


