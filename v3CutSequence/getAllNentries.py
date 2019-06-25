# outputs all the nentries for each of the bkgd and sig input files into new
# txt files
from ROOT import TFile, TTree

bkgdDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/DESY_pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/"
bkgdDataListFile = open("bkgd_TTDiLept_files")
bkgdNentriesOutfile = open("bkgd_TTDiLept_files_nentries", "w+")
for fileNum, line in enumerate(bkgdDataListFile):
    filename = line.rstrip()
    print filename
    inFile = TFile.Open(bkgdDataDir + filename, "READ")
    inTree = inFile.Get("AC1B")
    bkgdNentriesOutfile.write("%s %i\n" %(filename, inTree.GetEntries()))

sigDataDir = "/eos/user/a/alkaloge/HLLHC/Skims/v3/SingleStop/"
sigDataListFile = open("sig_SingleStop_files")
sigNentriesOutfile = open("sig_SingleStop_files_nentries", "w+")
for fileNum, line in enumerate(sigDataListFile):
    line = line.rstrip('\n')
    filename, xsec = line.split(" ")
    print filename
    inFile = TFile.Open(sigDataDir + filename, "READ")
    inTree = inFile.Get("AC1B")
    sigNentriesOutfile.write("%s %s %i\n" %(filename, xsec, inTree.GetEntries()))
