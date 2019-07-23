#!/usr/bin/env python

# Usage: 1 arg: directory of the root files to check.  Checks all root files in 
# the directory to see if they contain the tree AC1B (adapted from Alexis's code)

from ROOT import TFile, TTree
import sys
import glob

def checkDir(directory):
    for f in glob.glob(directory+"/*.root"):
        print f 
        f = TFile(f)
    
        try:
            t = f.Get("AC1B")
            #print 'the filename is ',file, myTree.GetEntries()
        except AttributeError:
            print "Problematic file:", f 

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Requires 1 argument, the dir to check"
    assert sys.argv[1][-1:] != "/", "Do not include the / at the end of the dir"
    checkDir(sys.argv[1])
