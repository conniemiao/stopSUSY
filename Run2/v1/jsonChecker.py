# Alexis's script for checking a json file and determining whether a given event is 
# good.

import json

class jsonChecker():

    def __init__(self, filein='/afs/cern.ch/work/c/cmiao/private/myDataSusy/Run2/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'):
        self.good, self.bad = 0, 0
        input_file = open(filein)
        self.json_array = json.load(input_file)

    # Given the LS (lumisection aka lumiblock) and run number, checks the json file 
    # and returns True if they are good, or False if not.
    def checkJSON(self, LS, run):
        try: 
            LSlist = self.json_array[str(run)]
            for LSrange in LSlist:
                if LS >= LSrange[0] and LS <= LSrange[1]:
                    self.good += 1
                    return True
        except KeyError: pass
        self.bad += 1
        return False

    def printJSONsummary(self):
        print "check JSON summary: nCalls=", self.good+self.bad, "nGood=", self.good, "nBad=", \
                self.bad
        return 
