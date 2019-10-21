# Checks if some ntuple file is a zombie and deletes it if it is.

import sys, os
from ROOT import TFile

if __name__ == "__main__":
    assert len(sys.argv) == 2, "need 1 command line arg: fileadr to check"

    fileadr = sys.argv[1]
    if not os.path.exists(fileadr):
        print fileadr, "does not exist."
    else:
        checkFile = TFile(fileadr)
        if checkFile.IsZombie():
            print "Zombie file, removing", fileadr
            os.remove(fileadr)
