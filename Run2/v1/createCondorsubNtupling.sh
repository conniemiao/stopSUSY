# Creates a condorsub file for makeNtuple.py with $1-4 as the arguments to 
# makeNtuple.py

fileRedirector="$2"_fileRedirector
ntupleListsDir="$2"_NtupleLists/$4

# cat > testcondorsub <<EOF
cat > condorsub_makeNtuple <<EOF
executable  = makeNtuple.py
arguments   = $1 $2 $3 $4

universe    = vanilla
output      = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).out
error       = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).err
log         = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).log

use_x509userproxy = true

+MaxRuntime = 100000

transfer_input_files = stopSelection.py, $fileRedirector, $ntupleListsDir, Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt

queue
EOF

echo "makeNtuple.py $1 $2 $3 $4"
