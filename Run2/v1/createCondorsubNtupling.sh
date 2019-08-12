# Creates a condorsub file for makeNtuple.py with $1-6 as the arguments to 
# makeNtuple.py

fileRedirector="$2"_fileRedirector
ntupleListsDir="$2"NtupleLists/$6

# cat > testcondorsub <<EOF
cat > condorsub_makeNtuple <<EOF
executable  = makeNtuple.py
arguments   = $1 $2 $3 $4 $5 $6

universe    = vanilla
initialdir  = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1
output      = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).out
error       = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).err
log         = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).log

+MaxRuntime = 100000

x509userproxy = /afs/cern.ch/user/c/cmiao/x509up_u112655

transfer_input_files = stopSelection.py, jsonChecker.py, $fileRedirector, $ntupleListsDir

queue
EOF

echo "makeNtuple.py $1 $2 $3 $4 $5 $6"
