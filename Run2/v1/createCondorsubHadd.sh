# NOTE: NEEDS 5 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# subprocess, process
# Process options: 
#   data: DoubleMuon, DoubleEG, MuonEG
#   bkgd: TTBar TT+X Diboson W-Jets Drell-Yan Single-Top
#   sig: Stop-Pair
# Subprocess: name of the dataset for data and sig, or name of the subprocess for bkgd
#
# Creates a condorsub file for haddSubprocess.sh with $1-5 as the arguments to 
# haddSubprocess.sh

# cat > testcondorsub <<EOF
cat > condorsub_haddSubprocess <<EOF
executable  = haddSubprocess.sh
arguments   = $1 $2 $3 $4 $5

universe    = vanilla
initialdir  = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1
output      = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/haddSubprocess.\$(ClusterId).out
error       = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/haddSubprocess.\$(ClusterId).err
log         = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/haddSubprocess.\$(ClusterId).log

+MaxRuntime = 100000

queue
EOF

echo "haddSubprocess.sh" $1 $2 $3 $4 $5
