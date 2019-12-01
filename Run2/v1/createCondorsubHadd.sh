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

# location of the github repo (base directory)
baseDir="/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY"
# location where condor logs will be placed
condorLogDir="/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs"

# cat > testcondorsub <<EOF
cat > condorsub_haddSubprocess <<EOF
executable  = haddSubprocess.sh
arguments   = $1 $2 $3 $4 $5

universe    = vanilla
initialdir  = $baseDir/Run2/v1
output      = $condorLogDir/haddSubprocess.\$(ClusterId).out
error       = $condorLogDir/haddSubprocess.\$(ClusterId).err
log         = $condorLogDir/haddSubprocess.\$(ClusterId).log

+MaxRuntime = 1140

transfer_input_files = checkRootFile.py

queue
EOF

echo "haddSubprocess.sh" $1 $2 $3 $4 $5
