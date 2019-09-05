# Creates a condorsub file for one of the plotting scripts (plot1D.py, 
# generateCutflows.py, plot2DMetPtGraphs.py) with $@[2:end] containing all the 
# arguments and $1 containing the name of the script.

# location of the github repo (base directory)
baseDir="/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY"
# location where condor logs will be placed
condorLogDir="/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs"

# cat > testcondorsub <<EOF
cat > condorsub_plotting <<EOF
executable  = $1
arguments   = ${@:2}

universe    = vanilla
initialdir  = $baseDir/Run2/v1
output      = $condorLogDir/${1%.py}.\$(ClusterId).out
error       = $condorLogDir/${1%.py}.\$(ClusterId).err
log         = $condorLogDir/${1%.py}.\$(ClusterId).log

+MaxRuntime = 129500

transfer_input_files = stopSelection.py, sig_fileRedirector, bkgd_fileRedirector, data_fileRedirector

queue
EOF

echo "$@"
