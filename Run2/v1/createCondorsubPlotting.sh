# Creates a condorsub file for one of the plotting scripts (plot1D.py, 
# generateCutflows.py, plot2DMetPtGraphs.py) with $@[2:end] containing all the 
# arguments and $1 containing the name of the script.

# cat > testcondorsub <<EOF
cat > condorsub_plotting <<EOF
executable  = $1
arguments   = ${@:2}

universe    = vanilla
output      = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/${1%.py}.\$(ClusterId).out
error       = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/${1%.py}.\$(ClusterId).err
log         = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/${1%.py}.\$(ClusterId).log

+MaxRuntime = 129500

transfer_input_files = stopSelection.py, sig_fileRedirector, bkgd_fileRedirector, data_fileRedirector

queue
EOF

echo "$@"
