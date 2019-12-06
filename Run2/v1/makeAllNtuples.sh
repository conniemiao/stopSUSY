# NOTE: requires 1 argument, {test, all} to indicate whether submitting jobs in 
# testing mode
#
# By calling makeAllNtuplesProcess.sh, creates a condorsub file for each data, bkgd, 
# and sig process and for each dilepton channel and submits the condorsub file to 
# condor.
#
# Args to makeNtuple.py: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# ntuple, subprocess, process (only required for bkgd)

testMode=$1
if [[ "$testMode" == "test" ]]; then 
    channels=("elel")
    bkgdProcesses=("QCD")
elif [[ "$testMode" == "all" ]]; then
    # channels=("muel" "elel" "mumu")
    channels=("mumu")
    bkgdProcesses=("W-Jets" "Drell-Yan" "Diboson" "Single-Top" "TTBar" "TT+X" "QCD")
    # bkgdProcesses=("QCD")
else
    echo "need {test, all} as arg to makeAllNtuples.sh"
    exit 1
fi


for channel in "${channels[@]}"
do
    echo ------------------ "$channel" ------------------
    # Bkgd
    echo --- bkgd ---
    for process in "${bkgdProcesses[@]}"
    do
        bash makeAllNtuplesProcess.sh $testMode bkgd "$channel" "$process"
        echo
    done

    # Sig
    echo --- sig ---
    bash makeAllNtuplesProcess.sh $testMode sig "$channel"
    echo

    # Data 
    echo --- data ---
    bash makeAllNtuplesProcess.sh $testMode data "$channel"
    
done
