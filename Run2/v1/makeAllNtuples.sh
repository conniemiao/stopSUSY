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

if [ "$1" != "test" ] && [ "$1" != "all" ]; then
    echo "need {test, all} as arg to makeAllNtuples.sh"
    exit 1
fi

declare -a channels=("mumu" "muel" "elel")

declare -a bkgdProcesses=("TTBar" "TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
# declare -a bkgdProcesses=("TTBar")
for channel in "${channels[@]}"
do
    # Bkgd
    for process in "${bkgdProcesses[@]}"
    do
        bash makeAllNtuplesProcess.sh $1 bkgd "$channel" "$process"
    done

    # Sig
    bash makeAllNtuplesProcess.sh $1 sig "$channel"

    # Data 
    bash makeAllNtuplesProcess.sh $1 data "$channel"
    
done
