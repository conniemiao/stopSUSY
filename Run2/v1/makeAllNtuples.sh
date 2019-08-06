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

# declare -a channels=("mumu" "muel" "elel")
declare -a channels=("elel")
# declare -a bkgdProcesses=("TTBar" "TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
declare -a bkgdProcesses=("Single-Top")

for channel in "${channels[@]}"
do
    echo ------------------ "$channel" ------------------
    # Bkgd
    echo --- bkgd ---
    for process in "${bkgdProcesses[@]}"
    do
        bash makeAllNtuplesProcess.sh $1 bkgd "$channel" "$process"
        echo
    done

    # Sig
    echo --- sig ---
    bash makeAllNtuplesProcess.sh $1 sig "$channel"
    echo

    # Data 
    echo --- data ---
    bash makeAllNtuplesProcess.sh $1 data "$channel"
    
done
