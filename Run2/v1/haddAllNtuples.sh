# NOTE: requires 1 argument, {test, all} to indicate whether submitting jobs in 
# testing mode
#
# By calling haddAllNtuplesProcess.sh, creates a condorsub file for each data, bkgd, 
# and sig process and for each dilepton channel and submits the condorsub file to 
# condor.
#
# Args to haddSubprocess.sh:
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# subprocess, process

testMode=$1
if [ "$testMode" != "test" ] && [ "$testMode" != "all" ]; then
    echo "need {test, all} as arg to haddAllNtuples.sh"
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
        bash haddAllNtuplesProcess.sh $testMode bkgd "$channel" "$process"
        echo
    done

    # Sig
    echo --- sig ---
    bash haddAllNtuplesProcess.sh $testMode sig "$channel"
    echo

    # Data 
    echo --- data ---
    bash haddAllNtuplesProcess.sh $testMode data "$channel"
    
done
