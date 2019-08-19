# NOTE: requires 1 argument, {test, all} to indicate whether submitting jobs in 
# testing mode
#
# By calling haddAllNtuplesProcess.sh, creates a condorsub file for each data, bkgd, 
# and sig subprocess and for each dilepton channel and submits the condorsub file to 
# condor.
#
# Args to haddSubprocess.sh:
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# subprocess, process

testMode=$1
if [[ "$testMode" == "test" ]]; then 
    channels=("elel")
    bkgdProcesses=("TTBar" "W-Jets")
elif [[ "$testMode" == "all" ]]; then
    # channels=("mumu" "muel" "elel")
    channels=("elel")
    # bkgdProcesses=("W-Jets" "Drell-Yan" "Diboson" "Single-Top" "TTBar" "TT+X" "QCD")
    bkgdProcesses=("W-Jets" "Drell-Yan" "Diboson" "Single-Top" "TTBar" "TT+X")
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
