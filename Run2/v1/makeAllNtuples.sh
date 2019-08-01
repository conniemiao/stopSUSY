# NOTE: requires 1 argument, {test, all} to indicate whether submitting jobs in 
# testing mode
# Creates a condorsub file for each data, bkgd, and sig process and for each dilepton
# channel and submits the condorsub file to condor.
#
# Args to makeNtuple.py: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# process

if [ "$1" != "test" ] && [ "$1" != "all" ]; then
    echo "need {test, all} as arg to makeAllNtuples.sh"
    exit 1
fi

declare -a channels=("mumu" "muel" "elel")

# Data
bash createCondorsubNtupling.sh $1 data mumu DoubleMuon
condor_submit condorsub_makeNtuple
bash createCondorsubNtupling.sh $1 data muel MuonEG
condor_submit condorsub_makeNtuple
bash createCondorsubNtupling.sh $1 data elel DoubleElectron
condor_submit condorsub_makeNtuple

declare -a bkgdProcesses=("TTBar" "TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
# declare -a bkgdProcesses=("TTBar")
for channel in "${channels[@]}"
do
    # Bkgd
    for process in "${bkgdProcesses[@]}"
    do
        bash createCondorsubNtupling.sh $1 bkgd "$channel" "$process"
        condor_submit condorsub_makeNtuple
    done

    # Sig
    bash createCondorsubNtupling.sh $1 sig "$channel" "Stop-Pair"
    condor_submit condorsub_makeNtuple
    
done
