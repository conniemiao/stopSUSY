# Creates a condorsub file for each data, bkgd, and sig process and for each dilepton
# channel and submits the condorsub file to condor.
#
# Args to makeNtuple.py: 
# testMode {test, real}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# process

declare -a channels=("mumu" "muel" "elel")

# Data
bash createCondorsubNtupling.sh real data mumu DoubleMuon
condor_submit condorsub_makeNtuple
bash createCondorsubNtupling.sh real data muel MuonEG
condor_submit condorsub_makeNtuple
bash createCondorsubNtupling.sh real data elel DoubleElectron
condor_submit condorsub_makeNtuple

declare -a bkgdProcesses=("TTBar" "TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
# declare -a bkgdProcesses=("TTBar")
for channel in "${channels[@]}"
do
    # Bkgd
    for process in "${bkgdProcesses[@]}"
    do
        bash createCondorsubNtupling.sh real bkgd "$channel" "$process"
        condor_submit condorsub_makeNtuple
    done
done
