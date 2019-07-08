# Creates a makeNtupleBkgd condorsub file for each bkgd process and each dilepton
# channel and submits the condorsub file to condor.
# Args to makeNtupleBkgd: test mode, same flav, mu pref, bkgd processes to stack

# declare -a processes=("TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
declare -a processes=("TT+X" "W-Jets" "Drell-Yan" "Single-Top")

for process in "${processes[@]}"
do
    ./createCondorsubBkgdNtuple.sh 0 0 0 "$process"
    condor_submit condorsub_makeNtupleBkgd
    ./createCondorsubBkgdNtuple.sh 0 1 0 "$process"
    condor_submit condorsub_makeNtupleBkgd
    ./createCondorsubBkgdNtuple.sh 0 1 1 "$process"
    condor_submit condorsub_makeNtupleBkgd
done
