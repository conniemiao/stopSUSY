# Creates a condorsub file for each process and each dilepton channel and submits
# the condorsub file to condor.

declare -a processes=("TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
# declare -a processes=("TT+X" "Diboson")

for process in "${processes[@]}"
do
    ./createCondorSub.sh 0 0 0 "$process"
    condor_submit condorsub_makeNtupleBkgd
    ./createCondorSub.sh 0 1 0 "$process"
    condor_submit condorsub_makeNtupleBkgd
    ./createCondorSub.sh 0 1 1 "$process"
    condor_submit condorsub_makeNtupleBkgd
done
