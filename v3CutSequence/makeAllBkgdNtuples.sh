declare -a processes=("TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")

for process in "${processes[@]}"
do
    ./createCondorSub.sh 0 0 0 "$process"
    condor_submit testcondorsub
    ./createCondorSub.sh 0 1 0 "$process"
    condor_submit testcondorsub
    ./createCondorSub.sh 0 1 1 "$process"
    condor_submit testcondorsub
done
