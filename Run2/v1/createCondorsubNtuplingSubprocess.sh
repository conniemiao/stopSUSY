# Creates one condorsub file which will submit makeNtuple.py for all ntuples in the 
# specified subprocess

testMode=$1
inputType=$2
channel=$3
subprocess=$4
process=$5

fileRedirector="$inputType"_fileRedirector
ntupleListsDir="$inputType"NtupleLists/$process

numNtuples=10000 # number of ntuples to loop on for each dataset
if [[ "$testMode" == "test" ]]; then 
    numNtuples=2 
fi

cat > condorsub_makeNtuple <<EOF
executable  = makeNtuple.py

universe    = vanilla
initialdir  = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1
output      = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).out
error       = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).err
log         = /afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/condorLogs/makeNtuple.\$(ClusterId).log

+MaxRuntime = 100000

x509userproxy = /afs/cern.ch/user/c/cmiao/x509up_u112655

transfer_input_files = stopSelection.py, jsonChecker.py, $fileRedirector, $ntupleListsDir

Queue arguments from (
EOF


ntuplesListFile=/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1/$inputType"NtupleLists/"$process/$subprocess
fileNum=0
while read -r ntupleFileName
do
    if [[ $((fileNum + 1)) -gt $numNtuples ]]; then
        break 
    fi

    # appending to condor sub
    echo $testMode $inputType $channel $ntupleFileName $subprocess $process | \
        cat - >> condorsub_makeNtuple
    # ./makeNtuple.py $testMode $inputType $channel $ntupleFileName $subprocess $process
    echo "makeNtuple.py $testMode $inputType $channel $ntupleFileName $subprocess $process"
    let "fileNum++"
done < $ntuplesListFile

# close the parentheses 
echo \) | cat - >> condorsub_makeNtuple
