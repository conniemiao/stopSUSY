# Creates one condorsub file which will submit makeNtuple.py for all ntuples in the 
# specified subprocess

# location of the github repo (base directory)
baseDir="/afs/cern.ch/user/a/alkaloge/work/Connie/CMSSW_10_2_9/src/stopSUSY/"
# location where condor logs will be placed
condorLogDir="/afs/cern.ch/user/a/alkaloge/work/Connie/CMSSW_10_2_9/src/stopSUSY/Run2/v1/Logs/"
# address to the GRID certificate created with voms-proxy-init
x509Adr="/afs/cern.ch/user/a/alkaloge/x509up_u13585"

testMode=$1
inputType=$2
channel=$3
subprocess=$4
process=$5

numNtuples=10000 # number of ntuples to loop on for each dataset
if [[ "$testMode" == "test" ]]; then 
    numNtuples=2 
fi

cat > condorsub_makeNtuple <<EOF
executable  = makeNtuple.py

universe    = vanilla
initialdir  = $baseDir/Run2/v1
output      = $condorLogDir/${channel}_${ntupleFileName}.\$(ClusterId).\$(ProcId).out
error       = $condorLogDir/${channel}_${ntupleFileName}.\$(ClusterId).\$(ProcId).err
log         = $condorLogDir/${channel}_${ntupleFileName}.\$(ClusterId).\$(ProcId).log

+MaxRuntime = 129600

x509userproxy = $x509Adr

transfer_input_files = stopSelection.py, jsonChecker.py

Queue arguments from (
EOF

ntuplesListFile=$baseDir/Run2/v1/$inputType"NtupleLists/"$process/$subprocess

fileNum=0
while read -r ntupleFileName
do
    if [[ $((fileNum + 1)) -gt $numNtuples ]]; then
        break 
    fi

    if [[ "$ntupleFileName" =~ \#.* ]]; then 
        continue
    fi

    # appending to condor sub
    echo $testMode $inputType $channel $ntupleFileName $subprocess $process | \
        cat - >> condorsub_makeNtuple
    if [[ "$testMode" == "test" ]]; then 
        ./makeNtuple.py $testMode $inputType $channel $ntupleFileName $subprocess \
            $process
    fi
    echo "makeNtuple.py $testMode $inputType $channel $ntupleFileName $subprocess \
        $process"
    let "fileNum++"
done < $ntuplesListFile

# close the parentheses 
echo \) | cat - >> condorsub_makeNtuple
