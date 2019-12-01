# NOTE: NEEDS 3-4 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# process (only required for bkgd)
#
# Creates an hadd submit file and submits it to condor for every dataset (aka 
# subprocess) in the selected input type, channel, and process. (This is the
# script that actually does the condor submission handling)
#
# Args to haddSubprocess.sh:
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# subprocess, process

testMode=$1
inputType=$2
channel=$3
process=$4

# note: must loop over all files and all datasets to have correct xsec
numNtuples=10000 # number of ntuples to loop on for each dataset
numDatasets=10000 # number of datasets to loop on for the process
if [[ "$testMode" == "test" ]]; then 
    numNtuples=2 
    numDatasets=2
fi

redirectorFileAdr=$inputType"_fileRedirector"

# the process arg is not necessary for data or sig, but just put it there to make
# the code simpler
if [[ "$inputType" == "data" ]]; then
    if [[ "$channel" == "mumu" ]]; then
        expectedSubfolder="DoubleMuon"
    elif [[ "$channel" == "elel" ]]; then
        expectedSubfolder="DoubleEG"
    elif [[ "$channel" == "muel" ]]; then
        expectedSubfolder="MuonEG"
    fi
elif [[ "$inputType" == "bkgd" ]]; then
    expectedSubfolder=$process
elif [[ "$inputType" == "sig" ]]; then
    expectedSubfolder="Stop-Pair"
fi 

# start submitting files
datasetNum=0
while read -r subprocess process xsec
do
    if [[ $((datasetNum + 1)) -gt $numDatasets ]]; then
        break 
    fi
    if [[ "$subprocess" =~ \#.* ]]; then 
        continue
    fi
    if [[ "$process" != "$expectedSubfolder" ]]; then
        continue
    fi

    bash createCondorsubHadd.sh $testMode $inputType $channel $subprocess $process
    #if [[ "$testMode" == "all" ]]; then 
    #    condor_submit condorsub_haddSubprocess
    #else
    #    ./haddSubprocess.sh $testMode $inputType $channel $subprocess $process
    #fi
    ./haddSubprocess.sh $testMode $inputType $channel $subprocess $process

    let "datasetNum++"
    echo
done < $redirectorFileAdr
