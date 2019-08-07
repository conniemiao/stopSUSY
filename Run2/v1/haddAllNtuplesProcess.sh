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

# the subfolder arg is not necessary for data or sig, but just put it there to make
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
while read -r datasetName subfolder xsec
do
    if [[ $datasetNum+1 > $numDatasets ]]; then
        break 
    fi
    if [[ "$datasetName" =~ \#.* ]]; then 
        continue
    fi
    if [[ "$subfolder" != "$expectedSubfolder" ]]; then
        continue
    fi

    bash createCondorsubHadd.sh $testMode $inputType $channel $datasetName $subfolder
    # condor_submit condorsub_haddSubprocess
    ./haddSubprocess.sh $testMode $inputType $channel $datasetName $subfolder

    let "datasetNum++"
    echo
done < $redirectorFileAdr
