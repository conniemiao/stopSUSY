# NOTE: NEEDS 3-4 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# process (only required for bkgd)
# Process options: 
#   (optional) data: DoubleMuon, DoubleEG, MuonEG
#   bkgd: TTBar TT+X Diboson W-Jets Drell-Yan Single-Top
#   (optional) sig: Stop-Pair
# 
# Creates a makeNtuple submit file and submits it to condor for every dataset (aka 
# subprocess) in the selected input type, channel, and process. (This is the
# script that actually does the condor submission handling)
#
# Args to makeNtuple.py: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# ntuple, subprocess, process (only required for bkgd)

testMode=$1
inputType=$2
channel=$3
process=$4 # only required for bkgd

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

# necessary for condor to work nicely with das
cp /tmp/x509up_u112655 /afs/cern.ch/user/c/cmiao
chmod 777 /afs/cern.ch/user/c/cmiao/x509up_u112655
export X509_USER_PROXY=/afs/cern.ch/user/c/cmiao/x509up_u112655

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
    echo Submitting $datasetName ntuples from $subfolder
    ntuplesListFile=/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/stopSUSY/Run2/v1/$inputType"NtupleLists/"$subfolder/$datasetName

    fileNum=0
    while read -r ntupleFileName
    do
        if [[ $fileNum+1 > $numNtuples ]]; then
            break 
        fi
        bash createCondorsubNtupling.sh $testMode $inputType $channel $ntupleFileName\
            $datasetName $subfolder
        condor_submit condorsub_makeNtuple
        # ./makeNtuple.py $testMode $inputType $channel $ntupleFileName $datasetName $subfolder
        
        let "fileNum++"
    done < $ntuplesListFile
    let "datasetNum++"
    echo
done < $redirectorFileAdr
