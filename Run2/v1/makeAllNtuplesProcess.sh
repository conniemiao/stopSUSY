# NOTE: NEEDS 3-4 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# process (only required for bkgd)
# Process options: 
#   (optional) data: DoubleMuon, DoubleEG, MuonEG
#   bkgd: TTBar TT+X Diboson W-Jets Drell-Yan Single-Top
#   (optional) sig: Stop-Pair
# 
# Creates a submit file and submits the job to condor for every dataset (aka 
# subprocess) in the selected input type, channel, and process. (This is the
# script that actually does the condor submission handling)
#
# Args to makeNtuple.py: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# ntuple, subprocess, process (only required for bkgd)

testMode=$1
# number of ntuples to loop on for each dataset
# note: must loop over all files and all datasets to have correct xsec
numNtuples=10000
if [[ "$testMode" == "test" ]]; then 
    numNtuples=2 
fi

inputType=$2

redirectorFileAdr=$inputType"_fileRedirector"

# the subfolder arg is not necessary for data or sig, but just put it there to make
# the code simpler
if [[ "$inputType" == "data" ]]; then
    if [[ "$3" == "mumu" ]]; then
        expectedSubfolder="DoubleMuon"
    elif [[ "$3" == "elel" ]]; then
        expectedSubfolder="DoubleEG"
    elif [[ "$3" == "muel" ]]; then
        expectedSubfolder="MuonEG"
    fi
elif [[ "$inputType" == "bkgd" ]]; then
    expectedSubfolder=$4
elif [[ "$inputType" == "sig" ]]; then
    expectedSubfolder="Stop-Pair"
fi 

# necessary for condor to work nicely with das
cp /tmp/x509up_u112655 .
chmod 777 x509up_u112655

# start submitting files
while read -r datasetName subfolder xsec
do
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
        bash createCondorsubNtupling.sh $1 $2 $3 $ntupleFileName $datasetName $subfolder
        condor_submit condorsub_makeNtuple
        # ./makeNtuple.py $1 $2 $3 $ntupleFileName $datasetName $subfolder
        
        echo
        let "fileNum++"
    done < $ntuplesListFile
done < $redirectorFileAdr
