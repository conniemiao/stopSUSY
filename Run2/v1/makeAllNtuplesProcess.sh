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
numDatasets=10000 # number of datasets to loop on for the process
if [[ "$testMode" == "test" ]]; then 
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

# necessary for condor to work nicely with das
cp /tmp/x509up_u112655 /afs/cern.ch/user/c/cmiao
if [ $? -ne 0 ]; then
    echo "You need to do voms-proxy init!"
    exit 1
fi
chmod 777 /afs/cern.ch/user/c/cmiao/x509up_u112655

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

    bash createCondorsubNtuplingSubprocess.sh $testMode $inputType $channel \
        $subprocess $process
    if [[ "$testMode" == "all" ]]; then 
        echo Submitting $subprocess ntuples from $process
        condor_submit condorsub_makeNtuple
    else
        echo Running on $subprocess ntuples from $process
    fi

    let "datasetNum++"
    echo
done < $redirectorFileAdr
