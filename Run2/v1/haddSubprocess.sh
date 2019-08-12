#!/bin/bash

# NOTE: NEEDS 5 CMD LINE ARGS with values: 
# testMode {test, all}, input type {data, bkgd, sig}, channel {mumu, elel, muel}, 
# subprocess, process
# Process options: 
#   data: DoubleMuon, DoubleEG, MuonEG
#   bkgd: TTBar TT+X Diboson W-Jets Drell-Yan Single-Top
#   sig: Stop-Pair
# Subprocess: name of the dataset for data and sig, or name of the subprocess for bkgd
# 
# For a particular subprocess, hadds together all the output root files from
# makeNtuple.py.

testMode=$1
inputType=$2
channel=$3
if [[ "$channel" == "mumu" ]]; then
    channel="MuMu"
elif [[ "$channel" == "elel" ]]; then
    channel="ElEl"
elif [[ "$channel" == "muel" ]]; then
    channel="MuEl"
fi

subprocess=$4
process=$5

myDataDir=/afs/cern.ch/work/c/cmiao/private/myDataSusy
inDir=$myDataDir/Run2/$inputType/$process/$subprocess
outDir=$myDataDir/Run2/$inputType/$process/$subprocess
if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

numInFiles="$(ls -1 $inDir/stopCut_$testMode*$channel".root" | wc -l)"

if [ $numInFiles  == 0 ]; then
    echo No files to hadd, skipping
    exit 1
fi

echo Num in files: $numInFiles

hadd -f -j -k $outDir/$subprocess"_"$testMode"_"$channel".root" \
    $inDir/stopCut_$testMode*$channel".root"