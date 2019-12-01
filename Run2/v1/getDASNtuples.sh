# For each dataset listed in data_fileRedirector, bkgd_datasets, and 
# sig_fileRedirector, creates the text file with the list of all ntuples for the 
# dataset. 
# For MC, if there is more than one dataset that corresponds to the same subprocess, 
# the ntuples from each of those datasets are all contained in 1 file.

IFS=' ' # built in bash variable, the line splitting delimiter

input="data_fileRedirector"
rm -rf dataNtupleLists/*
while read -r datasetName channel
do
    if [[ "$datasetName" =~ \#.* ]]; then
        continue
    fi

    if [ ! -d dataNtupleLists/$channel ]; then # check if dir already made
        mkdir dataNtupleLists/$channel
    fi

    dasgoclient --query="file dataset=/$channel/$datasetName/NANOAOD" | cat - >> \
        dataNtupleLists/$channel/$datasetName

    echo dataNtupleLists/$channel/$datasetName updated.
done < $input


echo
input="bkgd_datasets"
rm -rf bkgdNtupleLists/*
while read -r process subprocess datasetName 
do
    if [[ "$process" =~ \#.* ]]; then
        continue
    fi

    if [ ! -d bkgdNtupleLists/$process ]; then # check if dir already made
        mkdir bkgdNtupleLists/$process
    fi

    # Run dasgoclient and pipe the output to cat, which will take input from stdin and
    # append to the output file.
    dasgoclient --query="file dataset=$datasetName" | cat - >> \
        bkgdNtupleLists/$process/$subprocess

    echo bkgdNtupleLists/$process/$subprocess updated. 
done < $input

echo
input="sig_datasets"
rm -rf sigNtupleLists/*
while read -r process subprocess datasetName 
do
    if [[ "$process" =~ \#.* ]]; then
        continue
    fi

    if [ ! -d sigNtupleLists/$process ]; then # check if dir already made
        mkdir sigNtupleLists/$process
    fi

    # Run dasgoclient and pipe the output to cat, which will take input from stdin and
    # append to the output file.
    dasgoclient --query="file dataset=$datasetName" | cat - >> \
        sigNtupleLists/$process/$subprocess

    echo sigNtupleLists/$process/$subprocess updated. 
done < $input

echo
echo Done.
