#!/bin/bash

# For each channel, creates the sig ntuples.
# Args to makeNtupleSigs: test mode, same flav, mu pref

echo "python makeNtupleSigs.py 0 0 0"
python makeNtupleSigs.py 0 0 0
echo
echo
echo "python makeNtupleSigs.py 0 1 0"
python makeNtupleSigs.py 0 1 0
echo
echo
echo "python makeNtupleSigs.py 0 1 1"
python makeNtupleSigs.py 0 1 1
