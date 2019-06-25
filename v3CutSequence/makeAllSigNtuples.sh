#!/bin/bash

# For each channel, creates the sig ntuples.

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
