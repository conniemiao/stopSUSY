#!/bin/bash

# NOTE: THIS DOESN'T WORK. LOL. Right now just makes the sigs ntuples.
# for each setting combo for creating bkgd and sig ntuples, creates a condorsub
# file for it and then submits it.

echo "python makeNtupleSigs.py 0 0 0 0"
python makeNtupleSigs.py 0 0 0 0
echo
echo "python makeNtupleSigs.py 0 1 0 0"
python makeNtupleSigs.py 0 1 0 0
echo
echo "python makeNtupleSigs.py 0 0 1 0"
python makeNtupleSigs.py 0 0 1 0
echo
echo "python makeNtupleSigs.py 0 1 1 0"
python makeNtupleSigs.py 0 1 1 0
echo
echo "python makeNtupleSigs.py 0 0 1 1"
python makeNtupleSigs.py 0 0 1 1
echo
echo "python makeNtupleSigs.py 0 1 1 1"
python makeNtupleSigs.py 0 1 1 1

# echo
# echo "python makeNtupleBkgd.py 1 0 0 0"
# python makeNtupleBkgd.py 1 0 0 0
# echo
# echo "python makeNtupleBkgd.py 1 1 0 0"
# python makeNtupleBkgd.py 1 1 0 0
# echo
# echo "python makeNtupleBkgd.py 1 0 1 0"
# python makeNtupleBkgd.py 1 0 1 0
# echo
# echo "python makeNtupleBkgd.py 1 1 1 0"
# python makeNtupleBkgd.py 1 1 1 0
# echo
# echo "python makeNtupleBkgd.py 1 0 1 1"
# python makeNtupleBkgd.py 1 0 1 1
# echo
# echo "python makeNtupleBkgd.py 1 1 1 1"
# python makeNtupleBkgd.py 1 1 1 1
