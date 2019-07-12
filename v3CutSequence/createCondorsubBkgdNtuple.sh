# Creates a condorsub file for makeNtupleBkgd with $1-4 as the arguments to 
# makeNtupleBkgd.py

# cat > testcondorsub <<EOF
cat > condorsub_makeNtupleBkgd <<EOF
executable  = makeNtupleBkgd.py
arguments   = $1 $2 $3 $4

universe    = vanilla
output      = ../../condorLogs/makeNtupleBkgd.\$(ClusterId).out
error       = ../../condorLogs/makeNtupleBkgd.\$(ClusterId).err
log         = ../../condorLogs/makeNtupleBkgd.\$(ClusterId).log

+MaxRuntime = 100000

transfer_input_files = stopSelection.py, bkgd_files, bkgdProcesses/$4

queue
EOF

echo "makeNtupleBkgd.py $1 $2 $3 $4"
