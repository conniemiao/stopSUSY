# $1-4: arguments to makeNtupleBkgd.py

# cat > condorsub_makeNtupleBkgd <<EOF
cat > testcondorsub <<EOF
executable  = makeNtupleBkgd.py
arguments   = $1 $2 $3 $4

universe    = vanilla
output      = ../../condorLogs/makeNtupleBkgd.\$(ClusterId).\$(ProcId).out
error       = ../../condorLogs/makeNtupleBkgd.\$(ClusterId).\$(ProcId).err
log         = ../../condorLogs/makeNtupleBkgd.\$(ClusterId).log
+MaxRuntime = 10000

transfer_input_files = stopSelection.py, bkgd_files, bkgdProcesses/$4/*

queue
EOF

