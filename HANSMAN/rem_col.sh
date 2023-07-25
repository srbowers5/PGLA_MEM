#!/bin/bash 

#inFile = sys.argv[1]
#outFile = sys.argv[2]

tr=2
add_path=/SCHOOL/pyth_utils

echo "./rem_col.py Reps_all_${tr}.txt ran_walk${tr}.dat $add_path"

./rem_col.py Reps_all_${tr}.txt ran_walk${tr}.dat $add_path
