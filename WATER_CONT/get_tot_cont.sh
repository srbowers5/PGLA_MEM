#!/bin/bash

#inPrefix = sys.argv[1]
#outFile = sys.argv[2]

tr=1
rep=1
range="150001-200000"

util_path=/SCHOOL/pyth_utils

for rep in {1..20..1}; do
    inPrefix=DATA/pep_cont${tr}_${range}_${rep}_
    outFile=DATA/pep_cont_sum${tr}_${range}_${rep}.dat
    ./get_tot_cont.py $inPrefix $outFile $util_path
done
