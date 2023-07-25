#!/bin/bash
tr=1
range=150001-200000
dir=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT/interm_out

for rep in {2..20..1}; do
    ./get_hel_ulr ${dir}/dhWT${tr}_${rep}b_${range}.dat DATA/helUp${tr}_${rep}.dat DATA/helLo${tr}_${rep}.dat 
done

