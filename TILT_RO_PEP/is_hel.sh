#!/bin/bash
tr=1
rep=1
range=150001-200000
for rep in {1..20..1}; do
    ./is_hel.py DATA/helUp${tr}_${rep}.dat DATA/isHelUp${tr}_${range}_${rep}.dat /SCHOOL/pyth_utils
    ./is_hel.py DATA/helLo${tr}_${rep}.dat DATA/isHelLo${tr}_${range}_${rep}.dat /SCHOOL/pyth_utils
done
