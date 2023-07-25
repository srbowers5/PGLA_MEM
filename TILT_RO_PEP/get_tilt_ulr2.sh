#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = Use only helix
#    arg[8] = first_aa
#    arg[9] = last_aa
#    arg[10] = upper helix file
#    arg[11] = lower helix file

#    arg[7] = output_name
tr=1
make get_tilt_ulr
myDir=/SCHOOL//DO_REST_PCPG_HOT_H2O_1_CONT
prefix=pgl22_mem
first_str=1
last_str=150001


rep=1
for rep in {1..1..1}; do
  first_aa=15
  last_aa=20
  use_only_helix=0
  echo "./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} DATA/helUp${tr}_${rep}.dat DATA/helLo${tr}_${rep}.dat DATA/tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roUp${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roLo${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat"
  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} DATA/helUp${tr}_${rep}.dat DATA/helLo${tr}_${rep}.dat DATA/tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roUp${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roLo${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out



  first_aa=6
  last_aa=14
  use_only_helix=0
#  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} DATA/helUp${tr}_${rep}.dat DATA/helLo${tr}_${rep}.dat DATA/tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roUp${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roLo${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out

done




for rep in {1..1..1}; do
  first_aa=15
  last_aa=20
  use_only_helix=1
#  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} DATA/helUp${tr}_${rep}.dat DATA/helLo${tr}_${rep}.dat DATA/tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roUp_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roLo_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out

  first_aa=6
  last_aa=14
#  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} DATA/helUp${tr}_${rep}.dat DATA/helLo${tr}_${rep}.dat DATA/tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roUp_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat DATA/roLo_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt4.out

done

