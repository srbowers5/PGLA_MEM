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
make get_tilt_ulr
make get_tilt_ulr_neg
myDir=/SCHOOL/DO_REST_PCPG_22_INS3
tr=3
prefix=pgl22_mem
first_str=25001
last_str=75000
rep=1
for rep in {1..12..1}; do
  first_aa=15
  last_aa=20
  use_only_helix=0
  ./get_tilt_ulr_neg $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roUpNEG${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roLoNEG${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out


last_str=75000
#  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roUp${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roLo${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out





  first_aa=6
  last_aa=14
  use_only_helix=0
#  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roUp${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roLo${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out

done




for rep in {1..12..1}; do
  first_aa=15
  last_aa=20
  use_only_helix=1
#  ./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roUp_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat roLo_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out

done

