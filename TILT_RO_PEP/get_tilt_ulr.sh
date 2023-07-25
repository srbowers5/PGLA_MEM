
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

myDir=/SCHOOL/DO_REST_PCPG_22_INS3/
tr=3
prefix=pgl22_mem
first_str=25001
#last_str=75000
last_str=25020


rep=1
first_aa=15
last_aa=20
use_only_helix=0
./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#

#rep=2
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
 

#rep=3
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#

#rep=4
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
###./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#


#rep=5
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#


#rep=6
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#


#rep=7
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#
#
#rep=8
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#
#
#rep=9
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#
#
#rep=10
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#
#
#rep=11
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#
#
#rep=12
#first_aa=15
#last_aa=20
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#use_only_helix=1
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt_hel${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#first_aa=6
#last_aa=14
#use_only_helix=0
#./get_tilt_ulr $myDir $tr ${rep} ${prefix} ${first_str} ${last_str} ${use_only_helix} ${first_aa} ${last_aa} helUp${tr}_${rep}.dat helLo${tr}_${rep}.dat tilt${tr}_${first_str}-${last_str}_${first_aa}-${last_aa}_${rep}.dat > tilt.out
#
#
