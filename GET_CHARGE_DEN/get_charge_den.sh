#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = output_name
tr=1
rep=1
rm get_charge_den
make get_charge_den

start=150001
end=200000
dir=/External_5/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT
outDir=CHARGE_DEN/
echo "./get_charge_den $outDir $tr $rep pgl22_mem $start $end $dir > charge_den_hot${rep}.out"

for rep in {3..20..1}; do
    ./get_charge_den $outDir $tr $rep pgl22_mem $start $end $dir > charge_den_hot${rep}.out
done


#start=195001
#end=245000
#dir=/SCHOOL/DO_REST_MEMB1_CONT
#outDir=${dir}/CHARGE_DEN/
#echo "./get_charge_den $outDir $tr $rep pgl22_mem $start $end $dir > charge_den${rep}.out"
#
#for rep in {1..12..1}; do
#    ./get_charge_den $outDir $tr $rep pgl22_mem $start $end $dir > charge_den${rep}.out
##done
