#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = outFile
tr=1
rep=1
rm get_charge_den_z
make get_charge_den_z

start=150001
#end=150003
end=200000
dir=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT
outFile=${dir}/CHARGE_DEN/charge_Z${tr}_${rep}_${start}-${end}.dat
echo "./get_charge_den_z $dir $tr $rep pgl22_mem $start $end $outFile > charge_den${rep}.out"

for rep in {1..1..1}; do
    ./get_charge_den_z $dir $tr $rep pgl22_mem $start $end $outFile > charge_den${rep}.out
done


