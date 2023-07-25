#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = output_name

make get_lipid_den

tr=1
rep=1
dir=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT
for rep in {1..20..1}; do
    cmFile=${dir}/interm_out/cm${tr}_${rep}b_150001-200000.dat
    rm cmList.txt 2>/dev/null
    echo $cmFile > cmList.txt
    echo "./get_lipid_den ${dir} $tr $rep pgl22_mem 150001 200000 cmList.txt myOUTFILE"
    ./get_lipid_den ${dir} $tr $rep pgl22_mem 150001 200000 cmList.txt myOUTFILE > lip_den.out
done


