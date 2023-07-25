#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = COM output_name
#    arg[8] = PEP_CONT output_prefix
#    arg[9] = PEP_CONT_AA output_prefix
#    argv[10] = USE AA or WAT depth. 0=Water, 1=AA

tr=1
rep=1
firstStr=150001
lastStr=200000
dir=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT
use_aa=0

rep=1
lastRep=20
while [  $rep -le $lastRep ]; do
    xscOutFile=DATA/com_xsc${tr}_${firstStr}-${lastStr}_${rep}.dat
    pepContPrefix=DATA/pep_cont${tr}_${firstStr}-${lastStr}_${rep}
    pepAAContPrefix=DATA/pep_contAA${tr}_${firstStr}-${lastStr}_${rep}_
    echo "./get_water_cont $dir $tr $rep pgl22_mem $firstStr $lastStr  $xscOutFile $pepContPrefix \
      $pepAAContPrefix $use_aa > get_com${rep}.out"
    ./get_water_cont $dir $tr $rep pgl22_mem $firstStr $lastStr  $xscOutFile $pepContPrefix \
       $pepAAContPrefix $use_aa > get_com${rep}.out
    rep=`expr $rep + 1`
done

