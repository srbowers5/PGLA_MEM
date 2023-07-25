#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = output_name
#./get_scd ../SCHOOL/DO_REST_PCPG_22_INS2 2 1 pgl22_mem 31001 81000 > scd.out

tr=1
firstStr=150001
lastStr=200000
inputDir=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT
prefix=pgl22_mem

rep=1
for rep in {1..20..1}; do
  echo "START get_scd $inputDir $tr $rep $prefix $firstStr $lastStr"
  ./get_scd $inputDir $tr $rep $prefix $firstStr $lastStr > scd${rep}.out
done



