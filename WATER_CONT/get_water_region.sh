#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = pdb_file_name
#    arg[3] = tr
#    arg[4] = rep
#    arg[5] = prefix
#    arg[6] = first_str
#    arg[7] = last_str
#    arg[8] = COM output_name
#    arg[9] = PEP_CONT output_prefix
#    arg[10] = PEP_CONT_AA output_prefix

tr=1
rep=1
firstStr=150001
lastStr=200000
dir=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT
pdbname=spt.pdb
outFileUp=water_regionUp${tr}_${firstStr}-${lastStr}_${rep}.dat
outFileLo=water_regionLo${tr}_${firstStr}-${lastStr}_${rep}.dat

rm get_water_region
make get_water_region
rep=1
lastRep=20
while [  $rep -le $lastRep ]; do
    outFile=DATA/water_region_COM${tr}_${firstStr}-${lastStr}_${rep}.dat
    outFileUp=DATA/water_regionUp${tr}_${firstStr}-${lastStr}_${rep}.dat
    outFileLo=DATA/water_regionLo${tr}_${firstStr}-${lastStr}_${rep}.dat

    echo "./get_water_region $dir $tr $rep pgl22_mem $firstStr $lastStr $outFile $outFileUp $outFileLo > get_wire${rep}.out"
    ./get_water_region $dir $tr $rep pgl22_mem $firstStr $lastStr $outFile $outFileUp $outFileLo > get_wire${rep}.out
    rep=`expr $rep + 1`
done

