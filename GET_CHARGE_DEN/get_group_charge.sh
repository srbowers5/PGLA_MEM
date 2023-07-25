#!/bin/bash

#
#    arg[1] = Input_dir
#    arg[2] = tr
#    arg[3] = rep
#    arg[4] = prefix
#    arg[5] = first_str
#    arg[6] = last_str
#    arg[7] = output_name
rm get_group_charge
make get_group_charge


inPsf=/SCHOOL/DO_REST_MEMB1_CONT/mem22_pep.psf
outFile=Group_charge.dat

./get_group_charge $inPsf $outFile
