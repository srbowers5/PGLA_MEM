#!/bin/bash


#inDataPrefix1 = sys.argv[1]
#inDataPrefix2 = sys.argv[2]
#inDataSuffix1 = sys.argv[3]
#inDataSuffix2 = sys.argv[4]
#inScalePrefix = sys.argv[5]
#outFile = sys.argv[6]
#numReps = int(sys.argv[7])



inDataPrefix1=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT/interm_out/lip_coord_num_PC_group1_  
inDataSuffix1=_150001-200000.dat
inDataPrefix2=/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT/interm_out/lip_coord_num_PG_group1_  
inDataSuffix2=_150001-200000.dat
inScalePrefix=hist_scale_1
outFile=Prob_pep_trans_1_hot.dat
numReps=20

echo "./do_hist_num_trans.py $inDataPrefix1  $inDataPrefix2  $inDataSuffix1  $inDataSuffix2 $inScalePrefix $outFile $numReps /SCHOOL/pyth_utils"


./do_hist_num_trans.py $inDataPrefix1  $inDataPrefix2  $inDataSuffix1  $inDataSuffix2 $inScalePrefix $outFile $numReps /SCHOOL/pyth_utils


