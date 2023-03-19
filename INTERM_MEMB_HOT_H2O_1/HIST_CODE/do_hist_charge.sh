#!/bin/bash


#inDataPrefix = sys.argv[1]
#inDataSuffix = sys.argv[2]
#inScalePrefix = sys.argv[3]
#outFile = sys.argv[4]
#numReps = int(sys.argv[5])
#numCol = int(sys.argv[6])
#skipFirstCol = sys.argv[7]

numReps=20

./do_hist_data.py /SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT/CHARGE_DEN/charge_reg_wat_1_ _150001-200000.dat  hist_scale_1 charge_reg_wat_1_hot.dat $numReps 9 0 /SCHOOL/pyth_utils
