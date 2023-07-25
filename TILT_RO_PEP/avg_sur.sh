#!/bin/bash


#inRUPFile = sys.argv[1]
#inRLOFile = sys.argv[2]
#inZFile = sys.argv[2]
#outFile = sys.argv[4]

./avg_sur.py roLo3_25001-75000_15-20_1.dat roUp3_25001-75000_15-20_1.dat /SCHOOL/DO_REST_PCPG_22_INS3/interm_out/term_com3_all_1.dat roINS_avg_1.dat

./avg_sur.py roLo3_25001-75000_15-20_2.dat roUp3_25001-75000_15-20_2.dat /SCHOOL/DO_REST_PCPG_22_INS3/interm_out/term_com3_all_2.dat roINS_avg_2.dat

./avg_sur.py roLo3_25001-75000_15-20_3.dat roUp3_25001-75000_15-20_3.dat /SCHOOL/DO_REST_PCPG_22_INS3/interm_out/term_com3_all_3.dat roINS_avg_3.dat

rep=4
./avg_sur.py roLo3_25001-75000_15-20_$rep.dat roUp3_25001-75000_15-20_$rep.dat /SCHOOL/DO_REST_PCPG_22_INS3/interm_out/term_com3_all_$rep.dat roINS_avg_$rep.dat

