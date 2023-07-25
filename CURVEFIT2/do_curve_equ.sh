#!/bin/bash

#IN_DATA_PREFIX=$1
#OUT_DATA_PREFIX=$2
#valcol=$3

inFile=/home/dad/MERGE_DIR/MERGE_AUTOC/rep_autoc_rest_com_all_equ.dat
outFile=autoc_func_com_rest_equ.dat
numExp=2
valCol=1
OUT_DATA_DIR=TMP_DATA/


rm $outFile 2>/dev/null
touch $outFile
for ratioa0_1 in {-20..20..1}; do
    for ratio in {2..40..2}; do
        ./do_curve.py $inFile ${OUT_DATA_DIR}tmpOut 2 $ratio $ratio $ratioa0_1 0 $valCol /SCHOOL/pyth_utils
        cat $outFile ${OUT_DATA_DIR}tmpOut > ${OUT_DATA_DIR}outFileTmp
        mv ${OUT_DATA_DIR}outFileTmp $outFile
        rm ${OUT_DATA_DIR}tmpOut
    done
done


inFile=/home/dad/MERGE_DIR/MERGE_AUTOC/rep_autoc_hot_com_all_equ.dat
outFile=autoc_func_com_hot_equ.dat
rm $outFile 2>/dev/null
touch $outFile
for ratioa0_1 in {-20..20..1}; do
    for ratio in {2..40..2}; do
        ./do_curve.py $inFile ${OUT_DATA_DIR}tmpOut 2 $ratio $ratio $ratioa0_1 0 $valCol /SCHOOL/pyth_utils
        cat $outFile ${OUT_DATA_DIR}tmpOut > ${OUT_DATA_DIR}outFileTmp
        mv ${OUT_DATA_DIR}outFileTmp $outFile
        rm ${OUT_DATA_DIR}tmpOut
    done
done

