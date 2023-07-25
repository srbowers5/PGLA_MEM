#!/bin/bash

#inFile = sys.argv[1]
#funcFile = sys.argv[2]
#plotFile = sys.argv[3]
#plotSpace = sys.argv[4]
#plotMax = sys.argv[5]

plotSpace=1000
plotMax=200000
numExp=2

inFile=autoc_func_com_rest_equ.dat
funcFile=NEW_DATA/best_rest_autoc_func_COM_equ.dat
plotFile=NEW_DATA/best_rest_autoc_plot_COM_equ.dat
./get_best_curve.py $inFile $funcFile $plotFile $plotSpace $plotMax /SCHOOL/pyth_utils

inFile=autoc_func_com_hot_equ.dat
funcFile=NEW_DATA/best_hot_autoc_func_COM_equ.dat
plotFile=NEW_DATA/best_hot_autoc_plot_COM_equ.dat
./get_best_curve.py $inFile $funcFile $plotFile $plotSpace $plotMax /SCHOOL/pyth_utils




