#!/bin/bash


#clustStateFile = sys.argv[1]
#contFileUp = sys.argv[2]
#contFileLo = sys.argv[3]
#outFile = sys.argv[4]

tr=1
range=150001-200000

outFile=pep_cont_depth_state.dat

rep=1
lastRep=1
while [  $rep -le $lastRep ]; do
    clustStateFile=/home/dad/GET_CLUST_STATE/cont_state${tr}_${range}_${rep}.dat
    contFileUp=/home/dad/WATER_CONT/DATA/pep_cont${tr}_${range}_${rep}_up.dat
    contFileLo=/home/dad/WATER_CONT/DATA/pep_cont${tr}_${range}_${rep}_lo.dat

    ./pep_cont_state.py $clustStateFile $contFileUp $contFileLo $outFile /SCHOOL/pyth_utils3
    rep=`expr $rep + 1`
done

