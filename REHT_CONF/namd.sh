#!/bin/bash

rep=$1
corenum=`expr $rep - 1`

#./charmrun +p1 ./namd2 pgl_mem${rep}.namd > pgl_mem${rep}.out
./namd2 +p1 pgl_mem${rep}.namd > myTmpFs/pgl_mem${rep}.out
retcode1=$?
#echo "RETURNED $retcode1"
size_coor=`wc -c pgl_mem${rep}.coor | awk '{print $1}'`
size_vel=`wc -c pgl_mem${rep}.vel | awk '{print $1}'`
nenergylines=`cat myTmpFs/pgl_mem${rep}.out | grep '^ENERGY' | wc -l`
nwallclock=`cat myTmpFs/pgl_mem${rep}.out | grep '^WallClock' | wc -l`
if [ $size_coor -ne 604420 ] || [ $size_vel -ne 604420 ] || [ $nenergylines -ne 2 ] || [ $nwallclock -ne 1 ]
then
 echo "ERROR in namd", $size_coor, $size_vel, $nenergylines, $nwallclock
 retcode1=1
fi


#./charmrun +p1 ./namd2 pgl_mem${rep}_swap.namd > myTmpFs/pgl_mem${rep}_swap.out
./namd2 +p1  pgl_mem${rep}_swap.namd > myTmpFs/pgl_mem${rep}_swap.out
retcode2=$?
size_coor=`wc -c pgl_mem${rep}_swap.coor | awk '{print $1}'`
size_vel=`wc -c pgl_mem${rep}_swap.vel | awk '{print $1}'`
nenergylines=`cat myTmpFs/pgl_mem${rep}_swap.out | grep '^ENERGY' | wc -l`
nwallclock=`cat myTmpFs/pgl_mem${rep}_swap.out | grep '^WallClock' | wc -l`
if [ $size_coor -ne 604420 ] || [ $size_vel -ne 604420 ] || [ $nenergylines -ne 1 ] || [ $nwallclock -ne 1 ]
then
 echo "ERROR in swap namd", $size_coor, $size_vel, $nenergylines, $nwallclock
 retcode2=1
fi

if [ $retcode1 -ne 0 ] || [ $retcode2 -ne 0 ]
then
 touch ERROR_namd.txt
fi
retcode=0
if [ $retcode1 -ne 0 ] || [ $retcode2 -ne 0 ]
then
 retcode=1
fi
#  echo $retcode > namd${rep}.out
echo $retcode


