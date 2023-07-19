
#! /bin/sh
#
#   $1 is batch
#   $2 maxRep
#....attempting exchange.....

 rm RepEnergy.dat 2>/dev/null
 rm RepEnergySwap.dat 2>/dev/null
 rm RepBox.dat 2>/dev/null
 rm RepBoxSwap.dat 2>/dev/null
 rep=1
 echo 'Batch is ' ${1} 'Max rep is ' ${2}

 while [ $rep -le ${2} ]
 do
#    rm pgl_mem${rep}.namd
#    cat pgl_mem${rep}.out | grep '^TCL' |  cat -n >> forces${rep}_$2.dat
#    echo 'Get Energies' $rep
    cat myTmpFs/pgl_mem${rep}.out | grep '^ENERGY' | cut -c144-179 | cat -n >> energy${rep}_$1.dat
#    cp pgl_mem${rep}.out pgl_mem${rep}.out.SAV

    cat energy${rep}_$1.dat | tail -1 | cut -c8- >> RepEnergy.dat
    cat myTmpFs/pgl_mem${rep}_swap.out | grep '^ENERGY' | cut -c144-179 | cat -n >> energy${rep}_$1_swap.dat
    cat energy${rep}_$1_swap.dat | tail -1 | cut -c8- >> RepEnergySwap.dat

    tail -1 pgl_mem${rep}.xsc > temp1
    read dum x dum dum dum y dum dum dum z dum dum dum < temp1
    echo $x $y $z >> RepBox.dat
    cat temp1 >> xsc${rep}_$1.dat
    rm temp1

    tail -1 pgl_mem${rep}_swap.xsc > temp1
    read dum x dum dum dum y dum dum dum z dum dum dum < temp1
    echo $x $y $z >> RepBoxSwap.dat
    cat temp1 >> xsc${rep}_$1_swap.dat
    rm temp1

    rep=`expr $rep + 1`
 done
 ./exchange_npt.exe < status.dat > exchange_out
 retcode=$?
if [ $retcode -ne 0 ]
 then
  echo 'Exchange error.'
  touch 'ERROR_exchange'
  exit
 fi

 npairs=`cat exchange_out | wc -l`
 np=1
 while [ $np -le $npairs ]
 do
  head -$np exchange_out|tail -1 > exchange_out_temp
  read rep1 rep2 iexchange < exchange_out_temp
  if [ $iexchange -eq 1 ]
  then
    mv pgl_mem${rep1}.coor restart_coor1
    mv pgl_mem${rep2}.coor restart_coor2
    mv restart_coor1 pgl_mem${rep2}.coor
    mv restart_coor2 pgl_mem${rep1}.coor

    mv pgl_mem${rep1}s.vel pgl_mem${rep1}.vel
    mv pgl_mem${rep2}s.vel pgl_mem${rep2}.vel

    mv pgl_mem${rep1}.xsc restart_xsc1
    mv pgl_mem${rep2}.xsc restart_xsc2
    mv restart_xsc1 pgl_mem${rep2}.xsc
    mv restart_xsc2 pgl_mem${rep1}.xsc
  fi
  echo $np $rep1 $rep2 'exchange:' $iexchange >>  RunStatus$1
  np=`expr $np + 1`
 done

