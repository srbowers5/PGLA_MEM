tr=1
nrep=20
strmin=1
strmax=200000

str=$strmin
rm schedule.dat
while [ $str -le $strmax ]
do
 echo $str $tr $nrep > status.dat
 ./schedule.exe < status.dat >> schedule.dat
 str=`expr $str + 1`
done

