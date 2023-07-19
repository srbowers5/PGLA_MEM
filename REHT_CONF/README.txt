This folder contains the files needed to run the REHT (Replica Exchange with Hybrid Tempering).
The initial files are in initSTR2, initSTR3, or REHT_initSTR1.
initSTR2 and initSTR3 are identical to the files used for the pervious REST simulations. Replicas 13-20 repeat the use of replicas 1-8.
initSTR1 has 20 replicas taken from repica 1 of the REST simulation. The program reptamp.sh create the files ScaleFactor.dat and RepTemp.dat.
The program schedule.sh create the file schedule.dat.
THe file spt.pdb sets which atoms are tempered. 
To run the program run ./rem_mono.sh.


