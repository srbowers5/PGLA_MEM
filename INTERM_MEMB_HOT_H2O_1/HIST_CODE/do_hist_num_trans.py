#! /usr/bin/python3
#
#    This calculaes the energy from the dcd files created by namd.
#
import os,sys
MIN_ARGS=8
for i in range(0,MIN_ARGS):
    print("ARG ", i, sys.argv[i])

if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   do_hist_num_trans.py ")
    exit()

for i in range(MIN_ARGS,len(sys.argv)):
    print ("add path", sys.argv[i])
    sys.path.append(sys.argv[i])

import read_file

inDataPrefix1 = sys.argv[1]
inDataPrefix2 = sys.argv[2]
inDataSuffix1 = sys.argv[3]
inDataSuffix2 = sys.argv[4]
inScalePrefix = sys.argv[5]
outFile = sys.argv[6]
numReps = int(sys.argv[7])


totPf = 0.0
outVal1 = 0.0
outVal2 = 0.0
outVal = 0.0

for rep in range(1, numReps+1):
    inDataFile1 = inDataPrefix1 + str(rep) + inDataSuffix1
    inDataFile2 = inDataPrefix2 + str(rep) + inDataSuffix2
    inScaleFile = inScalePrefix + "_" + str(rep) + ".dat"

    inDataList1 = read_file.read_file_list(inDataFile1)
    print("READ ", inDataFile1, len(inDataList1))
    inDataList2 = read_file.read_file_list(inDataFile2)
    print("READ ", inDataFile2, len(inDataList2))
    inScaleList = read_file.read_file_list(inScaleFile)
    print("READ ", inScaleFile, len(inScaleList))

    for i in range(len(inDataList1)):
        line = read_file.comp_space(inDataList1[i])
        fields = line.split(' ')
        val1 = 0
        val = 0
        for j in range(5,8):
            if (int(fields[j]) != 0):
                val = 1
                val1 = 1

        line = read_file.comp_space(inDataList2[i])
        fields = line.split(' ')
        val2 = 0
        for j in range(5,8):
            if (int(fields[j]) != 0):
                val = 1
                val2 = 1


        scaleVal = float(inScaleList[i])
        totPf += scaleVal

        outVal1 += (val1*scaleVal)
        outVal2 += (val2*scaleVal)
        outVal += (val*scaleVal)


outVal1 = outVal1 / totPf
outVal2 = outVal2 / totPf
outVal = outVal / totPf

outFp = open(outFile,"w")
outStr = str(outVal1) + " " + str(outVal2) + " "  + str(outVal) + "\n"
outFp.write(outStr)
print ("VAL ", outVal1, outVal2, outVal)

outFp.close()
