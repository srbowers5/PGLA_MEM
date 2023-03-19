#! /usr/bin/python3
#
#    This calculaes the energy from the dcd files created by namd.
#
import os,sys
MIN_ARGS=8
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   fmt_cm.py <tr> <rep> <dir> <template> <prefix>")
    exit()

for i in range(0,MIN_ARGS):
    print("ARG ", i, sys.argv[i])

for i in range(MIN_ARGS,len(sys.argv)):
    print ("add path", sys.argv[i])
    sys.path.append(sys.argv[i])

import read_file

inDataPrefix = sys.argv[1]
inDataSuffix = sys.argv[2]
inScalePrefix = sys.argv[3]
outFile = sys.argv[4]
numReps = int(sys.argv[5])
numCol = int(sys.argv[6])
skipFirstCol = sys.argv[7]


outVals = []
for i in range(0, numCol):
    outVals.append(0.0)

totPf = 0.0




for rep in range(1, numReps+1):
    inDataFile = inDataPrefix + str(rep) + inDataSuffix
    inScaleFile = inScalePrefix + "_" + str(rep) + ".dat"

    inDataList = read_file.read_file_list(inDataFile)
    print("READ ", inDataFile, len(inDataList))
    inScaleList = read_file.read_file_list(inScaleFile)
    print("READ ", inScaleFile, len(inScaleList))

    for i in range(len(inDataList)):
        line = read_file.comp_space(inDataList[i])
        fields = line.split(' ')
        scaleVal = float(inScaleList[i])
        totPf += scaleVal
        if (skipFirstCol == 1):
            firstCol = 1
        else:
            firstCol = 0
        for j in range(numCol):
            outVals[j] += (float(fields[j+firstCol])*scaleVal)


for i in range(0, numCol):
    outVals[i] = outVals[i] / totPf

outFp = open(outFile,"w")
for i in range(0, numCol):
    outStr = str(outVals[i]) + "\n"
    print ("VALS ", outVals[i])
    outFp.write(outStr)

outFp.close()
