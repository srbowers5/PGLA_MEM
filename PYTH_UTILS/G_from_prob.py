#! /usr/bin/python3

import os,sys
MIN_ARGS=5
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   G_from_prob.py <input> <output>  <col> <temp> <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file
import math

inFile = sys.argv[1]
outFile = sys.argv[2]
col = int(sys.argv[3])
tempVal = float(sys.argv[4])

outFd = open(outFile, "w")

inputList = read_file.read_file_list(inFile)

RVal = 1.987    # cal/mol/K
# Find Max
maxVal1 = 0.0
maxVal2 = 0.0
sumVal = 0.0


for i in range(len(inputList)):
    line = read_file.comp_space(inputList[i])
    fields = line.split()

    val = float(fields[col])
    if (val > maxVal1):
        maxVal1 = val

outFd = open(outFile,"w")
for i in range(len(inputList)):
    line = read_file.comp_space(inputList[i])
    fields = line.split()

    val1 = float(fields[col])
    if (val1 <= 0.0):
        feVal1 = 99999.0
    elif (val1 == maxVal1):
        feVal1 = 0.0
    else:
        feVal1 = math.log(float(fields[col]) / maxVal1)
        feVal1 = -(feVal1 * RVal * tempVal) / 1000



    outStr = fields[0] + " " + str(round(feVal1,4)) + "\n"
    outFd.write(outStr)


outFd.close()


print ("DONE G_from_prob.py ", outFile)

