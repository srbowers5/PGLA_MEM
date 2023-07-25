#! /usr/bin/python3

import os,sys
MIN_ARGS=6
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   get_best_curve.py <input> <output>  <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file
import math
import numpy
from scipy.optimize import curve_fit


inFile =    sys.argv[1]
funcFile =  sys.argv[2]
plotFile =  sys.argv[3]
plotSpace = float(sys.argv[4])
plotMax =   float(sys.argv[5])


inList = read_file.read_file_list(inFile)
print ("READ ",inFile, len(inList))

funcList = []
for i in range(len(inList)):
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    print ("FIELDS ", len(fields))
    funcList.append( (float(fields[13]), line) )

funcList.sort()

goodList = []

for i in range(len(funcList)):
    line = read_file.comp_space(funcList[i][1])
    fields = line.split(' ')
    a0 = abs(float(fields[4]))
    a1 = abs(float(fields[5]))
    a2 = abs(float(fields[6]))
    a3 = abs(float(fields[7]))
    a4 = abs(float(fields[8]))
    totVal = a0 + a1 + a2 + a3 + a4
    diffVal = abs(1 - totVal)
    if (diffVal < 0.01):
        goodList.append(funcList[i])


for i in range(len(goodList)):
    print (" GOOD LIST ", goodList[i])

outFd = open(funcFile,"w")
for i in range(len(goodList)):
    outStr = goodList[i][1] + "\n"
    outFd.write(outStr)
outFd.close()
#
#  DO plot of best
line = read_file.comp_space(goodList[0][1])
print ("LINE ", line)
fields = line.split(' ')
a0 = abs(float(fields[4]))
a1 = abs(float(fields[5]))
a2 = abs(float(fields[6]))
a3 = abs(float(fields[7]))
a4 = abs(float(fields[8]))
t1 = abs(float(fields[9]))
t2 = abs(float(fields[10]))
t3 = abs(float(fields[11]))
t4 = abs(float(fields[12]))


outFd = open(plotFile,"w")

tVal = 0
print ("A0 = ", a0, "A1 = ", a1, "A2 = ", a2)

while (tVal < plotMax):
    yVal = a0 + (a1 * numpy.exp(-tVal/t1)) + (a2 * numpy.exp(-tVal/t2)) + (a3 * numpy.exp(-tVal/t3)) + (a4 * numpy.exp(-tVal/t4))
    outStr = str(tVal) + " " + str(round(yVal,6)) + "\n"
    outFd.write(outStr)
    tVal += plotSpace

outFd.close()



print ("DONE ")
