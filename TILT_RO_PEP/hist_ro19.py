#! /usr/bin/python

import os,sys
MIN_ARGS=3
if (len(sys.argv) < MIN_ARGS):
    print "NOT ENOUGH argv:\n   lin_reg.py <input> <output>  <lib paths>"
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print "ADD PATH", i, sys.argv[i]

import read_file
import math


#
xVals = []
yVals = []

inFileU = sys.argv[1]
inFileL = sys.argv[2]
outFileU = sys.argv[3]
outFileL = sys.argv[4]
outFileA = sys.argv[5]

roVals = []
roValsU = []
roValsL = []
for i in range(0,60):
    roVals.append(0.0)
    roValsU.append(0.0)
    roValsL.append(0.0)


totU = 0.0
totL = 0.0
inputList = read_file.read_file_list(inFileU)
print "INPUT ", inFileU, len(inputList)
for i in range(len(inputList)):
    line = read_file.comp_space(inputList[i])
    fields = line.split(' ')
    ro15Val = float(fields[0])
    ro19Val = float(fields[4])
    roBucket = int(ro19Val / 6.0)
    roVals[roBucket] += 1.0
    roValsU[roBucket] += 1.0
    totU += 1.0

inputList = read_file.read_file_list(inFileL)
print "INPUT ", inFileL, len(inputList)
for i in range(len(inputList)):
    line = read_file.comp_space(inputList[i])
    fields = line.split(' ')
    ro15Val = float(fields[0])
    ro19Val = float(fields[4])
    roBucket = int(ro19Val / 6.0)
    roVals[roBucket] += 1.0
    roValsL[roBucket] += 1.0
    totL += 1.0

outFd = open(outFileU, "w") 
for i in range(0,60):
    outStr = str((i*6) + 4) + " " + str(round(roValsU[i],5)) + "\n"
    outFd.write(outStr)
outFd.close()

outFd = open(outFileL, "w")
for i in range(0,60):
    outStr = str((i*6) + 4) + " " + str(round(roValsL[i],5)) + "\n"
    outFd.write(outStr)
outFd.close()

outFd = open(outFileA, "w")
for i in range(0,60):
    outStr = str((i*6) + 4) + " " + str(round(roVals[i],5)) + "\n"
    outFd.write(outStr)
outFd.close()


print "DONE"
