#! /usr/bin/python

import os,sys
MIN_ARGS=5
if (len(sys.argv) < MIN_ARGS):
    print "NOT ENOUGH argv:\n   avg_sur.py <input> <output>  <lib paths>"
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

inRUPFile = sys.argv[1]
inRLOFile = sys.argv[2]
inZFile = sys.argv[2]
outFile = sys.argv[4]

surcolVals = []
inscolVals = []
inputUPList = read_file.read_file_list(inRUPFile)
print "INPUT ", inRUPFile, len(inputUPList)
inputLOList = read_file.read_file_list(inRLOFile)
print "INPUT ", inRLOFile, len(inputLOList)
inputZList = read_file.read_file_list(inZFile)
print "INPUT ", inZFile, len(inputZList)


line = read_file.comp_space(inputUPList[0])
fields = line.split(' ')
num_col = len(fields)
for i in range(0,num_col):
    surcolVals.append(0.0)
    inscolVals.append(0.0)

insCnt = 0
surCnt = 0

for i in range(len(inputUPList)):
    line = read_file.comp_space(inputUPList[i])
    RUPfields = line.split(' ')
    line = read_file.comp_space(inputLOList[i])
    RLOfields = line.split(' ')
    line = read_file.comp_space(inputZList[i])
    Zfields = line.split(' ')

    if (float(Zfields[2]) < 17.553):
        for j in range(0,num_col):
            inscolVals[j] += float(RUPfields[j])
        insCnt += 1
    else:
        for j in range(0,num_col):
            surcolVals[j] += float(RUPfields[j])
        surCnt += 1

    if (-float(Zfields[4]) < 17.553):
        for j in range(0,num_col):
            inscolVals[j] += float(RLOfields[j])
        insCnt += 1
    else:
        for j in range(0,num_col):
            surcolVals[j] += float(RLOfields[j])
        surCnt += 1



outFd = open(outFile, "w")

for i in range(0, num_col):
     outStr = "INSERT  " + str(inscolVals[i]/ insCnt) + "\n"
     outFd.write(outStr)
     outStr = "SURFACE " + str(surcolVals[i]/ surCnt) + "\n"
     outFd.write(outStr)

outFd.close()

print "DONE"
