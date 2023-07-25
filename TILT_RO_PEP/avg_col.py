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

inFile = sys.argv[1]
outFile = sys.argv[2]

colVals = []
inputList = read_file.read_file_list(inFile)
print "INPUT ", inFile, len(inputList)


line = read_file.comp_space(inputList[0])
fields = line.split(' ')
num_col = len(fields)
for i in range(0,num_col):
    colVals.append(0.0)

for i in range(len(inputList)):
    line = read_file.comp_space(inputList[i])
    fields = line.split(' ')
    for j in range(0,num_col):
        colVals[j] += float(fields[j])

outFd = open(outFile, "w")
for i in range(0, num_col):
     outStr = str(colVals[i]/ len(inputList)) + "\n"
     outFd.write(outStr)

outFd.close()

print "DONE"
