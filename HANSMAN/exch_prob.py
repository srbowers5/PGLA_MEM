#! /usr/bin/python3
#   Sum the structures for each AA.
#   The input files should be the output from get_ss.sh and 
#   should have 1 letter for each AA for each structure

import os,sys
MIN_ARGS=4
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   exch_prob.py <inp_file> <output_file> <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file



inFile=sys.argv[1]
numRep=int(sys.argv[2])
outFile=sys.argv[3]

numExch = []
for i in range(0, numRep):
    numExch.append(0.0)

currVals = []
currVals.append('0')
for i in range(0, numRep):
    currVals.append(i+1)

inList=read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))
cnt = 0
for i in range (len(inList)):
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    for j in range(1, len(fields)):
        val = int(fields[j])
        if (val != currVals[j]):
            numExch[j-1] += 1.0
            currVals[j] = val
    cnt += 1

outFd = open(outFile,"w")
for i in range(0, numRep):
    if ( (i == 0) | (i == (numRep-1)) ):
        val = numExch[i] / (cnt/2.0)
    else:
        val = numExch[i] / cnt
    outStr = str(i+1) + " " + str(round(val, 4)) + "\n"
    outFd.write(outStr)
outFd.close()

print ("DONE exch_prob.py", outFile)

