#! /usr/bin/python3

import os,sys
MIN_ARGS=3
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   get_tot_cont.py  ")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file
import math

inPrefix = sys.argv[1]
outFile = sys.argv[2]

inUpFile=inPrefix + "up.dat"
inLoFile=inPrefix + "lo.dat"

upList = read_file.read_file_list(inUpFile)
loList = read_file.read_file_list(inLoFile)
print ("READ ", inUpFile, len(upList))
print ("READ ", inLoFile, len(loList))

outFd = open(outFile, "w")

for i in range(len(upList)):
    lineUp = read_file.comp_space(upList[i])
    fields = lineUp.split(' ')
    upVal = 0
    for j in range(len(fields)):
        upVal += int(fields[j])


    lineLo = read_file.comp_space(loList[i])
    fields = lineLo.split(' ')
    loVal = 0
    for j in range(len(fields)):
        loVal += int(fields[j])
    outStr = str(upVal) + " " + str(loVal) + "\n"
    outFd.write(outStr)

outFd.close()
   



print ("DONE get_tot_cont.py ", outFile)
