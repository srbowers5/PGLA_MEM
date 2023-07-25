#! /usr/bin/python3

import os,sys
import math
import time

num_files = int(sys.argv[1])
num_col = int(sys.argv[2])
outFile = sys.argv[3]
#  Infiles = 4....(num_file+3)
MIN_ARGS = num_files + 4


if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   ./merge_files.py ")
    exit
print ("Args", MIN_ARGS, len(sys.argv))
for i in range(MIN_ARGS,len(sys.argv)):
    print ("add path", sys.argv[i])
    sys.path.append(sys.argv[i])

import read_file

inList = []
outList = []
for i in range(0, num_files):
    inList.append(read_file.read_file_list(sys.argv[i+4]))
    print ("READ ", sys.argv[i+4], len(inList[i]))

for i in range(len(inList[0])):
    outList.append([])
    for j in range(0, num_col):
        outList[i].append(0.0)

for i in range(len(inList)):
    for j in range(len(inList[i])):
        line = read_file.comp_space(inList[i][j])
        fields = line.split(' ')
        if (len(fields) != num_col):
            print ("SHORT LINE ", len(fields), num_col)
        for k in range(0, num_col):
            outList[j][k] += float(fields[k])

outFd = open(outFile, "w")
for j in range(len(outList)):
    outStr = ""
    for k in range(0, num_col):
        outStr += str(round( (outList[j][k]/num_files), 10)) + " "
    outStr = outStr[:-1] + "\n"
    outFd.write(outStr)
outFd.close()
