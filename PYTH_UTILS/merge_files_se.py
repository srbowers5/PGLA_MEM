#! /usr/bin/python3

import os,sys
import math
import time

num_files = int(sys.argv[1])
num_col = int(sys.argv[2])
skip_first_col = sys.argv[3]
outFile = sys.argv[4]
#  Infiles = 5....(num_file+4)
MIN_ARGS = num_files + 5


if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   ./get_depth_aa.py ")
    exit
print ("Args", MIN_ARGS, len(sys.argv))
for i in range(MIN_ARGS,len(sys.argv)):
    print ("add path", sys.argv[i])
    sys.path.append(sys.argv[i])

import read_file
import get_sd

inList = []
outList = []
for i in range(0, num_files):
    inList.append(read_file.read_file_list_rem_blank(sys.argv[i+5]))

for i in range(len(inList[0])):
    outList.append([])
    for j in range(0, num_col):
        outList[i].append([])

for i in range(len(inList)):
    for j in range(len(inList[i])):
        line = read_file.comp_space(inList[i][j])
        fields = line.split(' ')
#        print ("LINE (", line, ")")
        if (len(fields) != num_col):
            print ("SHORT LINE ", len(fields), num_col)
            continue
        for k in range(0, num_col):
            if ((k==0) & (skip_first_col == '1')):
                outList[j][k].append(fields[k])
            else:
#                print("ADD ", j, k, fields[k]))
                outList[j][k].append(float(fields[k]))

outFd = open(outFile, "w")
for j in range(len(outList)):
    outStr = ""
    for k in range(0, num_col):
        if ((k==0) & (skip_first_col == '1')):
            outStr += str(outList[j][k][0]) + " "
        else:
            sdVal, seVal, meanVal, numSamp = get_sd.get_sd(outList[j][k])
            outStr += str(round( meanVal, 8)) + " " + str(round( seVal, 8)) + " "
    outStr = outStr[:-1] + "\n"
    outFd.write(outStr)
#    print ("WRITE ", outStr)
outFd.close()
