#! /usr/bin/python

import os,sys
import math
import time
MIN_ARGS=3
if (len(sys.argv) < MIN_ARGS):
    print "NOT ENOUGH argv:\n   ./get_exch.py <traj> <lib paths>"
    exit
print "Args", MIN_ARGS, len(sys.argv)
for i in range(MIN_ARGS,len(sys.argv)):
    print "add path", sys.argv[i]
    sys.path.append(sys.argv[i])

import read_file




inFile = sys.argv[1]
outFile = sys.argv[2]

inList = read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))

fdOut = open(outFile,"w")
for i in range(len(inList)):
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    outStr = ""
    for j in range(1, len(fields)):
        outStr = outStr + fields[j] + " "
    outStr = outStr[:-1] + "\n"
    fdOut.write(outStr)

fdOut.close()


