#! /usr/bin/python3

import os,sys
MIN_ARGS=3
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   is_hel.py <input> <output>  <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file
import math


#
#
inFile = sys.argv[1]
outFile = sys.argv[2]



#
# main
outFd = open(outFile, "w")
inputList = read_file.read_file_list(inFile)
print ("INPUT ", inFile, len(inputList))
for i in range(len(inputList)):
    line = read_file.comp_space(inputList[i])
    if (line[13:18] == "HHHHH"):
        outStr = "1\n"
    else:
        outStr = "0\n"
    outFd.write(outStr)
outFd.close()

print ("DONE ishel.py", outFile)
