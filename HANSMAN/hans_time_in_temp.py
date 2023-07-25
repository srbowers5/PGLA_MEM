#! /usr/bin/python3
#   Sum the structures for each AA.
#   The input files should be the output from get_ss.sh and 
#   should have 1 letter for each AA for each structure
#   For example one line could be CCCCCHHHHHBBBBBTTTTTT

import os,sys
MIN_ARGS=3
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   get_ss_aa.py <inp_file> <output_file> <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file


# ./mk_ran_walk.py $inFile $initFile $outFile /SCHOOL/pyth_utils

inFile=sys.argv[1]
outFile=sys.argv[2]

inList=read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))
outFd = open(outFile,"w")
for i in range (len(inList)):
    val1 = 0.0
    val2 = 0.0
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    for j in range(1, len(fields)):
        val2 += float(fields[j])
        val1 += (float(fields[j]) ** 2)
    val1 = val1 ** 0.5
    hansVal = 1 - (val1/val2)
    outStr = fields[0] + " " + str(round(hansVal,4)) + "\n"
    outFd.write(outStr)

outFd.close()

print ("DONE time_in_temp.py", outFile)

