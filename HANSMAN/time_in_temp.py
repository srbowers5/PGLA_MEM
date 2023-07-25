#! /usr/bin/python3
#   Sum the structures for each AA.
#   The input files should be the output from get_ss.sh and 
#   should have 1 letter for each AA for each structure
#   For example one line could be CCCCCHHHHHBBBBBTTTTTT

import os,sys
MIN_ARGS=4
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
numCol=int(sys.argv[2])
outFile=sys.argv[3]

time_in_temp_array = []
for i in range(0, numCol):
    time_in_temp_array.append([])
    for j in range(0, numCol):
        time_in_temp_array[i].append(0.0)

inList=read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))
cnt = 0
for i in range (len(inList)):
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    for j in range(1, len(fields)):
        off = int(fields[j]) - 1
        time_in_temp_array[j-1][off] += 1
    cnt += 1

outFd = open(outFile,"w")
for i in range(0, numCol):
    outStr = str(i+1) + " "
    for j in range(0, numCol):
        outStr += str(time_in_temp_array[i][j]) + " "  
    outStr = outStr[:-1] + "\n"
    outFd.write(outStr)
outFd.close()

print ("DONE time_in_temp.py", outFile)

