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

inFile = sys.argv[1]
initFile = sys.argv[2]
outFile = sys.argv[3]

currRep = []
outFd = open(outFile, "w")

inList=read_file.read_file_list(initFile)
print ("READ ", initFile, len(inList))
line = ""
for i in range (len(inList)):
    line = read_file.comp_space(inList[i])
    outStr = line + "\n"
    outFd.write(outStr)
    print ("LINE ", line)
print ("Last LINE ", line)
fields = line.split(' ')
step = int(fields[0])
for i in range(0, len(fields)):
    currRep.append(fields[i]) 
    
lookCompl = True
inList=read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))
for i in range (len(inList)):
    print("I ", i)
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    if (lookCompl == True):
        if (len(fields) >= 4):
            print ("Len compl ", len(fields))
            if (fields[3] == "completed"):
                lookCompl = False
                print ("COMPLETED ")
        continue
    else:
        if (len(fields) >= 3):
            if (fields[0] == "Replica"):
                step += 1
                outStr = str(step) + " " 
                for j in range(1, len(currRep)):
                    outStr += currRep[j] + " "
                outStr = outStr[:-1] + "\n"
                outFd.write(outStr)
                lookCompl = True
                print ("Write Rep")
            elif ( (len(fields) == 5) & (fields[3] == 'exchange:') & (fields[4] == '1') ):
                tmpVal = currRep[int(fields[1])]
                currRep[int(fields[1])] = currRep[int(fields[2])]
                currRep[int(fields[2])] = tmpVal
                print ("Do exch")

outFd.close()

print ("DONE mk_ran_walk.py", outFile)

