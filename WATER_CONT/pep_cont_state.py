#! /usr/bin/python3

import os,sys
MIN_ARGS=5
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   pep_cont_state.py  <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file
import math

clustStateFile = sys.argv[1]
contFileUp = sys.argv[2]
contFileLo = sys.argv[3]
outFile = sys.argv[4]

NUM_BUCKETS = 71
MAX_STATE = 5

cont_state = []
cont_state_cnt = []
for i in range(0,MAX_STATE+1):
    cont_state_cnt.append(0.0)
    cont_state.append([])
    for j in range(0,NUM_BUCKETS):
        cont_state[i].append(0.0)



clustStateList = read_file.read_file_list(clustStateFile)
print ("CLUST ", clustStateFile, len(clustStateList))

clust_state_up = []
clust_state_lo = []
for i in range(len(clustStateList)):
    line = read_file.comp_space(clustStateList[i])
    fields = line.split(' ')
    clust_state_up.append(int(fields[0]))
    clust_state_lo.append(int(fields[1]))

contListUp = read_file.read_file_list(contFileUp)
print ("COM ", contFileUp, len(contListUp))
contListLo = read_file.read_file_list(contFileLo)
print ("COM ", contFileLo, len(contListLo))

for i in range(len(contListUp)):
    currState = clust_state_up[i]
    line = read_file.comp_space(contListUp[i])
    fields = line.split(' ')
    for j in range(len(fields)):
        cont_state[currState][j] += int(fields[j])
        if ((j < 30) & (int(fields[j]) > 0)):
            print("UP NEG VAL ", j, fields[j], "Struct ", i, "state ", currState)
    cont_state_cnt[currState] += 1

for i in range(len(contListLo)):
    currState = clust_state_lo[i]
    line = read_file.comp_space(contListLo[i])
    fields = line.split(' ')
    for j in range(len(fields)):
        cont_state[currState][j] += int(fields[j])
        if ((j < 30) & (int(fields[j]) > 0)):
            print("LO NEG VAL ", j, fields[j], "Struct ", i, "state ", currState)
    cont_state_cnt[currState] += 1

for i in range(0, MAX_STATE):
    for j in range(0, len(cont_state[0])):
        if (cont_state_cnt[i] > 0.0):
            cont_state[i][j] = cont_state[i][j] / cont_state_cnt[i]



fd = open(outFile,"w")
depth = -35
for i in range(len(cont_state[0])):
    outStr = str(i + depth) + " "
    for j in range(0, MAX_STATE+1):
        outStr += str(round(cont_state[j][i],5)) + " "
    outStr = outStr[:-1] + "\n"
    fd.write(outStr)
fd.close()



print ("DONE get_clust_state.py ", outFile)
