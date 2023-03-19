#! /usr/bin/python3

import os,sys

MIN_ARGS=11
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   ")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file

FORTRAN_VERSION="9"


first_tr = int(sys.argv[1])
num_tr = int(sys.argv[2])
first_rep = int(sys.argv[3])
num_rep = int(sys.argv[4])
first_str = int(sys.argv[5])
num_str = int(sys.argv[6])
size_run = int(sys.argv[7])
intermPath = sys.argv[8]
first_run = int(sys.argv[9])
inFile = sys.argv[10]

num_runs = int(num_str / size_run)
print ("NUM RUNS ", num_runs, num_str, num_str, size_run, (num_str / size_run))
print ("FIRST RUN", first_run )

inList = read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))
outList = []
tmpOutFile = "interm_tmp.f90"
cmd = "rm " + tmpOutFile + " 2>/dev/null"
os.system(cmd)

outFp = open(tmpOutFile, "w")
for i in range(len(inList)):
    line = inList[i]
    outStr = line + "\n"
    offset = line.find("$SIZE_RUN")
    if (offset >= 0):
        outStr = line[:offset] + str(size_run) + line[offset+9:] + "\n"
    else:
        offset = line.find("$PATH")
        if (offset >= 0):
            outStr = line[:offset] + str(intermPath) + line[offset+5:] + "\n"
    outFp.write(outStr)

outFp.close()
cmd = "gfortran-" + FORTRAN_VERSION + " interm_tmp.f90 -g -o interm_tmp.exe"
print ("CMD ", cmd)
errCode = os.system(cmd)
if (errCode):
    print ("ERROR in compile ")
    exit(errCode)


for tr in range(first_tr, first_tr+num_tr):
    print ("TR ", tr)
    for rep in range(first_rep, first_rep+num_rep):
        print ("REP ", rep)
        run_cnt = 0
        print ("DO run ", first_run, num_runs, first_run+num_runs)
        for run in range(first_run, first_run+num_runs):
            print ("RUN ", run, num_str, size_run)
            currStr = first_str + (run_cnt * size_run)
            print ("RUN ", run, currStr)
            outStr = str(tr) + " " + str(rep) + " " + str(run) + " " + str(currStr) + ' "' + '/"\n'
            print ("CMD ", outStr[:-1])
            outFp = open("tmp_run", "w")
            outFp.write(outStr)
            outFp.close()
            cmd = "./interm_tmp.exe < tmp_run"
            print ("DO CMD ", cmd)
            run_cnt += 1
            os.system(cmd)
