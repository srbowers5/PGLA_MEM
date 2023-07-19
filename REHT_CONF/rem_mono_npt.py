#! /usr/bin/python3
#
#
#  Check if already running and quit if so
#
#  Check if Runstatus exists and quit if so.
#
#  Remove junk
#
#  if first step, Copy input structs form inputSTR
#
#  For each rep create namd template file
#
#  For each rep create namd_swap template file
#
#  For each step.
#      Write to Runstatus
#      for each rep
#          Create namd files
#          Run namd
#          Run namd swap
#
#  Wait until all are done
#  For each rep
#      Get energy
#      Get swap Energy
#      Get XSC data into RepBox.dat
#  call exchange_npt.exe <status.dat > exchange.out
#  If error
#    exit
#  Handle exchanges

import os,sys
import subprocess
import threading
import time

MIN_ARGS=7
if (len(sys.argv) < MIN_ARGS):
    print ("Args", MIN_ARGS, len(sys.argv))
    print ("NOT ENOUGH argv:\n   ./rem_mono_npt.py <tr> <runNum> <numRep> <repFirst> <repLast> <lib paths>")
    exit
for i in range(MIN_ARGS,len(sys.argv)):
    print ("add path", sys.argv[i])
    sys.path.append(sys.argv[i])

import read_file
import run_namd

tr = sys.argv[1]
runNum = sys.argv[2]
numRep = int(sys.argv[3])
repFirst = int(sys.argv[4])
stepMax = int(sys.argv[5])
maxCpu = int(sys.argv[6])

TEMP_MIN = "310"
cmd = "ln -s /run/user/1000 myTmpFs 2>/dev/null"
os.system(cmd)


cmd = "ps ax  | grep -v grep | grep python | grep rem_mono_npt | wc -l"
result = subprocess.check_output(cmd, shell=True)

if (int(result) > 1):
    print ("ALREADY running ", result)
    exit(99)

cmd = "rm ERROR_exchange 2>/dev/null"
os.system(cmd)
cmd = "rm ERROR_namd 2>/dev/null"
os.system(cmd)




def replace_str(input_string, prev_string, new_string):
    off = input_string.find(prev_string)
    if (off >= 0):
        return_string = input_string[:off] + new_string + input_string[len(prev_string)+off:]
        return return_string
    else:
        return input_string
    
    


no_namd_allowed = False

if (no_namd_allowed):
    cmd = "ps ax  | grep -v grep | grep namd | wc -l"
    result = subprocess.check_output(cmd, shell=True)

    if (int(result) > 0):
        print ("NAMD ALREADY running ", result)
        exit(98)


runFile = "RunStatus" + runNum
if (os.path.isfile(runFile)):
    print (runFile, " exists. Change batch or remove file ")
    exit(97)

if (runNum == '1'):
    cmd = "rm forces*.dat energy*.dat xsc*dat 2>/dev/null"
    os.system(cmd)
cmd = "rm pgl_mem*.out pgl_mem_template_*.namd pgl_mem_template_*_swap.namd RepEnergySwap.dat RepEnergy.dat 2>/dev/null"
os.system(cmd)

if (repFirst == 1):
    print ("INITIAL")
    for rep in range(1,numRep+1):
        cmd = "rm pgl_mem" + str(rep) + ".coor " + \
          "pgl_mem" + str(rep) + ".vel " + \
          "pgl_mem" + str(rep) + ".xsc 2>/dev/null"
        os.system(cmd)
        cmd = "rm pgl_mem" + str(rep) + "_i.coor " + \
          "pgl_mem" + str(rep) + "_i.vel " + \
          "pgl_mem" + str(rep) + "_i.xsc  2>/dev/null"
        os.system(cmd)
        cmd = "rm pgl_mem" + str(rep) + "_swap.coor " + \
          "pgl_mem" + str(rep) + "_swap.vel " + \
          "pgl_mem" + str(rep) + "_swap.xsc  2>/dev/null"
        os.system(cmd)

#        cmd = "cp inputSTR/memi_equil1_" + str(rep) + ".coor pgl22_mem" + str(rep) + "_0.coor"
#        os.system(cmd)
#        cmd = "cp inputSTR/memi_equil1_" + str(rep) + ".vel pgl22_mem" + str(rep) + "_0.vel"
###        os.system(cmd)
#        cmd = "cp inputSTR/memi_equil1_" + str(rep) + ".xsc pgl22_mem" + str(rep) + "_0.xsc"
#        os.system(cmd)


        cmd = "cp pgl22_mem" + str(rep) + "_0.coor pgl_mem" + str(rep) + "_i.coor"
        os.system(cmd)

        cmd = "cp pgl22_mem" + str(rep) + "_0.vel pgl_mem" + str(rep) + "_i.vel"
        os.system(cmd)

        cmd = "cp pgl22_mem" + str(rep) + "_0.xsc pgl_mem" + str(rep) + "_i.xsc"
        os.system(cmd)

temperList = read_file.read_file_list("RepTemp.dat")
scaleList = read_file.read_file_list("ScaleFactor.dat")
scheduleList = read_file.read_file_list("schedule.dat")


templateList = read_file.read_file_list("pgl_mem_template.namd")

tempFileList = []
tempSwapFileList = []
for rep in range (0,numRep):
    tempFileList.append([])
    tempSwapFileList.append([])


for rep in range (1,numRep+1):
    filename = "pgl_mem_template_" + str(rep) + ".namd"
    outFd = open(filename,"w")
    filename = "pgl_mem_template_" + str(rep) + "_swap.namd"
    outFdSwap = open(filename,"w")
    for i in range(len(templateList)):
        outStr = swapStr = templateList[i]


        outStr = replace_str(outStr, "FinalFile", ("pgl_mem" + str(rep)))
        swapStr = replace_str(swapStr, "FinalFile", ("pgl_mem" + str(rep) + "_swap"))
        outStr = replace_str(outStr, "RepTemp", temperList[rep-1])
        outStr = replace_str(outStr, "RunLength", "2000")
        swapStr = replace_str(swapStr, "RunLength", "0")
        outStr = replace_str(outStr,  "RestartFile", "pgl_mem" + str(rep) + "_i")
        swapStr = replace_str(swapStr, "RestartFile", "pgl_mem" + str(rep))
        outStr = replace_str(outStr, "OutputFrequency", "2000")
        swapStr = replace_str(swapStr, "OutputFrequency", "2000")
        outStr = replace_str(outStr, "StructureFile", "pgl_mem" + str(rep))
        outStr = replace_str(outStr, "SF", scaleList[rep-1])
        outStr = replace_str(outStr, "ConstraintsFile", "constraints")

        tempFileList[rep-1].append(outStr)
        tempSwapFileList[rep-1].append(swapStr)
        outStr = outStr + "\n"
        swapStr = swapStr + "\n"
        outFd.write(outStr)
        outFdSwap.write(swapStr)

    outFd.close()
    outFdSwap.close()

#
#  Start the namd threads
#  We will only start 12 at a time. When one finishes start a new one.
# 
for stepNum in range(repFirst, stepMax+1):
    outStr = "Replica step: " + str(stepNum) + " ' >>  RunStatus" + runNum
    cmd = "echo '" + outStr
    os.system(cmd)

    threadVals = []
    nextThread = 1
    while (nextThread <= numRep):
         newThread = threading.Thread(target=run_namd.run_namd, args=(tr, runNum, nextThread,stepNum, numRep, TEMP_MIN, scheduleList,temperList,scaleList,tempFileList,tempSwapFileList ))
         newThread.start()
         threadVals.append(newThread)
         nextThread += 1
         numActive = threading.activeCount()
         while (numActive > maxCpu):
             time.sleep(1)
             numActive = threading.activeCount()
#
#  All are started. Now check if all are finished
#
    for i in range(len(threadVals)):
        threadVals[i].join()
   
    if (os.path.isfile("./ERROR_namd") == True):
        print ("ERROR File in namd")
        break

#
    outStr = "Replica step: " + str(stepNum) + " completed ' >>  RunStatus" + runNum
    cmd = "echo '" + outStr
    os.system(cmd)


    line = read_file.comp_space(scheduleList[stepNum-1])
    fields = line.split(' ')
    lVal = fields[0]
    outFd = open("status.dat", "w")
    outStr = str(stepNum) + " " + tr + " " + str(numRep) + " " + TEMP_MIN + " " + lVal + "\n"
    outFd.write(outStr)
    outFd.flush()
    os.fsync(outFd)
    outFd.close()


    cmd = "./do_exchange.sh " + str(runNum) + " " + str(numRep)
    print ("CMD ", cmd)
    os.system(cmd)

    if (os.path.isfile("./ERROR_exchange") == True):
        print ("ERROR File in exchange")
        break


    


print ("DONE")
