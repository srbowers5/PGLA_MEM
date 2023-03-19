#! /usr/bin/python3

import os,sys
MIN_ARGS=8
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   do_interm.py <tr> <num of rep> <first str> <sec str> <prefix> <inputDir> <outDir>")
    print ("   <prefix> <inputDir> <outputFile>")
    exit()
if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])



traFirst = int(sys.argv[1])
traLast = int(sys.argv[2])
rangeVal = sys.argv[3]
DO_1 = int(sys.argv[4])
DO_2 = int(sys.argv[5])
DO_3 = int(sys.argv[6])
DO_4 = int(sys.argv[7])

DO_2=0
import read_file

inFile1 = "hist1_10.f90"
inList1 = read_file.read_file_list(inFile1)
print ("FILE ", inFile1, "len ", len(inList1))

inFile2 = "hist2_3.f90"
inList2 = read_file.read_file_list(inFile2)
print ("FILE ", inFile2, "len ", len(inList2))

inFile3 = "hist3_16.f90"
inList3 = read_file.read_file_list(inFile3)
print ("FILE ", inFile3, "len ", len(inList3))

inFile4 = "hist4_5.f90"
inList4 = read_file.read_file_list(inFile4)
print ("FILE ", inFile4, "len ", len(inList4))



for traNum in range(traFirst, traLast+1):

    if (DO_1 == 1):
        outList = []
        tmpOutFile = "shist1z_aedit05_tmp.f90"
        outFp = open(tmpOutFile, "w")
        for i in range(len(inList1)):
            line = inList1[i]
            outStr = line + "\n"
            offset = line.find("in_tr=1")
            if (offset >= 0):
                outStr = "integer, parameter :: in_tr=" + str(traNum) + "\n"
            offset = line.find("trmin=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmin=" + str(traNum) + "\n"
            offset = line.find("trmax=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmax=" + str(traNum) + "\n"

            offset = line.find("$RANGE")
            if (offset >= 0):
                outStr = line[:offset] + rangeVal + line[offset+6:] + "\n"

            outFp.write(outStr)
        outFp.close()
        cmd = "gfortran shist1z_aedit05_tmp.f90 -o shist1z_aedit05_tmp.exe"
        print ("CMD ", cmd)
    
        os.system(cmd)
        cmd = "./shist1z_aedit05_tmp.exe > 1_" + str(traNum) + ".out"
        print ("CMD ", cmd)
        os.system(cmd)

        outList = []
        tmpOutFile = "hist4_2_tmp.f90"
        outFp = open(tmpOutFile, "w")
        for i in range(len(inList1)):
            line = inList1[i]
            outStr = line + "\n"
            offset = line.find("in_tr=1")
            if (offset >= 0):
                outStr = "integer, parameter :: in_tr=" + str(traNum) + "\n"
            offset = line.find("trmin=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmin=" + str(traNum) + "\n"
            offset = line.find("trmax=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmax=" + str(traNum) + "\n"

            offset = line.find("$RANGE")
            if (offset >= 0):
                outStr = line[:offset] + rangeVal + line[offset+6:] + "\n"

            outFp.write(outStr)


    if (DO_2 == 1):
        outList = []
        tmpOutFile = "shist2edit6_tmp.f90"
        outFp = open(tmpOutFile, "w")
        for i in range(len(inList2)):
            line = inList2[i]
            outStr = line + "\n"
            offset = line.find("in_tr=1")
            if (offset >= 0):
                outStr = "integer, parameter :: in_tr=" + str(traNum) + "\n"
            offset = line.find("trmin=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmin=" + str(traNum) + "\n"
            offset = line.find("trmax=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmax=" + str(traNum) + "\n"
            outFp.write(outStr)
        outFp.close()
        cmd = "gfortran shist2edit6_tmp.f90 -o shist2edit6_tmp.exe"
        print ("CMD ", cmd)

        os.system(cmd)
        cmd = "./shist2edit6_tmp.exe > 2_" + str(traNum) + ".out"
        print ("CMD ", cmd)
        os.system(cmd)

    if (DO_3 == 1):
        outList = []
        tmpOutFile = "shist7edit02_tmp.f90"
        outFp = open(tmpOutFile, "w")
        for i in range(len(inList3)):
            line = inList3[i]
            outStr = line + "\n"
            offset = line.find("in_tr=1")
            if (offset >= 0):
                outStr = "integer, parameter :: in_tr=" + str(traNum) + "\n"
            offset = line.find("trmin=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmin=" + str(traNum) + "\n"
            offset = line.find("trmax=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmax=" + str(traNum) + "\n"
            outFp.write(outStr)
        outFp.close()
        cmd = "gfortran shist7edit02_tmp.f90 -o shist7edit02_tmp.exe"
        print ("CMD ", cmd)
    
        os.system(cmd)
        cmd = "./shist7edit02_tmp.exe > 3_" + str(traNum) + ".out"
        print ("CMD ", cmd)
        os.system(cmd)

    if (DO_4 == 1):

        outList = []
        tmpOutFile = "hist4_5_tmp.f90"
        tmpExeFile = "hist4_5_tmp.exe"
        outFp = open(tmpOutFile, "w")
        for i in range(len(inList4)):
            line = inList4[i]
            outStr = line + "\n"
            offset = line.find("in_tr=1")
            if (offset >= 0):
                outStr = "integer, parameter :: in_tr=" + str(traNum) + "\n"
            offset = line.find("trmin=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmin=" + str(traNum) + "\n"
            offset = line.find("trmax=1")
            if (offset >= 0):
                outStr = "integer, parameter :: trmax=" + str(traNum) + "\n"

            offset = line.find("$RANGE")
            if (offset >= 0):
                outStr = line[:offset] + rangeVal + line[offset+6:] + "\n"

            outFp.write(outStr)
        outFp.close()
        cmd = "gfortran " + tmpOutFile + " -o " + tmpExeFile
        print ("CMD ", cmd)
    
        os.system(cmd)
        cmd = "./" + tmpExeFile + " > 4_" + str(traNum) + ".out"
        print ("CMD ", cmd)
        os.system(cmd)




print ("DONE")

