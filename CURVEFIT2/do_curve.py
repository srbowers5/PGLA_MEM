#! /usr/bin/python3


import os,sys
MIN_ARGS=9
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   cnt_gene.py <input> <output>  <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

for i in range (len(sys.argv)):
    print ("AFG", i, sys.argv[i])

sys.stdout.flush()

import read_file
import math
import numpy
from scipy.optimize import curve_fit


inFile = sys.argv[1]
outFile = sys.argv[2]
numExp = int(sys.argv[3])
ratio = float(sys.argv[4])
ratio2 = float(sys.argv[5])
ratioa0_1 = float(sys.argv[6])
ratioa1_2 = float(sys.argv[7])
valCol = int(sys.argv[8])

if (ratioa0_1 < 0):
    ratioa0_1 = -(1/ratioa0_1)

if (ratioa1_2 < 0):
    ratioa1_2 = -(1/ratioa1_2)

ratioa=0
if (ratioa == 0):
    ratio2a = 1
    ratio3a = 0
elif (ratioa < 0):
    ratioa = ratioa
    ratioa = -ratioa
    ratio3a = 1.0/(ratioa+1)
    ratio2a = (ratioa)/(ratioa+1)
else:
    ratio2a = 1.0/(ratioa+1)
    ratio3a = (ratioa)/(ratioa+1)

#
#
#a1 = ((1-a0) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2))))
#a2 = (((1-a0) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) / ratioa0_1)
#a3 = (((1-a0) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) / (ratioa0_1*ratioa1_2))

# y = a0 + (a1 * exp(-t/t1)) + (a2 * exp(-t/t2)) + (a3 * exp(-t/t3))

def objective3(t, a0, t1, t2, t3):
    return a0 + \
      (((1-a0) /  (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) * (numpy.exp(-(t/t1)))) + \
      ((((1-a0) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) / ratioa0_1) * (numpy.exp(-(t/t2))) ) + \
      ((((1-a0) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) / (ratioa0_1*ratioa1_2)) * (numpy.exp(-(t/t2))))


def objective2(t, a0, t1, t2):
#    return a0 + (a1 * numpy.exp(-(t/t1))) + ((1 - (a0+a1)) * numpy.exp(-(t/(t1/ratio))))
    return a0 + (((1-a0) * (ratioa0_1 / (ratioa0_1 + 1))) * numpy.exp(-t/t1)) + \
      (((1-a0) / (ratioa0_1 + 1)) * numpy.exp(-t/t2))

def objective1(t, a0, t1):
    return a0 +( (1.0-a0) * numpy.exp(-(t/t1)))

def objective0(t, a0, a1):
    return a0 + (a1 * t)


#
#   main()
inList = read_file.read_file_list(inFile)
print ("READ ", inFile, len(inList))

t=[]
y=[]
for i in range(len(inList)):
    line = read_file.comp_space(inList[i])
    fields = line.split(' ')
    t.append(fields[0])
    y.append(fields[valCol])

#for i in range(len(t)):
    #print ("OFF", i, t[i], y[i])
tVec = numpy.array(t,dtype=numpy.float128)
yVec = numpy.array(y,dtype=numpy.float128)
#print ("TVEC", tVec)
#print ("YVEC", yVec)


if (numExp == 0):
    fpOut = open(outFile,"w")
    a2Val=0
    a3Val=0
    a4Val=0
    t1Val=0
    t2Val=0
    t3Val=0
    t4Val=0
    iniA0 = 0.25
    while(iniA0 < 0.99999):
        iniA1 = 1 - iniA0
        maxA0 = iniA0 + 0.05
        if (maxA0 > 1.0):
            maxA0 = 1.0
        minA0 = iniA0 - 0.05
        if (minA0 < 0):
            minA0 = 0.01
        popt, _ = curve_fit(objective0, tVec, yVec, p0=[iniA0, -1], bounds=( [minA0,-1000], [maxA0, 1] ))
        a0Val, a1Val = popt
        iniA0 += 0.10


        outVals = []
        for i in range(len(tVec)):
            inVal = tVec[i]
            val = a0Val * (a1Val * inVal)
            outVals.append(val)

        sqSum = 0.0
        for i in range(len(yVec)):
            val = (yVec[i] - outVals[i]) ** 2
            sqSum += val
        sqSum = sqSum/len(tVec)
#        print("SQ SUM ", sqSum)
    
        ratStr = str(round(ratioa,6))
        outStr = str(ratio) + "\t" + str(ratio2) + "\t" + ratStr.ljust(7,'0') + " " + \
            str(round(abs(a0Val),6)) + "\t" + \
            str(round(abs(a1Val),6)) + "\t" + \
            str(round(abs(a2Val),6)) + "\t" + \
            str(round(abs(a3Val),6)) + "\t" + \
            str(round(abs(a4Val),6)) + "\t" + \
            str(round(t1Val,6)) + "\t" + \
            str(round(t2Val,6)) + "\t" + \
            str(round(t3Val,6)) + "\t" + \
            str(round(t4Val,6)) + "\t" + \
            str(round(sqSum,6)) + "\n"
        fpOut.write(outStr)
        iniA0 += 0.25
    fpOut.close()

    exit(0)



fpOut = open(outFile,"w")
iniA0 = 0.34
while(iniA0 < 1.0):
    iniA1 = 1 - iniA0

    if (numExp == 1):
        popt, _ = curve_fit(objective1, tVec, yVec, p0=[iniA0, 20], bounds=( [0,0], [1, numpy.inf] ))
        a0Val, t1Val = popt
        a1Val = (1-a0Val)
        a2Val = 0
        t2Val = 1
        a3Val = 0
        t3Val = 1
        a4Val = 0
        t4Val = 1
    elif (numExp == 2):
        popt, _ = curve_fit(objective2, tVec, yVec, p0=[0.5, 20, 20*ratio], bounds=( [0.0, 0.0, 0.0], [1.0, numpy.inf, numpy.inf] ))
#        popt, _ = curve_fit(objective2, tVec, yVec)
#        a0Val, a1Val,t1Val = popt
        a0Val, t1Val, t2Val = popt


        a1Val = ((1-a0Val) * (ratioa0_1 / (ratioa0_1 + 1)))
        a2Val = ((1-a0Val) * (1 /(ratioa0_1 + 1)))

        a3Val = 0
        t3Val = 1
        a4Val = 0
        t4Val = 1


    elif (numExp == 3):
        popt, _ = curve_fit(objective3, tVec, yVec, p0=[0.5, 20, 20*ratio, 20*ratio2], bounds=( [0.0, 0.0, 0.0, 0.0], [1.0, numpy.inf, numpy.inf, numpy.inf] )) 
        a0Val, t1Val, t2Val, t3Val = popt

        a1Val = ((1-a0Val) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2))))
        a2Val = (((1-a0Val) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) / ratioa0_1)
        a3Val = (((1-a0Val) / (1 + (1.0/ratioa0_1) + (1.0/(ratioa0_1 * ratioa1_2)))) / (ratioa0_1*ratioa1_2))

        a4Val = 0
        t4Val = 1
    else:
        popt, _ = curve_fit(objective4, tVec, yVec)
        a0Val, a1Val, t1Val, a2Val, t2Val, a3Val, t3Val, a4Val, t4Val = popt
    
    
    #print ("VALS", popt)
    #print ("VALS2", a0Val, a1Val, a2Val, t1Val)
    outVals = []
    for i in range(len(t)):
        inVal = tVec[i]
        val =  a0Val + (a1Val * numpy.exp(-inVal/t1Val)) + (a2Val * numpy.exp(-inVal/t2Val)) + (a3Val * numpy.exp(-inVal/t3Val)) + (a4Val * numpy.exp(-inVal/t4Val))
        outVals.append(val)
    
    
    #for i in range(len(t)):
        #print( "FINAL ", tVec[i], "   ", round(yVec[i],4), "    ", round(outVals[i],4) )
    
    sqSum = 0.0
    for i in range(len(t)):
        val = (yVec[i] - outVals[i]) ** 2
        sqSum += val
        #print ("SqSuM, VAL", val, sqSum)
    
    
    
    sqSum = sqSum/len(t)
#    print("SQ SUM ", sqSum)
    
    ratStr = str(round(ratioa0_1,6))
    rat2Str = str(round(ratioa1_2,6))
    outStr = str(ratio) + "\t" + str(ratio2) + "\t" + \
        str(ratStr.ljust(7,'0')) + "\t" + \
        str(rat2Str.ljust(7,'0')) + "\t" + \
        str(round(abs(a0Val),6)) + "\t" + \
        str(round(abs(a1Val),6)) + "\t" + \
        str(round(abs(a2Val),6)) + "\t" + \
        str(round(abs(a3Val),6)) + "\t" + \
        str(round(abs(a4Val),6)) + "\t" + \
        str(round(t1Val,6)) + "\t" + \
        str(round(t2Val,6)) + "\t" + \
        str(round(t3Val,6)) + "\t" + \
        str(round(t4Val,6)) + "\t" + \
        str(round(sqSum,6)) + "\n"
    fpOut.write(outStr)
    iniA0 += 0.34
fpOut.close()


print ("DONE ")
