#! /usr/bin/python

import os,sys
import math
MIN_ARGS=7
if (len(sys.argv) < MIN_ARGS):
    print "NOT ENOUGH argv:\n   reptemp <low temp> <high temp> <num rep> <outFile>"
    exit

lowTemp = float(sys.argv[1])
highTemp = float(sys.argv[2])
lowSolTemp = float(sys.argv[3])
highSolTemp = float(sys.argv[4])
numRep = int(sys.argv[5])
outName = sys.argv[6]

lowLog =  math.log(lowTemp)
highLog =  math.log(highTemp)
lowSolLog =  math.log(lowSolTemp)
highSolLog =  math.log(highSolTemp)
print "low log", lowLog, 

logDiff = (highLog - lowLog) / (numRep-1)
logSolDiff = (highSolLog - lowSolLog) / (numRep-1)

outFp = open(outName,"w")
out2Fp = open("RepTemp.dat","w")
out3Fp = open("ScaleFactor.dat","w")
out4Fp = open("ScaleFactor2.dat","w")
print "OPEN FILE"
for i in range (0, numRep):
    logVal = lowLog + (logDiff * i)
    logSolVal = lowSolLog + (logSolDiff * i)
    valTemp = math.exp(logVal)
    valSolTemp = math.exp(logSolVal)
    intValTemp = int(valTemp+0.5)
    intSolValTemp = int(valSolTemp+0.5)
    valTemp2 = float(intValTemp)
    valSolTemp2 = float(intSolValTemp)
    valScale = valTemp2 / valSolTemp2
    valScale2 =  valSolTemp2 / valTemp2
    print "Temp ", i, "    ", valTemp, intValTemp, valTemp2/lowTemp
    scaleStr = "%.15f" % valScale
    scaleStr2 = "%.15f" % valScale2
    outStr = str(intValTemp) + " " + scaleStr + " " + scaleStr2 + "\n"
    outStrScale = scaleStr + "\n"
    outStrScale2 = scaleStr2 + "\n"
    outFp.write(outStr)
    outStr = str(intValTemp) + "\n"
    out2Fp.write(outStr)
    out3Fp.write(outStrScale)
    out4Fp.write(outStrScale2)
outFp.close()    
out2Fp.close()    
out3Fp.close()    
out4Fp.close()    
