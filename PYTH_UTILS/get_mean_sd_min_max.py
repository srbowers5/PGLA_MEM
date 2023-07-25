#!/usr/bin/python

import os, sys

def get_mean_sd(inList):
    minVal = 999999999999.0
    maxVal = -999999999999.0

    valArray = []
    for i in range(len(inList)):
        try:
            val = float(inList[i])
            if (val < minVal):
                minVal = val
            if (val > maxVal):
                maxVal = val
            valArray.append(val)
        except:
            continue

    total = 0.0
    for i in range(len(valArray)):
        total += valArray[i]

    myMean = total/len(valArray)

    sqDiff = 0.0
    for i in range(len(valArray)):
        sqDiff += ((valArray[i] - myMean)**2)

    myVar = sqDiff/(len(valArray))

    mySd = myVar**0.5

    print "Mean, Var, SD, min, max ", myMean, myVar, mySd, minVal, maxVal
    return myMean, myVar, mySd, minVal, maxVal

