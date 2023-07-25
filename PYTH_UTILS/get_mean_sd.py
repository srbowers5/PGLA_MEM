#!/usr/bin/python

import os, sys

def get_mean_sd(inList):
    valArray = []
    for i in range(len(inList)):
        try:
            val = float(inList[i])
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

#    print ("Mean, Var, SD ", myMean, myVar, mySd)
    return myMean, myVar, mySd

