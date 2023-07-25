#!/usr/bin/python

import os, sys
import numpy
from scipy.stats import shapiro
from scipy.stats import anderson

def get_sig_val(critValList, sigValList):
#  Look for 5
    useVal = len(sigValList) -1
    for i in range(len(sigValList)):
        if (sigValList[i] <= 5):
            useVal = i
            break
    return critValList[useVal]

def get_mean_sd(inList, maxArray=11000, quiet=True):
    
    valList = []
    for i in range(len(inList)):
        try:
            val = float(inList[i])
            valList.append(val)
        except:
            continue

    total = 0.0
    for i in range(len(valList)):
        total += valList[i]

    myMean = total/len(valList)

    sqDiff = 0.0
    for i in range(len(valList)):
        sqDiff += ((valList[i] - myMean)**2)

    myVar = sqDiff/(len(valList))

    mySd = myVar**0.5

    valArray = numpy.zeros(len(valList))
    arrayCnt = 0
    for i in range(len(valList)):
        valArray[i] = valList[i]
#    statis, p_val = shapiro(valArray)

    statis, critVal, sigVal = anderson(valArray)
    if (quiet == False):
        print "STATIS ", statis, "CRITIVAL ", critVal[2], "SIG val ", sigVal[2]
    critValue = critVal[2]



    arrayCnt = 0
    if (quiet == False):
        print "Mean, Var, SD ", myMean, myVar, mySd
    if (critValue > statis):
        if (quiet == False):
            print "NORMAL distribution "
        isNorm = True
    else:
        if (quiet == False):
            print "not NORMAL distribution "
        isNorm = False
#
#  Check other distributions
#
    statis, critVal, sigVal = anderson(valArray, dist='expon')
    critValue = get_sig_val(critVal, sigVal)
    if (quiet == False):
        if (critValue > statis):
            print "FOUND DIST expon", statis, sigVal

    statis, critVal, sigVal = anderson(valArray, dist='logistic')
    critValue = get_sig_val(critVal, sigVal)
    if (quiet == False):
        if (critValue > statis):
            print "FOUND DIST logistic", statis, sigVal

    statis, critVal, sigVal = anderson(valArray, dist='gumbel')
    critValue = get_sig_val(critVal, sigVal)
    if (quiet == False):
        if (critValue > statis):
            print "FOUND DIST gumbel", statis, sigVal

    statis, critVal, sigVal = anderson(valArray, dist='gumbel_l')
    critValue = get_sig_val(critVal, sigVal)
    if (quiet == False):
        if (critValue > statis):
            print "FOUND DIST gumbel_l", statis, sigVal

    statis, critVal, sigVal = anderson(valArray, dist='gumbel_r')
    critValue = get_sig_val(critVal, sigVal)
    if (quiet == False):
        if (critValue > statis):
            print "FOUND DIST gumbel_r", statis, sigVal

    statis, critVal, sigVal = anderson(valArray, dist='extreme1')
    critValue = get_sig_val(critVal, sigVal)
    if (quiet == False):
        if (critValue > statis):
            print "FOUND DIST extreme1", statis, sigVal



    return myMean, myVar, mySd, isNorm

