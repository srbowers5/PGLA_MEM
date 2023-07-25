#! /usr/bin/python3

import os,sys
MIN_ARGS=6
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   G_2D_from_prob.py <input> <output> <temp> <num_col> <lib paths>")
    exit()

if (len(sys.argv) > MIN_ARGS):
    for i in (MIN_ARGS, len(sys.argv)-1):
        sys.path.append(sys.argv[i])
        print ("ADD PATH", i, sys.argv[i])

import read_file
import math

inFile = sys.argv[1]
outFile = sys.argv[2]
tempVal = float(sys.argv[3])
num_col = int(sys.argv[4])
sig_digits = int(sys.argv[5])

outFd = open(outFile, "w")

inputList = read_file.read_file_list(inFile)

RVal = 1.987    # cal/mol/K
# Find Max
maxVal = 0.0
sumVal = 0.0
floatFields = []
energyArray = []


maxVal = 0.0
sumVal = 0.0
minVal = 1.0
for i in range(len(inputList)):
#   Field1 = label, field2 = value
    fields = inputList[i].split()
    if (len(fields) != num_col):
        print ("short ", len(fields), num_col)
        continue
    for j in range(0, len(fields)):
        newVal = float(fields[j])
        if (newVal > maxVal):
            maxVal = newVal
        if ((newVal < minVal) & (newVal != 0.0)):
            minVal = newVal
        sumVal += newVal
minVal = minVal /10.0
print ("MAX/MIN", maxVal, minVal, sumVal)

cnt = 0
maxFE = 0.0
for i in range(len(inputList)):
    fields = inputList[i].split()
    if (len(fields) != num_col):
        continue

    outStr = ""
    for j in range(0, len(fields)):
        val = fields[j]
        floatFields = []
        val = float(fields[j])
        # No zero prob, so set to lowest/10.
        if (val < minVal):
            val = minVal
    
        val = math.log(val / maxVal)
        if (val != 0.0):
            val = -(val * RVal * tempVal) / 1000
        if (val > maxFE):
            maxFE = val

#  (f'{value:.6f}')
        if (sig_digits == 2):
            outStr += (f'{val:.2f}') + " "
        elif (sig_digits == 3):
            outStr += (f'{val:.3f}') + " "
        elif (sig_digits == 4):
            outStr += (f'{val:.4f}') + " "
        elif (sig_digits == 5):
            outStr += (f'{val:.5f}') + " "
        elif (sig_digits == 6):
            outStr += (f'{val:.6f}') + " "
        elif (sig_digits == 7):
            outStr += (f'{val:.7f}') + " "
        elif (sig_digits == 8):
            outStr += (f'{val:.8f}') + " "
        elif (sig_digits == 9):
            outStr += (f'{val:.9f}') + " "
        elif (sig_digits == 10):
            outStr += (f'{val:.10f}') + " "
        else:
            outStr += str(round(val,sig_digits)) + " "
    outStr = outStr[:-1] + "\n"
    outFd.write(outStr)

print ("SUM ", sumVal, " MAX ", maxVal, " MIN ", minVal, "MAX FE ", maxFE)


outFd.close()


print ("DONE G_2D_from_prob.py", outFile)

