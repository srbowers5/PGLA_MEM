
import math

def factorial(val):
    tot_val = 0.0
    for j in range (1,val+1):
        float_val = float(j)
        tot_val += math.log10(float_val)
    out_10_val = int(tot_val)
    rem_val = tot_val - out_10_val
    out_val = math.pow(10,rem_val)
    return(out_val, out_10_val)
