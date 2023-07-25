
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

def big_power(x_val, y_val):
#
    if (x_val == 0):
        return(0.0, 0.0)
    if (y_val == 0):
        return(1.0, 0.0)
    elif (y_val == 1):
        out_log_10 = math.log10(x_val)
    else:
        out_log_10 = y_val / math.log(10,x_val)

    out_10_val = int(out_log_10)
    rem_val = out_log_10 - out_10_val
    out_val = math.pow(10,rem_val)

    return(out_val, out_10_val)

#
#   Adds a list of tupals number = [0] * (10 to the [1]power)
#
def big_add(add_list):
#
#  sort list smallest to largest.
#
    good = 0
    cnt = 0
    while (good == 0):
        try:
            tot = float(add_list[cnt][0])
            tot_p10 = float(add_list[cnt][1])
            good = 1
            cnt += 1
        except:
            cnt += 1
    for j in range(cnt,len(add_list)):
        try:
            curr_val = float(add_list[j][0])
            curr_p10 = float(add_list[j][1])
        except:
            continue
        if (curr_p10 > tot_p10):
            diff = curr_p10 - tot_p10
            tot = tot / math.pow(10,diff)
            tot_p10 = curr_p10
        elif (tot_p10 > curr_p10):
            diff = tot_p10 - curr_p10
            curr_val = curr_val / math.pow(10,diff)
        tot += curr_val
        while (tot > 10):
            tot = tot / 10.0
            tot_p10 += 1
    return(tot, tot_p10)
