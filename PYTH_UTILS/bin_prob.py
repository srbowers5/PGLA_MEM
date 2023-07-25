import math
import bmath

def binomial_prob_ge(success, tot, success_prob):

    tot_prob_10 = 0.0
    tot_prob_rem = 0.0

    tot_prob_rem, tot_prob_10 = binomial_prob(tot, tot, success_prob)

    for j in range (tot-1, success-1, -1):
        prob_rem, prob_10 = binomial_prob(j, tot, success_prob)
# convert to same power of 10 (use higher power)
        prob_10_diff = (tot_prob_10 - prob_10)
        if (prob_10_diff > 307):
            continue
        elif (prob_10_diff < -307):
            tot_prob_10 = prob_10
            tot_prob_rem = prob_rem
        elif (prob_10_diff > 0):
            prob_rem = prob_rem / math.pow(10,prob_10_diff)
        elif (prob_10_diff < 0):
            tot_prob_10 = prob_10
            tot_prob_rem = tot_prob_rem * math.pow(10,prob_10_diff)
            
        tot_prob_rem += prob_rem
        while (tot_prob_rem > 10):
            tot_prob_10 += 1
            tot_prob_rem = tot_prob_rem / 10
    return(tot_prob_rem, tot_prob_10)

def binomial_prob_le(success, tot, success_prob):

    tot_prob_10 = 0.0
    tot_prob_rem = 0.0

    tot_prob_rem, tot_prob_10 = binomial_prob(0, tot, success_prob)

    for j in range (1, success+1):
        prob_rem, prob_10 = binomial_prob(j, tot, success_prob)
# convert to same power of 10 (use higher power)
        prob_10_diff = (tot_prob_10 - prob_10)
        if (prob_10_diff > 307):
            continue
        elif (prob_10_diff < -307):
            tot_prob_10 = prob_10
            tot_prob_rem = prob_rem
        elif (prob_10_diff > 0):
            prob_rem = prob_rem / math.pow(10,prob_10_diff)
        elif (prob_10_diff < 0):
            tot_prob_10 = prob_10
            tot_prob_rem = tot_prob_rem * math.pow(10,prob_10_diff)

        tot_prob_rem += prob_rem
        while (tot_prob_rem > 10):
            tot_prob_10 += 1
            tot_prob_rem = tot_prob_rem / 10
    return(tot_prob_rem, tot_prob_10)


def binomial_prob(success, tot, success_prob):
    rem1, fact_10_1 = bmath.factorial(tot-success)
    rem2, fact_10_2 = bmath.factorial(success)
    rem3, fact_10_3 = bmath.factorial(tot)

    fact_10 = fact_10_3 - (fact_10_1 + fact_10_2)
    rem =  rem3 / (rem1 * rem2)

    rem_p1, fact_10_p1 = bmath.big_power(success_prob, success)
    rem_p2, fact_10_p2 = bmath.big_power(1.0 - success_prob, tot - success)
    
    rem = rem * rem_p1 * rem_p2
    fact_10 = fact_10 + fact_10_p1 + fact_10_p2

    while (rem > 10):
        fact_10 += 1
        rem = rem / 10
    while (rem < 1):
        fact_10 -= 1
        rem = rem * 10

    return(rem, fact_10)
