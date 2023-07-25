#!/usr/bin/python

import math
#
#  First arg is probibility at 0.
def get_deltaG(zero_prob, other_prob, temper_val):
    RVAL = 1.987    # cal/mol/K
    val = math.log(other_prob / zero_prob)
    val = -(val * RVAL * 310) / 1000
    return val  # Kcal/mole
