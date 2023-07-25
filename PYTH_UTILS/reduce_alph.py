
import numpy

res_map = { 'A' : 0, 'V' : 0, 'F' : 0, 'P' : 0, 'M' : 0, 
'I' : 0, 'L' : 0, 'D' : 1, 'E' : 1, 'K' : 1, 'R' : 1, 
'S' : 2, 'T' : 2, 'Y' : 2, 'H' : 2, 'C' : 2, 'N' : 2, 
'Q' : 2, 'W' : 2, 'G' : 2 }

def get_simp_type(simp_string):
    if (len(simp_string) != 4):
        return(-1)

    type_cnt_array = numpy.zeros( 3, dtype=int)
    for j in range(0,4):
        my_type = res_map.get(simp_string[j],4)
        if (my_type == 4):
            return(-1)
        type_cnt_array[my_type] += 1

    val = type_cnt_array[0] + (type_cnt_array[1] * 5) + \
      (type_cnt_array[2] * 25) 

    return(val)

def get_simp_string(val):
    num_polar = val / 25
    num_charged = (val - (num_polar * 25)) / 5
    num_phobic = (val - ((num_polar * 25) + (num_charged * 5)))
    out_string = ""
    for j in range(0,num_phobic):
        out_string += "H"
    for j in range(0,num_charged):
        out_string += "C"
    for j in range(0,num_polar):
        out_string += "P"

    return(out_string)
