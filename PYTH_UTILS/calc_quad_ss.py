#include "quad_ss.h"
#include "simp_type.h"

def count_ss(sec_stru):
    num_helix = 0
    num_turn = 0
    num_strand = 0
    num_bridge = 0
    if ((sec_stru[0] != 'H') & \
    (sec_stru[0] != 'E') & \
    (sec_stru[0] != 'T') & \
    (sec_stru[0] != 'E') & \
    (sec_stru[0] != 'B') & \
    (sec_stru[0] != 'G') & \
    (sec_stru[0] != 'I') & \
    (sec_stru[0] != 'S') & \
    (sec_stru[0] != 'c')):
       return(-1, -1, -1, -1)
    for j in range (0,4):
        if (sec_stru[j] == 'H'):
            num_helix +=1
        elif (sec_stru[j] == 'E'):
            num_strand +=1
        elif (sec_stru[j] == 'S'):
            num_strand +=1
        elif (sec_stru[j] == 'T'):
            num_turn +=1
        elif (sec_stru[j] == 'B'):
            num_bridge +=1
    return(num_helix, num_turn, num_strand, num_bridge)

def get_quad_ss(res1, res2, res3, res4, type, sec_stru):

    t4000 = 4
    t3100 = 3
    t2200 = 2
    t2110 = 1
    t1111 = 0
    t_bad = -1

    res_ss_list = []
    res_ss_list.append((int(res1),sec_stru[0]))
    res_ss_list.append((int(res2),sec_stru[1]))
    res_ss_list.append((int(res3),sec_stru[2]))
    res_ss_list.append((int(res4),sec_stru[3]))
    res_ss_list.sort()
    num_helix, num_turn, num_strand, num_bridge = count_ss(sec_stru)
    if (num_helix < 0):
        return("BAD SS")

    if (type == t4000):
        if (num_turn >= 2):
            return("TURN_SS")
        if (num_helix == 4):
            return("FULL_HELIX")
        if (num_helix == 3):
            return("PARTIAL_HELIX")
        if (num_helix == 2):
            if ((res_ss_list[0][1] == 'H') & (res_ss_list[1][1] == 'H')):
                return("PARTIAL_HELIX")
            if ((res_ss_list[2][1] == 'H') & (res_ss_list[3][1] == 'H')):
                return("PARTIAL_HELIX")
        if (num_strand >= 2):
            return("PARTIAL_BETA")

    if (type == t3100):
        if (num_helix == 2):
            if ((res_ss_list[0][0] + 1) == res_ss_list[1][0]):
                if ((res_ss_list[0][1] == 'H') & (res_ss_list[1][1] == 'H')):
                    return("PARTIAL_HELIX")
                if ((res_ss_list[1][1] == 'H') & (res_ss_list[2][1] == 'H')):
                    return("PARTIAL_HELIX")
            else:
                if ((res_ss_list[1][1] == 'H') & (res_ss_list[2][1] == 'H')):
                    return("PARTIAL_HELIX")
                if ((res_ss_list[2][1] == 'H') & (res_ss_list[3][1] == 'H')):
                    return("PARTIAL_HELIX")
        elif (num_helix >= 3):
                    return("PARTIAL_HELIX")
        elif (num_bridge == 2):
            return("BRIDGE_BETA")
        if (num_strand == 2):
            return("PARTIAL_BETA")
        if (num_strand == 3):
            return("PARTIAL_BETA")
        if (num_strand == 4):
            return("FULL_BETA")
        if (num_turn >= 2):
            return("TURN_SS")

    if (type == t2200):
        if (num_helix == 4):
# Helix i, i+1, i+3, i+4 
            if ((res_ss_list[0][0] + 4) == res_ss_list[3][0]):
                return("FULL_HELIX")
            else:
                return("PARTIAL_HELIX")
        if (num_helix >= 2):
            if ((res_ss_list[0][1] == 'H') & (res_ss_list[1][1] == 'H')):
                return("PARTIAL_HELIX")
            if ((res_ss_list[2][1] == 'H') & (res_ss_list[3][1] == 'H')):
                return("PARTIAL_HELIX")
        if (num_bridge == 2):
            return("BRIDGE_BETA")
        if (num_strand == 2):
            return("PARTIAL_BETA")
        if (num_strand == 3):
            return("PARTIAL_BETA")
        if (num_strand == 4):
            return("FULL_BETA")
        if (num_turn >= 2):
            return("TURN_SS")

    if (type == t2110):
        if (num_helix == 4):
            if (((res_ss_list[0][0] + 3) == res_ss_list[1][0]) & \
            ((res_ss_list[0][0] + 4) == res_ss_list[2][0]) & \
            ((res_ss_list[0][0] + 7) == res_ss_list[3][0]) ):
                return("FULL_HELIX")

        if (num_helix >= 2):
            if ((res_ss_list[0][0] + 1) == (res_ss_list[1][0])):
                if ((res_ss_list[0][1] == 'H') & (res_ss_list[1][1] == 'H')):
                    return("PARTIAL_HELIX")
            elif ((res_ss_list[1][0] + 1) == (res_ss_list[2][0])):
                if ((res_ss_list[1][1] == 'H') & (res_ss_list[2][1] == 'H')):
                    return("PARTIAL_HELIX")
            elif ((res_ss_list[2][0] + 1) == (res_ss_list[3][0])):
                if ((res_ss_list[2][1] == 'H') & (res_ss_list[3][1] == 'H')):
                    return("PARTIAL_HELIX")
        if (num_bridge >= 2):
            return("BRIDGE_BETA")
        if (num_strand == 2):
            return("PARTIAL_BETA")
        if (num_strand == 3):
            return("PARTIAL_BETA")
        if (num_strand == 4):
            return("FULL_BETA")
        if (num_turn >= 2):
            return("TURN_SS")

    if (type == t1111):
        if (num_helix == 4):
            if ((((res_ss_list[0][0]) + 4) == (res_ss_list[1][0])) & \
            (((res_ss_list[0][0]) + 7) == (res_ss_list[2][0])) & \
            (((res_ss_list[0][0]) + 11) == (res_ss_list[3][0]))):
                return("FULL_HELIX")
        if (num_strand >= 2):
            return("PARTIAL_BETA")
        if (num_bridge >= 2):
            return("BRIDGE_BETA")
        if (num_turn >= 2):
            return("TURN_SS")
    return("NO_SS")
    
