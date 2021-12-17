# addgroup needs atomid
#

set upperlayer [ addgroup {13 131 249 360 471 582 700 818 936 1054 1172 1283 1394 1512 1630 1748 1859 1970 2081 2199 2317 2428 2546 2664 2782 2900 3018 3136 3247 3365 3476 3594 3705 3816 3927 4045 4163 4274 4392 4503 4621 4739 4857 4968 5086 5204 5322 5433 5551}]

set lowerlayer [ addgroup {5669 5787 5898 6016 6134 6245 6363 6474 6585 6703 6821 6939 7057 7175 7293 7411 7529 7640 7758 7869 7987 8105 8223 8341 8459 8577 8688 8799 8910 9028 9146 9264 9375 9493 9611 9729 9840 9958 10069 10180 10291 10402 10520 10631 10749 10860 10971 11082 11200}]


# spring constant
set k 6.5

# set distance between P-P
set d 35.11127


# required function
proc calcforces {} {
 # load in atom coordinates (add group computes COM)
 global upperlayer lowerlayer k d
 loadcoords coord

 # constraint between cm1 and cm2 of phosphate atoms in bilayer
 set cm1 [split $coord($upperlayer) { }]
 set cm2 [split $coord($lowerlayer) { }]
 set cm1_z [lindex $cm1 2]
 set cm2_z [lindex $cm2 2]
 set r_cm1_d [expr $cm1_z-($d/2.0)]
 set r_cm2_d [expr $cm2_z+($d/2.0)]
 set f_cm1_z [expr -$k*$r_cm1_d]
 set f_cm2_z [expr -$k*$r_cm2_d]
 addenergy [expr 0.5*$k*$r_cm1_d**2.0]
 addenergy [expr 0.5*$k*$r_cm2_d**2.0]
 addforce $upperlayer [list 0.0 0.0 $f_cm1_z]
 addforce $lowerlayer [list 0.0 0.0 $f_cm2_z]


# print cm1=$cm1 cm1_z=$cm1_z r_cm1_d=$r_cm1_d f_cm1_z=$f_cm1_z
# print cm2=$cm2 cm2_z=$cm2_z r_cm2_d=$r_cm2_d f_cm2_z=$f_cm2_z
# print cm1_z=$cm1_z

}
