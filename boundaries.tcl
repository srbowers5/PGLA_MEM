# wrapmode input
# wrapmode cell
# wrapmode nearest
# wrapmode patch ;# the default

#   Need to set the start and end of each peptide in the code.
#   zw = distance from top/bottom of cell where we start to add forces.
#   k = spring constant
proc calcforces {step unique start_pgl1 end_pgl1 start_pgl2 end_pgl2 start_sod1 end_sod1 start_sod2 end_sod2 start_cla1 end_cla1 start_cla2 end_cla2 zdist k} {

 set cell [getcell]
 set cellz [lindex [lindex $cell 3] 2]

# print  "Cell Z $cellz"

 set cellZ [expr $cellz / 2.0]
 set zw [expr $cellZ-$zdist]

 while {[nextatom]} {
  set atomnum [getid]

  # PGL 1, (PGL1 placed above membrane in +Z)
  if {$atomnum >= $start_pgl1 && $atomnum <= $end_pgl1} {
   set coord [getcoord]
   set zk [lindex $coord 2]
   if {$zk > $zw} {
    set rz [expr $zk-$zw]
    set fz [expr -$k*$rz]
    addenergy [expr 0.5*$k*$rz**2.0]
    addforce [list 0.0 0.0 $fz]
   }
  # SOD top, (SOD placed above membrane in +Z)
  } elseif {$atomnum >= $start_sod1 && $atomnum <= $end_sod1} {
   set coord [getcoord]
   set zk [lindex $coord 2]
   if {$zk > $zw} {
    set rz [expr $zk-$zw]
    set fz [expr -$k*$rz]
    addenergy [expr 0.5*$k*$rz**2.0]
    addforce [list 0.0 0.0 $fz]
   }
  # CLA top, (CLA placed above membrane in +Z)
  } elseif {$atomnum >= $start_cla1 && $atomnum <= $end_cla1} {
   set coord [getcoord]
   set zk [lindex $coord 2]
   if {$zk > $zw} {
    set rz [expr $zk-$zw]
    set fz [expr -$k*$rz]
    addenergy [expr 0.5*$k*$rz**2.0]
    addforce [list 0.0 0.0 $fz]
   }
  # PGL 2, (PGL placed below membrane in -Z)
  } elseif {$atomnum >= $start_pgl2 && $atomnum <= $end_pgl2} {
   set coord [getcoord]
   set zk [lindex $coord 2]
   if {$zk < [expr $zw*-1.0]} {
    set rz [expr $zk+$zw]
    set fz [expr -$k*$rz]
#    print "calcforces Add neg energy  $atomnum, $zdist, $cellZ, $zw, $zk, $k, $rz, $fz "
    addenergy [expr 0.5*$k*$rz**2.0]
    addforce [list 0.0 0.0 $fz]
   }
  # SOD 2, (SOD placed below membrane in -Z)
  } elseif {$atomnum >= $start_sod2 && $atomnum <= $end_sod2} {
   set coord [getcoord]
   set zk [lindex $coord 2]
   if {$zk < [expr $zw*-1.0]} {
    set rz [expr $zk+$zw]
    set fz [expr -$k*$rz]
#    print "calcforces Add neg energy  $atomnum, $zdist, $cellZ, $zw, $zk, $k, $rz, $fz "
    addenergy [expr 0.5*$k*$rz**2.0]
    addforce [list 0.0 0.0 $fz]
   }
  # CLA 2, (CLA placed below membrane in -Z)
  } elseif {$atomnum >= $start_cla2 && $atomnum <= $end_cla2} {
   set coord [getcoord]
   set zk [lindex $coord 2]
   if {$zk < [expr $zw*-1.0]} {
    set rz [expr $zk+$zw]
    set fz [expr -$k*$rz]
#    print "calcforces Add neg energy  $atomnum, $zdist, $cellZ, $zw, $zk, $k, $rz, $fz "
    addenergy [expr 0.5*$k*$rz**2.0]
    addforce [list 0.0 0.0 $fz]
   }
  # drop all other atoms from consideration
  } else {
   dropatom
   continue
 }
 }
}


