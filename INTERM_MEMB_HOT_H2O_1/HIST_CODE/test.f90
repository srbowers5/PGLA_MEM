program test1

implicit none
integer, parameter :: long = selected_real_kind(20,1400)

integer :: i,j
real (long), dimension(0:21,0:119) :: betaLen, betaRg

do j=0,21
 do i=0,119
   betaLen(j,i) = 0.0
   print*,j,i
 enddo
enddo


end program test1

