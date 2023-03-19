!
!   Need input nrep
module parameters
 implicit none
 integer, parameter :: long = selected_real_kind(15,307)
 integer :: mseed
end  module parameters

program schedule
 use parameters
 implicit none
! integer, dimension(1:nrep) :: partner
 integer, dimension(1:20) :: partner
 integer :: run,repstep,nrep,l,i
 character(len=3) :: ch
 real(long) :: ran2

 read*,repstep,run,nrep
 mseed=-run*181800 - repstep
 l = int(2*ran2(mseed)) + 1
 do i=1,nrep,1
  partner(i)=i
 end do
 do i=l,nrep-1,2
  partner(i)=i+1
  partner(i+1)=i
 end do
 write(ch,'(i0)') nrep+1
! write(*,'('//trim(ch)//'i4)') l, partner
 write(*,'('//trim(ch)//'i4)') l, (partner(i), i=1,nrep)
end program schedule

!--------------------------------------------------------------------

function  ran2(idum)

 use parameters, only : long
 implicit none
 integer, intent(inout)  :: idum
 integer, parameter :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1, &
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
    NTAB=32,NDIV=1+IMM1/NTAB
 real (long), parameter :: AM=1.D0/IM1,ESP1=1.2D-7,RNMX=1.D0-ESP1
 integer :: idum2,j,k,iy
 integer, dimension(NTAB) :: iv
 real (long) :: ran2

 save iv,iy,idum2
 data idum2/123456789/, iv/NTAB*0/, iy/0/

 if(idum <= 0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum < 0) idum=idum+IM1
        if (j <= NTAB) iv(j)=idum
     enddo
     iy=iv(1)
 endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum < 0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if(idum2 < 0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2  
    iv(j)=idum
    if(iy < 1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)

end function  ran2 


