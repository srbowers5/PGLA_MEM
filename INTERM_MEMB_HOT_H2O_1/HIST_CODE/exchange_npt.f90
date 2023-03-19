!
!  Need to create parameters temp0, na, nrep
!    Now read temp0 and nrep
module parameters
 implicit none

 integer, parameter :: long = selected_real_kind(15,307)
 integer, parameter :: na= 25184
 real(long), parameter :: RT=0.616d0, PN=1.458397d-05
 real(long), dimension(na,3) :: vel1, vel2
 real(long), dimension(20) :: temp, energy, energyswap
 real(long), dimension(2) :: energy_i
 real(long), dimension(3) :: box_i
 integer :: mseed
 character(len=40) :: file1,file2
end  module parameters

program exchange
 use parameters
 implicit none
 integer :: i,j,m,n,run,repstep,nrep,l,l0,na0
 real(long) :: ran2,scale,r,expd,delta,temp0
 real(long) :: penergy 
 logical :: Accept

 read*,repstep,run,nrep,temp0,l
 mseed=-run*181800 - repstep
100 format(a6,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,a26)

 open(20,file='RepTemp.dat')
  do i=1,nrep
   read(20,*) temp(i)
   temp(i) = RT*temp(i)/temp0
  enddo
 close(20)

 open(20,file='RepEnergy.dat')
 open(21,file='RepEnergySwap.dat')
  do i=1,nrep
   read(20,*) energy_i(:)
   energy(i) = energy_i(2) - energy_i(1)
   read(21,*) energy_i(:)
   energyswap(i) = energy_i(2) - energy_i(1)
  enddo
 close(20)
 close(21)

 open(20,file='RepBox.dat')
 open(21,file='RepBoxSwap.dat')
  do i=1,nrep
   read(20,*) box_i(:)
   penergy = PN*product(box_i(:))
   energy(i) = energy(i) + PN*product(box_i(:))
   read(21,*) box_i(:)
   energyswap(i) = energyswap(i) + PN*product(box_i(:))
  enddo
 close(20)
 close(21)

 l0 = int(2*ran2(mseed)) + 1
 if(l.ne.l0) then
  print*, 'Error: l /= l0.'
  stop 1
 end if 
 do m=l,nrep-1,2
  Accept=.false.
  n = m+1
  delta = (energyswap(n)-energy(m))/temp(m) + (energyswap(m)-energy(n))/temp(n)
  if(delta <= 0.d0) then 
    Accept=.true.
  else
    r = ran2(mseed)
    expd = exp(-delta)
    if(r < expd)  Accept=.true.
  endif
  if(Accept)then
    print*,m,n,1
    if(m < 10) file1 = 'pgl_mem'//char(48+m)
    if(m >= 10) file1 = 'pgl_mem'//char(48+m/10)//char(48+m-m/10*10)
    if(n < 10) file2 = 'pgl_mem'//char(48+n)
    if(n >= 10) file2 = 'pgl_mem'//char(48+n/10)//char(48+n-n/10*10)

    !open(21,file=trim(file1)//'.vel',form='binary')
    !open(22,file=trim(file2)//'.vel',form='binary')
    open(21,file=trim(file1)//'.vel',access='stream',form='unformatted')
    open(22,file=trim(file2)//'.vel',access='stream',form='unformatted')

    read(21) na0
    if(na.ne.na0) then
     print*, 'Error: na /= na0.'
     stop 2
    end if
    read(22) na0
    if(na.ne.na0) then
     print*, 'Error: na /= na0.'
     stop 3
    end if
    read(21) ((vel1(i,j),j=1,3),i=1,na)
    read(22) ((vel2(i,j),j=1,3),i=1,na)
    close(21)
    close(22)

    vel1=vel1*20.45482706 ! convert from NAMD internal units
    vel2=vel2*20.45482706 ! convert from NAMD internal units

    scale=sqrt(temp(n)/temp(m))
    !open(21,file=trim(file1)//'s.vel',form='binary')
    !open(22,file=trim(file2)//'s.vel',form='binary')
    open(21,file=trim(file1)//'s.vel',access='stream',form='unformatted')
    open(22,file=trim(file2)//'s.vel',access='stream',form='unformatted')

    write(21) na
    write(22) na
    write(21) ((vel2(i,j)/scale/20.45482706,j=1,3),i=1,na)
    write(22) ((vel1(i,j)*scale/20.45482706,j=1,3),i=1,na)
    close(21)
    close(22)

  else
    print*,m,n,0
  endif
 enddo

end program exchange

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


