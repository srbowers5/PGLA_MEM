program hist2
!    INPUT params
!      trmin, trmax
!      nrep, repmin, repmax
!      itmax = max iterations
!      be1, be2 not used
!      tsim = number of structs
!      teq = number needed to equilibrate
!      eps = max error
!      e1, e2, edist not used
!
!   OUTPUT:
!      ffinal.dat (free energy at each step)
implicit none 
integer, parameter :: long = selected_real_kind(20,1400)
!integer, parameter :: long = selected_real_kind(15,700)
!integer, parameter :: long = selected_real_kind(15,307)
integer, parameter :: itmax=10000
!integer, parameter :: be1=0,be2=7000
integer, parameter :: nrep=20, tp=nrep, repmin=1,repmax=20
integer, parameter :: in_tr=1
integer, parameter :: trmin=1
integer, parameter :: trmax=1
integer, dimension(trmin:trmax), parameter :: tsim=(/50000/), teq=(/0/)

real (long), parameter :: eps=1.d-7, e1=-4000.d0, e2=10000.d0, edist = 2.d0,&
                          RT=0.636d0, temp0=320.d0, edelta=0.d0, p=1.458397d-05
real (long), dimension(:,:), allocatable :: en
real (long), dimension(tp) :: fe=0.d0,fenew,dfe=0.d0,temp,restScale
!integer, dimension(tp) :: nk
!integer, dimension(repmin:repmax), parameter :: nk=(/50000,50000,50000,50000,50000, 50000,50000,50000,50000,50000,50000,50000/)
integer, dimension(repmin:repmax), parameter :: nk=(/50000,50000,50000,50000,50000, &
 50000,50000,50000,50000,50000, 50000,50000,50000,50000,50000, 50000,50000,50000,50000,50000/)

real (long) :: f1, t0, maxdf = 10000.d0, hi, hr


integer :: i,j,k,k0,l,m, icount,eh,ehb,ehtots,istr,nstr,hibin,hrbin

allocate(en(5,sum(tsim-teq)*(repmax-repmin+1)))

 open(20,file='RepTemp.dat')
  do i=1,nrep
   read(20,*)temp(i)
   temp(i) = RT*temp(i)/temp0
  enddo
 close(20)


 open(20,file='ScaleFactor.dat')
  do i=1,nrep
   read(20,*)restScale(i)
  enddo
 close(20)

 open(20,file='finit.dat')
  do i=1,nrep
   read(20,*)fe(i)
  enddo
 close(20)


! nk is total number of structs at each temp
 nstr = sum(nk)

 ! read*,in_tr
 do i=1,tp
   print*,temp(i),nk(i)  
 enddo
 
  if (in_tr < 10) then
      open(34,file='energyi'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='energyi1'//char(38+in_tr)//'.dat',form='unformatted')
   endif

  do i=1,nstr
    read(34)en(:,i)
!    print*, "ENERGY ", i, en(:,i)
  enddo
  close(34)


! Creating energy histogram ............

!  do i=1,istr
!    ehtots = hcm(ncont+4,i)
!    he(ehtots) = he(ehtots) + 1 
!  enddo

!  print*,sum(he),sum(nk)

!  Starting the cycles of iteration .....
  if (in_tr < 10) then
      open(36,file='f'//char(48+in_tr)//'.dat')
      open(35,file='df'//char(48+in_tr)//'.dat')
  else
      open(36,file='f1'//char(38+in_tr)//'.dat')
      open(35,file='df1'//char(38+in_tr)//'.dat')
  endif

 k0=0
 do while(eps < maxdf)
   k0 = k0 + 1

   fenew = 0.d0  
   
   do l=1,tp
   
    do j=1,nstr
     
     f1 = 0.d0
! en(1 = protein-protein)
! en(2 = protein water)
! en(3 = water water)
! en(5 = energy Volume)
    hi = en(1,j) + &
          sqrt(restScale(l))*en(2,j) + &
          (restScale(l))*en(3,j) + &
          p*en(5,j) + edelta
!     if (k0 < 3) then
!         print*, "EN1 ", k0, j, en(1,j), en(2,j), en(3,j), en(5,j), hi
!    endif
     do m=1,tp
      hr = en(1,j) + &
          sqrt(restScale(m))*en(2,j) + &
          (restScale(m))*en(3,j) + &
          p*en(5,j) + edelta
!         if (k0 < 3) then
!             print*, "EN2 ", k0, j, tp,  hr
!         endif
      f1 = f1 + real(nk(m),long)*exp( fe(m)-hr/temp(m)+hi/temp(l))
     enddo
     
     fenew(l) = fenew(l) + 1.d0 / f1
    enddo
    
   enddo

   dfe =- log(fenew) - fe
   print*,"DFE ",dfe
   print*,"FE NEW",fenew
   print*,"FE ",fe
   print*,"LOG FENEW" ,log(fenew)
   maxdf = maxval( abs(dfe) )

   print*,"MAX DF ", k0,maxdf
   fe = -log(fenew)

   write(35,*)k0,dfe
   write(36,*)k0,fe
   
  if(k0 == itmax) exit
 enddo

close(36)
close(35)

  if (in_tr < 10) then
      open(38,file='ffinal'//char(48+in_tr)//'.dat')
  else
      open(38,file='ffinal1'//char(38+in_tr)//'.dat')
  endif
   write(38,*)k0,fe
 close(38)

end program hist2

