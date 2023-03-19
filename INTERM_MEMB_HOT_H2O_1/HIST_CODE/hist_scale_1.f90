
program hist7
!    INPUT params
!      trmin, trmax
!      nrep, repmin, repmax
!      tsim = number of structs
!      teq = number needed to equilibrate
!      nk = tsim - teq
!      numStates
!
!    OUTPUT:
!      th_cm.dat - average contact map for each trajectory
!      th_eppw.dat - average energy (from p-->p and p-->w contacts, no w-->w)
!      th_tt.dat - average turns
!      th_hel.dat - average hel
!      th_str.dat - Average SS
!      th_ran.dat - Ramach plot
!      th_psi.dat - psi hist
!      th_phi.dat - phi hist

implicit none
integer, parameter :: long = selected_real_kind(20,1400)
!integer, parameter :: long = selected_real_kind(15,307)
integer, parameter :: ilong = selected_int_kind(10)
integer, parameter :: repmin=1, repmax=20, nrep=20, dnt=2, nmol=1
integer, parameter :: in_tr=1
integer, parameter :: trmin=1
integer, parameter :: trmax=1
!integer, dimension(trmin:trmax), parameter :: tsim=(/9999,9999,9999,9999,9999,9999,9999,9999,9999,9999/),&
!                                              teq=(/0,0,0,0,0,0,0,0,0,0/)
integer, parameter :: tsim=50000, teq=0
REAL, PARAMETER :: M_PI = 3.1415927

integer, parameter :: nmon=21, ncont=210, ncont_leaf=105, ncont_pep=210
real (long), parameter :: eps=1.d-7, e1=-4000.d0, e2=10000.d0, edist = 2.d0, &
                          p=1.458397d-05, RT=0.636d0, edelta=0.d0, &
                          temp0=320.d0, emid=0.d0, tmin=310.d0,tmax=310.d0, z21_clust=17.556, &
                          zSurface=17.556, zUnbound=24.056

real (long), dimension(nrep) :: restScale

integer, dimension(repmin:repmax), parameter :: nk=(/50000,50000,50000,50000,50000, 50000, 50000,50000,50000,50000,&
                                                    50000,50000,50000,50000,50000, 50000, 50000,50000,50000,50000/)
real (long), dimension(:,:), allocatable :: en



real (long), dimension(nrep) :: templist
real (long) :: f1, tempi, pf,hi,hr, minPf, maxPf, maxF1
real (long), dimension(repmin:repmax) :: fe, temp
real (long), dimension(repmin:repmax) :: totUsed
real (long) :: usedTotal




integer :: currRep, currRepCnt
integer ::  nt1=nint(tmin), nt2=nint(tmax)

integer :: i,j,k,l,m,n,aa, itermax, nt, istr, helCnt, betaCnt, helCnt2, betaCnt2

character(len=16):: trString
character(len=16):: repString

integer :: vals, tens

allocate( en(5,(tsim-teq)*(repmax-repmin+1)))
!
!   Create trajectory string
 if (in_tr < 10) then
     trString = char(48+in_tr)
 else
     tens = in_tr/10
     vals = in_tr-(tens*10)
     trString = char(48+tens)//char(48+vals)
 endif




 open(20,file='RepTemp.dat')
  do i=1,nrep
   read(20,*)templist(i)
   templist(i) = RT*templist(i)/temp0
  enddo
 close(20)
 temp=templist(repmin:repmax)




 open(20,file='ScaleFactor.dat')
  do i=1,nrep
   read(20,*)restScale(i)
  enddo
 close(20)





 istr=sum(nk)

 open(34,file='energyi'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)en(:,i)
!    print*, 'READ EN ', i, en(1,i), en(2,i)
  enddo
 close(34)



 do i=repmin,repmax
   print*,i,real(temp(i)),nk(i)  
 enddo
 


 open(36,file='ffinal'//trim(trString)//'.dat')
 read(36,*)itermax,fe
 print*, 'FE ', fe
 close(36)

!  do nt=nt1,nt2,dnt (do temperatures)
   nt=nt1

   pf=0.d0

   tempi = real(nt,long)*RT/temp0


!....computing partition function...........
   currRep = repmin
   currRepCnt = 0

   if (currRep < 10) then
       repString = char(48+currRep)
   else
       tens = currRep/10
       vals = currRep-(tens*10)
       repString = char(48+tens)//char(48+vals)
    endif


  open(136,file='hist_scale_'//trim(trString)//'_'//trim(repString)//'.dat')

  minPf = 1.0
  maxPf = 0.0
  maxf1 = 0.0
  do i = 1,istr

     f1 = 0.d0
     hi = en(1,i) + en(2,i) + en(3,i) + p*en(5,i) + edelta
!     print*,'ENE ', hi, en(1,i), en(2,i), en(3,i), p*en(5,i)
     do m=repmin,repmax
      hr = en(1,i) + &
          sqrt(restScale(m))*en(2,i) + &
          (restScale(m))*en(3,i) + &
          p*en(5,i) + edelta
      f1 = f1 + real(nk(m),long)*exp( fe(m) - hr/temp(m)+hi/tempi)
!      print*, 'REP ', i, ' F1 ', f1, fe(m), hr, temp(m), hi, tempi
     enddo
     if(f1 > huge(real(100,long)))then
       print*,'f1=inf, f1, i and m=',f1,i,m !; stop
     endif
     pf = pf + real(1.d0/f1,long)

!     print*,'F1 ', real(f1,long), real(pf, long)


    write(136,*)real((1.d0/f1),long)
    if (f1 > maxF1) then
        maxF1 = f1
    endif
    if (real((1.d0/f1),long) > real(maxPf,long)) then
        maxPf = real((1.d0/f1),long)
    endif
    totUsed(currRep) = totUsed(currRep) + (1.0/f1)
    currRepCnt = currRepCnt + 1
     if (currRepCnt >= nk(currRep)) then
         currRepCnt = 0
         close(136)
         if (currRep < repmax) then
             currRep = currRep + 1
             if (currRep < 10) then
                 repString = char(48+currRep)
             else
                 tens = currRep/10
                 vals = currRep-(tens*10)
                 repString = char(48+tens)//char(48+vals)
              endif
             open(136,file='hist_scale_'//trim(trString)//'_'//trim(repString)//'.dat')
         endif
     endif
  enddo

  usedTotal = sum(totUsed)

  open(49,file='th_rep_used_'//trim(trString)//'.dat')

  do currRep = repmin,repmax
      write(49,*)currRep,real(totUsed(currRep)/usedTotal, long)
  enddo
  close(49)


  minPf = real((1.d0/maxF1), long)
  print*,'MIN PF', real((minPf),long), 'MAX PF ', real((maxPf),long), 'MAX F1 ', maxF1



end program hist7


