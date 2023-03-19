program hist4

!    INPUT params
!      in_tr
!      nrep, repmin, repmax
!      nk = array of length of data
!      teq = length of data to skip
!      tsim = length of data

implicit none
integer, parameter :: long = selected_real_kind(20,1400)
integer, parameter :: ilong = selected_int_kind(10)

integer, parameter :: repmin=1, repmax=20, nrep=20
integer, parameter :: nmon=21
integer, parameter :: in_tr=1
integer, parameter :: NUM_Z_BUCK=71
integer, parameter :: NUM_AA=21


real (long), parameter ::  p=1.458397d-05, RT=0.636d0, edelta=0.d0, temp0=320.d0, tmin=310.d0

integer, dimension(repmin:repmax), parameter :: nk=(/50000,50000,50000,50000,50000, 50000, 50000,50000,50000,50000,&
                                                    50000, 50000,50000,50000,50000, 50000, 50000,50000,50000,50000/)
integer, parameter :: tsim=50000, teq=0


real (long), dimension(nrep) :: restScale
real (long), dimension(repmin:repmax) :: fe, temp
real (long), dimension(nrep) :: templist

integer :: i, j, k, istr, repVal, nt, itermax, tprod, totStructs, depth, aa, xx
integer :: rep
integer :: currRep, currRepCnt
integer :: rep_val

real (long) :: f1, tempi, pf, hi, hr
character(len=16):: trString
character(len=16):: aaString
character(len=16):: repString
character(len=80), dimension(1:3) :: watContPath, statePath
real (long), dimension(200) :: dummy
integer :: vals, tens
real (long), dimension(:,:), allocatable :: en, pep_wat_cont_up, pep_wat_cont_lo
integer , dimension(:,:), allocatable :: pep_wat_cont_up_aa, pep_wat_cont_lo_aa
integer , dimension(:), allocatable :: pep_wat_cont_up_sum, pep_wat_cont_lo_sum
real (long), dimension(:,:), allocatable :: zsc, zsc2
integer, dimension(:,:), allocatable :: zsc_buck, zsc2_buck


integer, dimension(:), allocatable :: clust_up, clust_lo

real (long), dimension(repmin:repmax) :: pf_rep
real (long), dimension(0:5) :: pf_state
real (long), dimension(0:5, NUM_Z_BUCK) ::  pep_wat_cont_state
real (long), dimension(0:5, NUM_AA, NUM_Z_BUCK) ::  pep_wat_cont_state_aa, pep_wat_cont_state_aa_cnt
real (long), dimension(0:5, NUM_AA) ::  pep_wat_cont_state_aa_sum, pep_wat_cont_state_aa_sum_cnt
real (long), dimension(0:5) ::  pep_wat_cont_state_sum, pep_wat_cont_state_sum_cnt
integer :: stateUp, stateLo
integer :: bucket



type trajectory
    real, dimension(NUM_Z_BUCK) :: pep_wat_cont_up
    real, dimension(NUM_Z_BUCK) :: pep_wat_cont_lo
    integer, dimension(NUM_AA) :: pep_wat_cont_up_aa
    integer, dimension(NUM_AA) :: pep_wat_cont_lo_aa
    integer :: pep_wat_cont_sum_up
    integer :: pep_wat_cont_sum_lo
    integer :: clust_up
    integer :: clust_lo
end type trajectory
type (trajectory), dimension(:), allocatable  :: trj


statePath(1) = '/home/dad/GET_CLUST_STATE_HOT_1/'
statePath(2) = '/home/dad/GET_CLUST_STATE/'
statePath(3) = '/home/dad/GET_CLUST_STATE/'

watContPath(1) = '/home/dad/WATER_CONT_HOT_1/DATA/'
watContPath(2) = '/home/dad/WATER_CONT/DATA/'
watContPath(3) = '/home/dad/WATER_CONT/DATA/'

 totStructs=sum(nk)
 allocate( en(5,(tsim-teq)*(repmax-repmin+1)), &
  pep_wat_cont_up(NUM_Z_BUCK,(tsim-teq)*(repmax-repmin+1)), &
  pep_wat_cont_lo(NUM_Z_BUCK,(tsim-teq)*(repmax-repmin+1)), &
  pep_wat_cont_up_aa(NUM_AA, (tsim-teq)*(repmax-repmin+1)), &
  pep_wat_cont_lo_aa(NUM_AA, (tsim-teq)*(repmax-repmin+1)), &
  pep_wat_cont_up_sum((tsim-teq)*(repmax-repmin+1)), &
  pep_wat_cont_lo_sum((tsim-teq)*(repmax-repmin+1)), &
  zsc(nmon,(tsim-teq)*(repmax-repmin+1)),   &
  zsc2(nmon,(tsim-teq)*(repmax-repmin+1)),   &
  zsc_buck(nmon,(tsim-teq)*(repmax-repmin+1)),   &
  zsc2_buck(nmon,(tsim-teq)*(repmax-repmin+1)),   &
  clust_up((tsim-teq)*(repmax-repmin+1)), &
  clust_lo((tsim-teq)*(repmax-repmin+1)))

 tprod=tsim-teq
 allocate(trj(tprod))

!
!   Create trajectory string
 if (in_tr < 10) then
     trString = char(48+in_tr)
 else
     tens = in_tr/10
     vals = in_tr-(tens*10)
     trString = char(48+tens)//char(48+vals)
 endif

!
!  Read in replica temps.
 open(20,file='RepTemp.dat')
 do i=1,nrep
     read(20,*)templist(i)
     templist(i) = RT*templist(i)/temp0
 enddo
 close(20)
 temp=templist(repmin:repmax)

!
!  Read in scaling factors
 open(20,file='ScaleFactor.dat')
 do i=1,nrep
     read(20,*)restScale(i)
 enddo
 close(20)

!
!  Read in energy file
 open(34,file='energyi'//trim(trString)//'.dat',form='unformatted')

 do i=1,totStructs
     read(34)en(:,i)
 enddo
 close(34)
!
!  Read in FE data
 open(36,file='ffinal'//trim(trString)//'.dat')
 read(36,*)itermax,fe
 close(36)
!
!   Read the ZSC files
!   And save the bucket
 open(34,file='zsc1_'//trim(trString)//'.dat',form='unformatted')
 do i=1,totStructs
     read(34)zsc(:,i)
     do j=1,nmon
         bucket = zsc(j,i) + 35.0
         if (bucket < 1) then
             bucket = 1
         endif
         if (bucket > 71) then
             bucket = 71
         endif
         zsc_buck(j,i) = bucket
     enddo
 enddo
 close(34)

 open(34,file='zsc2_'//trim(trString)//'.dat',form='unformatted')
 do i=1,totStructs
     read(34)zsc2(:,i)
     do j=1,nmon
         bucket = zsc2(j,i) + 35.0
         if (bucket < 1) then
             bucket = 1
         endif
         if (bucket > 71) then
             bucket = 71
         endif
         zsc2_buck(j,i) = bucket
     enddo

 enddo
 close(34)




!
!   Read in data files
 istr=0
 do rep=repmin,repmax
     if (rep < 10) then
         repString = char(48+rep)
     else
         tens = rep/10
         vals = rep-(tens*10)
         repString = char(48+tens)//char(48+vals)
     endif
!
! Water contacts per peptide


     open(21,file=trim(watContPath(in_tr))//'pep_cont'//trim(trString)//'_$RANGE_'//trim(repString)//'_up.dat')
     open(22,file=trim(watContPath(in_tr))//'pep_cont'//trim(trString)//'_$RANGE_'//trim(repString)//'_lo.dat')
     open(23,file=trim(statePath(in_tr))//'cont_state'//trim(trString)//'_$RANGE_'//trim(repString)//'.dat')
     open(24,file=trim(WatContPath(in_tr))//'pep_cont_sum'//trim(trString)//'_$RANGE_'//trim(repString)//'.dat')


!    Read in any header data

!   Data to Skip
     do k=1,teq
       read(21,*)dummy(1:71)
       read(22,*)dummy(1:71)
       read(23,*)dummy(1:2)
     enddo
!   Read in good data
     do k=teq+1,tsim
       print*,"pep_cont ", rep, k
       flush(6)
       read(21,*)trj(k)%pep_wat_cont_up(:)
       print*, 'Added ', trj(k)%pep_wat_cont_up(:)
       print*,"DONE ", rep, k
       flush(6)
       read(22,*)trj(k)%pep_wat_cont_lo(:)
       read(23,*)trj(k)%clust_up, trj(k)%clust_lo
       read(24,*)trj(k)%pep_wat_cont_sum_up, trj(k)%pep_wat_cont_sum_lo
     enddo

    close(21)
    close(22)
    close(23)
    close(24)
!
!  Open AA cont files
    print*,'READ AA UP files', rep
    flush(6)
    do aa=1,NUM_AA
        if (aa < 10) then
            aaString = char(48+aa)
        else
            tens = aa/10
            vals = aa-(tens*10)
            aaString = char(48+tens)//char(48+vals)
        endif

        open(21,file=trim(watContPath(in_tr))//'pep_contAA'//trim(trString)//&
          '_$RANGE_'//trim(repString)//'_'//trim(aaString)//'.dat')

!   Data to Skip
        do k=1,teq
            read(21,*)dummy(1:2)
        enddo
!   Read in good data
        do k=teq+1,tsim
           read(21,*)trj(k)%pep_wat_cont_up_aa(aa), trj(k)%pep_wat_cont_lo_aa(aa)
        enddo

        close(21)
    enddo


!   print*,'READ AA LO files', rep
!   flush(6)
!   do aa=1,NUM_AA
!        if (aa < 10) then
!            aaString = char(48+aa)
!        else
!            tens = aa/10
!            vals = aa-(tens*10)
!            aaString = char(48+tens)//char(48+vals)
!        endif
!
!        open(21,file=trim(watContPath(in_tr))//'pep_contAA'//trim(trString)//&
!          '_$RANGE_'//trim(repString)//'_'//trim(aaString)//'_up.dat')
!
!!   Data to Skip
!        do k=1,teq
!            read(21,*)dummy(1:71)
!        enddo
!!   Read in good data
!        do k=teq+1,tsim
!           read(21,*)trj(k)%pep_wat_cont_up_aa(:,aa)
!        enddo
!
!        close(21)
!    enddo



    print*,'DONE READ AA files',rep
    flush(6)
    
!
! Get data for one replica into variable
    tprod=tsim-teq
    do k=1, tprod
      istr = istr+1
      pep_wat_cont_up(1:NUM_Z_BUCK,istr) =  trj(k)%pep_wat_cont_up
      pep_wat_cont_lo(1:NUM_Z_BUCK,istr) =  trj(k)%pep_wat_cont_lo
      pep_wat_cont_up_aa(1:NUM_AA,istr) =  trj(k)%pep_wat_cont_up_aa
      pep_wat_cont_lo_aa(1:NUM_AA,istr) =  trj(k)%pep_wat_cont_lo_aa
      pep_wat_cont_up_sum(istr) = trj(k)%pep_wat_cont_sum_up
      pep_wat_cont_lo_sum(istr) = trj(k)%pep_wat_cont_sum_lo
      clust_up(istr) =  trj(k)%clust_up
      clust_lo(istr) =  trj(k)%clust_lo
    enddo
 enddo   ! enddo for each replica

    print*,'DONE reading. zero vars'
    flush(6)

    do i=0,5
        pf_state(i) = 0.d0
        pep_wat_cont_state_sum(i) =  0.d0
        pep_wat_cont_state_sum_cnt(i) =  0.d0
        do j=1,NUM_Z_BUCK
            pep_wat_cont_state(i,j) = 0.d0
            do aa=1,NUM_AA
                pep_wat_cont_state_aa(i,aa,j) = 0.d0
                pep_wat_cont_state_aa_cnt(i,aa,j) = 0.d0
                pep_wat_cont_state_aa_sum(i,aa) = 0.d0
                pep_wat_cont_state_aa_sum_cnt(i,aa) = 0.d0
            enddo
        enddo
    enddo
    do i=1,20
        pf_rep(i) = 0.d0
    enddo

    pf=0.d0
    nt=nint(tmin)
    tempi = real(nt,long)*RT/temp0

!....computing partition function...........
    print*,'Compute partition func'
    flush(6)
    currRep = repmin
    currRepCnt = 0
    do i = 1,totStructs
        f1 = 0.d0
        hi = en(1,i) + en(2,i) + en(3,i) + p*en(5,i) + edelta
        do repVal=repmin,repmax
            hr = en(1,i) + &
              sqrt(restScale(repVal))*en(2,i) + &
              (restScale(repVal))*en(3,i) + &
              p*en(5,i) + edelta
            f1 = f1 + real(nk(repVal),long)*exp( fe(repVal) - hr/temp(repVal)+hi/tempi)
!            print*,i, repVal, f1, exp( fe(repVal) - hr/temp(repVal)+hi/tempi),  fe(repVal), hr/temp(repVal), hi/tempi
        enddo
        if(f1 > huge(real(100,long)))then
            print*,'f1=inf, repVal, i and m=',f1,i,repVal !; stop
            flush(6)
        endif
        


        stateUp = clust_up(i)
        stateLo = clust_lo(i)
        print*,'Done compute ', i,  1.d0/f1, ' stateup ', stateUp, ' statelo ', stateLo
        flush(6)
        pf_state(stateUp) = pf_state(stateUp) + 1.d0/f1
        pf_state(stateLo) = pf_state(stateLo) + 1.d0/f1
        pep_wat_cont_state(stateUp,:) = pep_wat_cont_state(stateUp,:) + pep_wat_cont_up(:,i)/f1
        pep_wat_cont_state(stateLo,:) = pep_wat_cont_state(stateLo,:) + pep_wat_cont_lo(:,i)/f1

        pep_wat_cont_state_sum(stateUp) = pep_wat_cont_state_sum(stateUp) + pep_wat_cont_up_sum(i)/f1
        pep_wat_cont_state_sum(stateLo) = pep_wat_cont_state_sum(stateLo) + pep_wat_cont_up_sum(i)/f1
        do k=1,NUM_AA
            bucket = zsc_buck(k,i)
!            print*, 'Bucket', k, i, bucket, stateUp
            flush(6)
            pep_wat_cont_state_aa(stateUp,k,bucket) = pep_wat_cont_state_aa(stateUp,k,bucket) +&
              (real(pep_wat_cont_up_aa(k,i),(long))/f1)
            pep_wat_cont_state_aa_cnt(stateUp,k,bucket) = pep_wat_cont_state_aa_cnt(stateUp,k,bucket) +&
              (1.d0/f1)
            pep_wat_cont_state_aa_sum(stateUp,k) = pep_wat_cont_state_aa_sum(stateUp,k) +&
              (real(pep_wat_cont_up_aa(k,i),(long))/f1)
            pep_wat_cont_state_aa_sum_cnt(stateUp,k) = pep_wat_cont_state_aa_sum_cnt(stateUp,k) +&
              (1.d0/f1)
        enddo
        do k=1,NUM_AA
            bucket = zsc2_buck(k,i)
!            print*, 'Bucket2', k, i, bucket, stateLo
            flush(6)
            pep_wat_cont_state_aa(stateLo,k,bucket) = pep_wat_cont_state_aa(stateLo,k,bucket) +&
              (real(pep_wat_cont_lo_aa(k,i),(long))/f1)
            pep_wat_cont_state_aa_cnt(stateLo,k,bucket) = pep_wat_cont_state_aa_cnt(stateLo,k,bucket) +&
              (1.d0/f1)
            pep_wat_cont_state_aa_sum(stateLo,k) = pep_wat_cont_state_aa_sum(stateLo,k) +&
              (real(pep_wat_cont_lo_aa(k,i),(long))/f1)
            pep_wat_cont_state_aa_sum_cnt(stateLo,k) = pep_wat_cont_state_aa_sum_cnt(stateLo,k) +&
              (1.d0/f1)
        enddo
        rep_val =  (i / tprod) + 1
        pf_rep(rep_val) = pf_rep(rep_val) +  1.d0/f1
        pf = pf + 1.d0/f1

    enddo

    print*, 'pf state', pf_state, pf
    flush(6)
    do i = 0, 5
        pep_wat_cont_state(i,:) = pep_wat_cont_state(i,:) / pf_state(i)
        pep_wat_cont_state_sum(i) = pep_wat_cont_state_sum(i) / pf_state(i)
        do aa=1,NUM_AA
            do j=1,NUM_Z_BUCK
                if (pep_wat_cont_state_aa_cnt(i,aa,j) == 0) then
                    pep_wat_cont_state_aa(i,aa,j) = -1
                else
!                    pep_wat_cont_state_aa(i,aa,j) = pep_wat_cont_state_aa(i,aa,j) /pep_wat_cont_state_aa_cnt(i,aa,j)
                    pep_wat_cont_state_aa(i,aa,j) = pep_wat_cont_state_aa(i,aa,j) /pf_state(i)
                endif
            enddo
            pep_wat_cont_state_aa_sum(i,aa) = pep_wat_cont_state_aa_sum(i,aa) /pf_state(i)
        enddo
    enddo
    pf_rep = pf_rep / pf

    do i = 0, 5
        open(100,file='th_pep_wat_cont'//trim(trString)//'_$RANGE_state'//char(48+i)//'.dat')
        do j = 1,71
            depth = j - 36
            write(100,*)depth, real(pep_wat_cont_state(i,j), long)
        enddo
        close(100)
    enddo

    open(100,file='th_pep_wat_cont_sum'//trim(trString)//'_$RANGE.dat')
    do i = 0, 5
         write(100,*)i, pep_wat_cont_state_sum(i)
    enddo
    close(100)


    open(100,file='th_pf_reps'//trim(trString)//'_$RANGE.dat')
    do i = 1, nrep
        write(100,*)real(pf_rep(i),long)
    enddo
    close(100)

    do i = 0, 5
        open(100,file='th_pep_wat_cont_aa'//trim(trString)//'_$RANGE_state'//char(48+i)//'.dat')
        do aa = 1, 21
            write(100,*)real(pep_wat_cont_state_aa(i,aa,:), long)
        enddo
        close(100)
    enddo

    open(100,file='th_pep_wat_cont_aa_sum'//trim(trString)//'_$RANGE.dat')

    do aa = 1, 21
        write(100,*)real(pep_wat_cont_state_aa_sum(0,aa), long),&
          real(pep_wat_cont_state_aa_sum(1,aa), long),&
          real(pep_wat_cont_state_aa_sum(2,aa), long),&
          real(pep_wat_cont_state_aa_sum(3,aa), long),&
          real(pep_wat_cont_state_aa_sum(4,aa), long),&
          real(pep_wat_cont_state_aa_sum(5,aa), long)
    enddo
    close(100)

end program hist4

