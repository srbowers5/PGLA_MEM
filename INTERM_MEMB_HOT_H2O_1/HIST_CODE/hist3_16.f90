
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
integer, parameter :: numStates=5
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
                                                    50000, 50000,50000,50000,50000, 50000, 50000,50000,50000,50000/)
integer, dimension(:,:), allocatable :: hcm_pc_up, hcm_pg_up, hcm_pc_low, hcm_pg_low, heppw
integer, dimension(:,:), allocatable :: scdPgRegCnt, scdPcRegCnt
real (long), dimension(:,:), allocatable :: scdPgReg, scdPcReg
real (long), dimension(:,:), allocatable :: regLipDen
integer, dimension(:,:), allocatable :: hcm_pep1, hcm_pep2
integer, dimension(:,:), allocatable :: hst, hhel, hran, hstr
real (long), dimension(:,:), allocatable :: zsc, zsc2
real (long), dimension(:,:), allocatable :: ro_12, ro2_12, ro_15, ro2_15, ro_19, ro2_19

integer, dimension(:,:), allocatable :: hst2, hhel2, hran2, hstr2

real (long), dimension(:,:), allocatable :: en, psiVal, phiVal
real (long), dimension(:), allocatable :: protLen
real (long), dimension(:), allocatable :: sasa, emm, bond, angle, dihed, imprp, elect, vdw, bondP, angleP, &
  dihedP, imprpP, electP, vdwP, kinetic, vol, emmP
real (long), dimension(:), allocatable :: emmIon, bondIon, angleIon, dihedIon, imprpIon, electIon, vdwIon, &
  bondPIon, anglePIon, dihedPIon, imprpPIon, electPIon, vdwPIon,  emmPIon
real (long), dimension(:), allocatable :: rg1Inter, rg2Inter, rg1, rg2
real (long), dimension(:), allocatable :: tilt1_1_8, tilt1_8_14, tilt1_14_20, tilt1_1_11, tilt1_11_21
real (long), dimension(:), allocatable :: tilt2_1_8, tilt2_8_14, tilt2_14_20, tilt2_1_11, tilt2_11_21
real (long), dimension(:), allocatable :: utilt1_15_20, utiltH1_15_20, utilt1_6_14
real (long), dimension(:), allocatable :: utilt2_15_20, utiltH2_15_20, utilt2_6_14
real (long), dimension(:), allocatable :: termC_Z1, termN_Z1, termC_Z2, termN_Z2, termM_Z1, termM_Z2
integer , dimension(:), allocatable :: lrc_up, lrc_lo
integer , dimension(:), allocatable :: isHel_up, isHel_lo
real (Long) , dimension(:,:), allocatable :: lip_coord

real (long) :: pEneVal
real (long) :: roVal
real (long) :: radVal
real (long) :: sinVal
real (long) :: cosVal
real (long) :: ro_12_1_vec, ro_12_2_vec, ro_12_ins_vec, ro_12_sur_vec
real (long) :: ro_15_1_vec, ro_15_2_vec, ro_15_ins_vec, ro_15_sur_vec
real (long) :: ro_19_1_vec, ro_19_2_vec, ro_19_ins_vec, ro_19_sur_vec
real (long), dimension(0:numStates) :: clustHelC, clustHelN, clustCnt, clustCntSaddle
real (long), dimension(0:numStates) :: clustRo12_vec, clustRo15_vec, clustRo19_vec
real (long), dimension(0:numStates) :: clustRo12_sin, clustRo15_sin, clustRo19_sin
real (long), dimension(0:numStates) :: clustRo12_cos, clustRo15_cos, clustRo19_cos
real (long), dimension(0:numStates) :: clustTiltN_vec, clustTiltC_vec
real (long), dimension(0:numStates) :: clustTiltN_sin, clustTiltC_sin
real (long), dimension(0:numStates) :: clustTiltN_cos, clustTiltC_cos
integer :: clustHelCnt, currClust, currClust1, currClust2, currClustSaddle

real (long), dimension(0:21,0:70) :: zVal_2d, zVal_2d_ins, zVal_2d_sur
real (long), dimension(0:numStates,0:21,0:70) :: clustZVal_2d
real (long), dimension(0:29,0:34) :: z_ro15, z_ro19
real (long), dimension(0:34,0:34) :: nTerm_cTerm, mid_cTerm
real (long), dimension(0:29,0:34) :: z_ro15_hel, z_ro19_hel
real (long), dimension(0:20,0:20) :: lrc_hel
real (long), dimension(1:3) :: scdPgRegAvg
real (long), dimension(1:3) :: scdPcRegAvg
real (long), dimension(1:3) :: scdPgRegCntAvg
real (long), dimension(1:3) :: scdPcRegCntAvg
real (long), dimension(1:3) :: PcRegLipDenArea
real (long), dimension(1:3) :: PcRegLipDenVol
real (long), dimension(1:3) :: PgRegLipDenArea
real (long), dimension(1:3) :: PgRegLipDenVol
real (long) :: zVal_2d_sum
real (long) :: zscVal

real (long) :: ro_12_1, ro_12_2, ro_15_1, ro_15_2, ro_19_1, ro_19_2
real (long) :: ro_12_1_sin, ro_12_2_sin, ro_15_1_sin, ro_15_2_sin, ro_19_1_sin, ro_19_2_sin
real (long) :: ro_12_1_cos, ro_12_2_cos, ro_15_1_cos, ro_15_2_cos, ro_19_1_cos, ro_19_2_cos
real (long) :: ro_15_1_hel, ro_15_2_hel, ro_19_1_hel, ro_19_2_hel
real (long) :: ro_12_1_avg, ro_12_2_avg, ro_15_1_avg, ro_15_2_avg, ro_19_1_avg, ro_19_2_avg
real (long) :: tilt1_8, tilt8_14, tilt14_20, tilt1_11, tilt11_21
real (long) :: utilt15_20, utiltH15_20, utilt6_14
real (long) :: utilt15_20_ins, utiltH15_20_ins, utilt6_14_ins
real (long) :: utilt15_20_sur, utiltH15_20_sur, utilt6_14_sur
real (long) :: ro_12_sur, ro_15_sur, ro_19_sur, ro_15_sur_hel, ro_19_sur_hel
real (long) :: ro_12_sur_sin, ro_15_sur_sin, ro_19_sur_sin
real (long) :: ro_12_sur_cos, ro_15_sur_cos, ro_19_sur_cos
real (long) :: ro_12_ins, ro_15_ins, ro_19_ins, ro_15_ins_hel, ro_19_ins_hel
real (long) :: ro_12_ins_sin, ro_15_ins_sin, ro_19_ins_sin
real (long) :: ro_12_ins_cos, ro_15_ins_cos, ro_19_ins_cos
real (long) :: ro_12_sur_avg, ro_15_sur_avg, ro_19_sur_avg
real (long) :: ro_12_ins_avg, ro_15_ins_avg, ro_19_ins_avg
real (long) :: tilt1_8_ins, tilt8_14_ins, tilt14_20_ins, tilt1_11_ins, tilt11_21_ins
real (long) :: tilt1_8_sur, tilt8_14_sur, tilt14_20_sur, tilt1_11_sur, tilt11_21_sur
integer :: zBucket
integer :: nBucket, cBucket, mBucket
integer :: ro_bucket
integer :: zBuck
integer :: termInt
real (long) :: termVal

integer, dimension(0:34,0:29) :: clust_all, clust_all_saddle
integer :: val


real (long), dimension(:), allocatable :: contCnt
real (long) :: lenVal, rg1Val, rg2Val, rg1InterVal, rg2InterVal, contCntVal
integer :: lenValInt, rgValInt, contCntValInt
real (long), dimension(nrep) :: templist
real (long), dimension(repmin:repmax) :: fe, temp
real (long) :: f1, factor,tempi, pf, hi, hr, pf_ins, pf_sur, pf_ulr, pf_ulr_sur, pf_ulr_ins, isHel_pf, pf_sur_hel, pf_ins_hel
real (long), dimension(repmin:repmax) :: totUsed
real (long) :: usedTotal
real (long), dimension(ncont_leaf) :: cm_pc1, cm_pg1, cm_pc2, cm_pg2
real (long), dimension(ncont_leaf) :: cm_pc_sur, cm_pc_ins, cm_pg_sur, cm_pg_ins
real (long), dimension(0:numStates,ncont_leaf) :: clustCmPC, clustCmPG
real (long), dimension(ncont_pep) :: cm_pep1, cm_pep2
real (long) :: eppw
real (long), dimension(nmon) :: st, str, myRan, hel, totVal, zscAvg, hel_sur, hel_ins, zscAvg_sur, zscAvg_ins
real (long), dimension(nmon) :: st2, str2, myRan2, hel2, totVal2, zscAvg2
real (long), dimension(nmon-1) :: myPsi, myPhi
real (long) :: avgIns1, avgSur1, avgUnb1 
real (long) :: avgInsR1_1, avgSurR1_1, avgUnbR1_1, avgInsR2_1, avgSurR2_1, avgUnbR2_1
real (long) :: avgIns2, avgSur2, avgUnb2 
real (long) :: avgInsR1_2, avgSurR1_2, avgUnbR1_2, avgInsR2_2, avgSurR2_2, avgUnbR2_2
real (long) :: phi, psi
integer :: phiBuck, psiBuck
integer :: pEneInt


real (long) :: sasaAvg, emmAvg, &
   bondAvg, angleAvg, dihedAvg, imprpAvg, electAvg, vdwAvg, kineticAvg, volAvg, &
   bondPAvg, anglePAvg, dihedPAvg, imprpPAvg, electPAvg, vdwPAvg, emmPAvg
real (long) :: emmAvgIon, &
   bondAvgIon, angleAvgIon, dihedAvgIon, imprpAvgIon, electAvgIon, vdwAvgIon, &
   bondPAvgIon, anglePAvgIon, dihedPAvgIon, imprpPAvgIon, electPAvgIon, vdwPAvgIon, emmPAvgIon

real (long) :: rgInterAvg, rgAvg


real (long), dimension(0:21,0:160) :: helLen
real (long), dimension(0:21,0:119) :: betaLen, betaRg
real (long), dimension(0:21,0:119) :: contLen, contRg
real (long), dimension(0:21,0:21) :: betaCont
real (long), dimension(0:72,0:72) :: ramach
real (long), dimension(0:20,0:72) :: psiHist
real (long), dimension(0:20,0:72) :: phiHist
real (long), dimension(0:21) :: contHist
real (long), dimension(0:119) :: lenHist
real (long), dimension(0:119) :: rgHist
real (long), dimension(0:2999) :: pEnergyHist
real (long), dimension(0:70) :: cTermHist1
real (long), dimension(0:70) :: cTermHist2
real (long), dimension(0:70) :: nTermHist1
real (long), dimension(0:70) :: nTermHist2

real (long), dimension(1:4) :: lipid_coord

character(len=16):: trString
integer :: vals, tens



integer ::  nt1=nint(tmin), nt2=nint(tmax)
integer :: i,j,k,l,m,n,aa, itermax, nt, istr, helCnt, betaCnt, helCnt2, betaCnt2
integer :: currRep, currRepCnt
integer :: lenProb
integer :: off
real (long) :: ramachSum

 do i=0,21
   do j=0,160
     helLen(i,j) = 0.0
   enddo
   do j=0,119
     betaLen(i,j) = 0.0
     betaRg(i,j) = 0.0
     contLen(i,j) = 0.0
     contRg(i,j) = 0.0
   enddo
   do j=0,21
     betaCont(i,j) = 0.0
   enddo
 enddo

 do i=0,72
   do j=0,72
     ramach(i,j) = 0.0
   enddo
   do j=0,20
     psiHist(j,i) = 0.0
     phiHist(j,i) = 0.0
   enddo
 enddo

 do i=0,21
  do j=0, 70
    zVal_2d(i,j) = 0.0
    zVal_2d_sur(i,j) = 0.0
    zVal_2d_ins(i,j) = 0.0
    do k=0,numStates
      clustZVal_2d(k,i,j) = 0.0
    enddo
  enddo
 enddo


 do i=0,numStates
  clustCnt(i) = 0.0
  clustCntSaddle(i) = 0.0
  clustHelC(i) = 0.0
  clustHelN(i) = 0.0
  clustRo12_vec(i) = 0.0
  clustRo15_vec(i) = 0.0
  clustRo19_vec(i) = 0.0
  clustRo12_sin(i) = 0.0
  clustRo15_sin(i) = 0.0
  clustRo19_sin(i) = 0.0
  clustRo12_cos(i) = 0.0
  clustRo15_cos(i) = 0.0
  clustRo19_cos(i) = 0.0
  clustTiltN_vec(i) = 0.0
  clustTiltC_vec(i) = 0.0
  clustTiltN_sin(i) = 0.0
  clustTiltC_sin(i) = 0.0
  clustTiltN_cos(i) = 0.0
  clustTiltC_cos(i) = 0.0
  do j=1, ncont_leaf
      clustCmPC(i,j) = 0.0
      clustCmPG(i,j) = 0.0
  enddo
 enddo

 do i=0, 34
  do j=0, 34
    nTerm_cTerm(i,j) = 0.0
    mid_cTerm(i,j) = 0.0
  enddo
 enddo

 do i=0, 29
  do j=0, 34
    z_ro15(i,j) = 0.0
    z_ro19(i,j) = 0.0
    z_ro15_hel(i,j) = 0.0
    z_ro19_hel(i,j) = 0.0
  enddo
 enddo

 do i=0, 20
  do j=0, 20
    lrc_hel(i,j) = 0.0
  enddo
 enddo

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



 open(20,file='clust_saddle_hot.dat')
 do i=0,34
   read(20,*) clust_all_saddle(i,0:29)
 enddo
 close(20)

 open(20,file='clust_saddle_6Kc_hot.dat')
 do i=0,34
   read(20,*) clust_all(i,0:29)
 enddo
 close(20)



! do i=0,34
!   do j=0,29
!     val = clust_all(i,j)
!     if (val > 0) then
!        print*, "CLUST", i, j, clust_all(i,j)
!     endif
!   enddo
! enddo
! do i=0,34
!   print*,(34.5-i),'  ',clust_all(i,:)
! enddo
    

 istr=sum(nk)
 allocate( hcm_pc_up(ncont_leaf,(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pc_low(ncont_leaf,(tsim-teq)*(repmax-repmin+1)),   &
          scdPgRegCnt(3,(tsim-teq)*(repmax-repmin+1)), &
          scdPcRegCnt(3,(tsim-teq)*(repmax-repmin+1)), &
          scdPgReg(3,(tsim-teq)*(repmax-repmin+1)), &
          scdPcReg(3,(tsim-teq)*(repmax-repmin+1)), &
          regLipDen(12,(tsim-teq)*(repmax-repmin+1)), &
          hcm_pg_up(ncont_leaf,(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pg_low(ncont_leaf,(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pep1(ncont_pep,(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pep2(ncont_pep,(tsim-teq)*(repmax-repmin+1)),   &
          hstr(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hhel(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hran(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hst(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hstr2(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hhel2(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hran2(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          hst2(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          heppw(1,(tsim-teq)*(repmax-repmin+1)),   &
          protLen((tsim-teq)*(repmax-repmin+1)),   &
          rg1((tsim-teq)*(repmax-repmin+1)),   &
          rg2((tsim-teq)*(repmax-repmin+1)),   &
          rg1Inter((tsim-teq)*(repmax-repmin+1)),   &
          rg2Inter((tsim-teq)*(repmax-repmin+1)),   &
          tilt1_1_8((tsim-teq)*(repmax-repmin+1)),   &
          tilt1_8_14((tsim-teq)*(repmax-repmin+1)),   &
          tilt1_14_20((tsim-teq)*(repmax-repmin+1)),   &
          tilt1_1_11((tsim-teq)*(repmax-repmin+1)),   &
          tilt1_11_21((tsim-teq)*(repmax-repmin+1)),   &
          tilt2_1_8((tsim-teq)*(repmax-repmin+1)),   &
          tilt2_8_14((tsim-teq)*(repmax-repmin+1)),   &
          tilt2_14_20((tsim-teq)*(repmax-repmin+1)),   &
          tilt2_1_11((tsim-teq)*(repmax-repmin+1)),   &
          tilt2_11_21((tsim-teq)*(repmax-repmin+1)),   &
          utilt1_15_20((tsim-teq)*(repmax-repmin+1)),   &
          utiltH1_15_20((tsim-teq)*(repmax-repmin+1)),   &
          utilt1_6_14((tsim-teq)*(repmax-repmin+1)),   &
          utilt2_15_20((tsim-teq)*(repmax-repmin+1)),   &
          utiltH2_15_20((tsim-teq)*(repmax-repmin+1)),   &
          utilt2_6_14((tsim-teq)*(repmax-repmin+1)),   &
          phiVal(nmon-1,(tsim-teq)*(repmax-repmin+1)),   &
          psiVal(nmon-1,(tsim-teq)*(repmax-repmin+1)),   &
          sasa((tsim-teq)*(repmax-repmin+1)),     &
          emm((tsim-teq)*(repmax-repmin+1)),     &
          bond((tsim-teq)*(repmax-repmin+1)),     &
          angle((tsim-teq)*(repmax-repmin+1)),     &
          dihed((tsim-teq)*(repmax-repmin+1)),     &
          imprp((tsim-teq)*(repmax-repmin+1)),     &
          elect((tsim-teq)*(repmax-repmin+1)),     &
          vdw((tsim-teq)*(repmax-repmin+1)),     &
          vol((tsim-teq)*(repmax-repmin+1)),     &
          bondP((tsim-teq)*(repmax-repmin+1)),     &
          angleP((tsim-teq)*(repmax-repmin+1)),     &
          dihedP((tsim-teq)*(repmax-repmin+1)),     &
          imprpP((tsim-teq)*(repmax-repmin+1)),     &
          electP((tsim-teq)*(repmax-repmin+1)),     &
          vdwP((tsim-teq)*(repmax-repmin+1)),     &
          emmP((tsim-teq)*(repmax-repmin+1)),     &
          emmIon((tsim-teq)*(repmax-repmin+1)),     &
          bondIon((tsim-teq)*(repmax-repmin+1)),     &
          angleIon((tsim-teq)*(repmax-repmin+1)),     &
          dihedIon((tsim-teq)*(repmax-repmin+1)),     &
          imprpIon((tsim-teq)*(repmax-repmin+1)),     &
          electIon((tsim-teq)*(repmax-repmin+1)),     &
          vdwIon((tsim-teq)*(repmax-repmin+1)),     &
          bondPIon((tsim-teq)*(repmax-repmin+1)),     &
          anglePIon((tsim-teq)*(repmax-repmin+1)),     &
          dihedPIon((tsim-teq)*(repmax-repmin+1)),     &
          imprpPIon((tsim-teq)*(repmax-repmin+1)),     &
          electPIon((tsim-teq)*(repmax-repmin+1)),     &
          vdwPIon((tsim-teq)*(repmax-repmin+1)),     &
          emmPIon((tsim-teq)*(repmax-repmin+1)),     &
          kinetic((tsim-teq)*(repmax-repmin+1)),     &
          zsc(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          zsc2(nmon,(tsim-teq)*(repmax-repmin+1)),   &
          termC_Z1((tsim-teq)*(repmax-repmin+1)),   &
          termC_Z2((tsim-teq)*(repmax-repmin+1)),   &
          termN_Z1((tsim-teq)*(repmax-repmin+1)),   &
          termN_Z2((tsim-teq)*(repmax-repmin+1)),   &
          termM_Z1((tsim-teq)*(repmax-repmin+1)),   &
          termM_Z2((tsim-teq)*(repmax-repmin+1)),   &
          ro_12(2,(tsim-teq)*(repmax-repmin+1)),   &
          ro2_12(2,(tsim-teq)*(repmax-repmin+1)),   &
          ro_15(2,(tsim-teq)*(repmax-repmin+1)),   &
          ro2_15(2,(tsim-teq)*(repmax-repmin+1)),   &
          ro_19(2,(tsim-teq)*(repmax-repmin+1)),   &
          ro2_19(2,(tsim-teq)*(repmax-repmin+1)),   &
          lrc_up((tsim-teq)*(repmax-repmin+1)),   &
          lrc_lo((tsim-teq)*(repmax-repmin+1)),   &
          isHel_up((tsim-teq)*(repmax-repmin+1)),   &
          isHel_lo((tsim-teq)*(repmax-repmin+1)),   &
          lip_coord(4,(tsim-teq)*(repmax-repmin+1)),   &
          en(5,(tsim-teq)*(repmax-repmin+1)))

 allocate(contCnt((tsim-teq)*(repmax-repmin+1)))

!
!   Create trajectory string
 if (in_tr < 10) then
     trString = char(48+in_tr)
 else
     tens = in_tr/10
     vals = in_tr-(tens*10)
     trString = char(48+tens)//char(48+vals)
 endif


 do i=repmin,repmax
   print*,i,real(temp(i)),nk(i)  
 enddo
 

  open(34,file='scdPgReg_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)scdPgReg(:,i)
  enddo
 close(34)

  open(34,file='scdPcReg_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)scdPcReg(:,i)
  enddo
 close(34)


  open(34,file='scdPgRegCnt_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)scdPgRegCnt(:,i)
  enddo
 close(34)

  open(34,file='scdPcRegCnt_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)scdPcRegCnt(:,i)
  enddo
 close(34)

  open(34,file='regionDen_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)regLipDen(:,i)
  enddo
 close(34)


  open(34,file='histcm_pc1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hcm_pc_up(:,i)
  enddo
 close(34)


   open(34,file='LRC1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)lrc_up(i)
  enddo
  close(34)

  open(34,file='LRC2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)lrc_lo(i)
  enddo
  close(34)


     open(34,file='isHel1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)isHel_up(i)
  enddo
 close(34)

     open(34,file='isHel2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)isHel_lo(i)
  enddo
 close(34)

     open(34,file='histcm_pc2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hcm_pc_low(:,i)
  enddo
 close(34)



     open(34,file='histcm_pg1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hcm_pg_up(:,i)
  enddo
 close(34)

     open(34,file='histcm_pg2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hcm_pg_low(:,i)
  enddo
 close(34)


     open(34,file='histcm_pep1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hcm_pep1(:,i)
  enddo
 close(34)

     open(34,file='histcm_pep2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hcm_pep2(:,i)
  enddo
 close(34)

     open(34,file='histst'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hst(:,i)
  enddo
 close(34)
     open(34,file='histst2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hst2(:,i)
  enddo
 close(34)

     open(34,file='histstr'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hstr(:,i)
  enddo
 close(34)
     open(34,file='histstr2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)hstr2(:,i)
  enddo
 close(34)

     open(34,file='histhel'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)hhel(:,i)
  enddo
 close(34)
     open(34,file='histhel2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)hhel2(:,i)
  enddo
 close(34)


    open(34,file='zsc1_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)zsc(:,i)
  enddo
  close(34)

    open(34,file='zsc2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)zsc2(:,i)
  enddo
  close(34)


    open(34,file='ro12_1_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)ro_12(:,i)
  enddo
  close(34)

    open(34,file='ro12_2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)ro2_12(:,i)
  enddo
  close(34)

    open(34,file='ro15_1_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)ro_15(:,i)
  enddo
  close(34)

    open(34,file='ro15_2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)ro2_15(:,i)
  enddo
  close(34)
    open(34,file='ro19_1_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)ro_19(:,i)
  enddo
  close(34)

    open(34,file='ro19_2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)ro2_19(:,i)
  enddo
  close(34)


!     open(34,file='psiVal'//trim(trString)//'.dat',form='unformatted')
!
!  do i=1,istr
!    read(34)psiVal(:,i)
!  enddo
! close(34)
!
!     open(34,file='phiVal'//trim(trString)//'.dat',form='unformatted')
!  
!  do i=1,istr
!    read(34)phiVal(:,i)
!  enddo
! close(34)

     open(34,file='histran'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)hran(:,i)
  enddo
 close(34)
     open(34,file='histran2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)hran2(:,i)
  enddo
 close(34)

     open(34,file='sasa'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)sasa(i)
  enddo
 close(34)



     open(34,file='cTerm1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)termC_Z1(i)
  enddo
 close(34)

     open(34,file='cTerm2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)termC_Z2(i)
  enddo
 close(34)


     open(34,file='midCOM1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)termM_Z1(i)
  enddo
 close(34)

     open(34,file='midCOM2_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)termM_Z2(i)
  enddo
 close(34)

     open(34,file='nTerm1_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)termN_Z1(i)
  enddo
 close(34)

     open(34,file='nTerm2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)termN_Z2(i)
  enddo
 close(34)





  open(34,file='lipid_coord_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)lip_coord(:,i)
  enddo
 close(34)



     open(34,file='emm'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr 
    read(34)emm(i)
  enddo
 close(34)

     open(34,file='bond'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)bond(i)
  enddo
 close(34)
     open(34,file='bondP'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)bondP(i)
  enddo
 close(34)


     open(34,file='angle'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)angle(i)
  enddo
 close(34)
     open(34,file='angleP'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)angleP(i)
  enddo
 close(34)

     open(34,file='dihed'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)dihed(i)
  enddo
 close(34)
     open(34,file='dihedP'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)dihedP(i)
  enddo
 close(34)

     open(34,file='imprp'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)imprp(i)
  enddo
 close(34)
     open(34,file='imprpP'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)imprpP(i)
  enddo
 close(34)

     open(34,file='elect'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)elect(i)
  enddo
 close(34)

     open(34,file='electP'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)electP(i)
  enddo
 close(34)

     open(34,file='vdw'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)vdw(i)
  enddo
 close(34)


     open(34,file='vdwP'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)vdwP(i)
  enddo
 close(34)

     open(34,file='emmP'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)emmP(i)
  enddo
 close(34)


     open(34,file='emmIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr 
    read(34)emmIon(i)
  enddo
 close(34)

     open(34,file='bondIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)bondIon(i)
  enddo
 close(34)
     open(34,file='bondPIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)bondPIon(i)
  enddo
 close(34)


     open(34,file='angleIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)angleIon(i)
  enddo
 close(34)
     open(34,file='anglePIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)anglePIon(i)
  enddo
 close(34)

     open(34,file='dihedIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)dihedIon(i)
  enddo
 close(34)
     open(34,file='dihedPIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)dihedPIon(i)
  enddo
 close(34)

     open(34,file='imprpIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)imprpIon(i)
  enddo
 close(34)
     open(34,file='imprpPIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)imprpPIon(i)
  enddo
 close(34)

     open(34,file='electIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)electIon(i)
  enddo
 close(34)

     open(34,file='electPIon'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)electPIon(i)
  enddo
 close(34)

     open(34,file='vdwIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)vdwIon(i)
  enddo
 close(34)


     open(34,file='vdwPIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)vdwPIon(i)
  enddo
 close(34)

     open(34,file='emmPIon'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)emmPIon(i)
  enddo
 close(34)




     open(34,file='vol'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)vol(i)
  enddo
 close(34)

 flush(6)

     open(34,file='kinetic'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)kinetic(i)
  enddo
 close(34)

 
      open(34,file='energyi'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)en(:,i)
  enddo
 close(34)
 flush(6)

!     open(34,file='protLen'//trim(trString)//'.dat',form='unformatted')
!
!  do i=1,istr
!    read(34)protLen(i)
!  enddo
! close(34)


     open(34,file='rg1_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)rg1(i)
  enddo
 close(34)

     open(34,file='rg2_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)rg2(i)
  enddo
 close(34)

     open(34,file='rg1Int_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)rg1Inter(i)
  enddo
 close(34)

     open(34,file='rg2Int_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)rg2Inter(i)
  enddo
 close(34)

     open(34,file='utilt1_15_20_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)utilt1_15_20(i)
  enddo
 close(34)

!     open(34,file='utiltH1_15_20_'//trim(trString)//'.dat',form='unformatted')
!
!  do i=1,istr
!    read(34)utiltH1_15_20(i)
!  enddo
! close(34)

     open(34,file='utilt1_6_14_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)utilt1_6_14(i)
  enddo
 close(34)


     open(34,file='utilt2_15_20_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)utilt2_15_20(i)
  enddo
 close(34)

!     open(34,file='utiltH2_15_20_'//trim(trString)//'.dat',form='unformatted')
!
!  do i=1,istr
!    read(34)utiltH2_15_20(i)
!  enddo
! close(34)

     open(34,file='utilt2_6_14_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)utilt2_6_14(i)
  enddo
 close(34)

     open(34,file='tilt1_1_8_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt1_1_8(i)
  enddo
 close(34)




     open(34,file='tilt1_8_14_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt1_8_14(i)
  enddo
 close(34)


     open(34,file='tilt1_14_20_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt1_14_20(i)
  enddo
 close(34)


     open(34,file='tilt1_1_11_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt1_1_11(i)
  enddo
 close(34)



     open(34,file='tilt1_11_21_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt1_11_21(i)
  enddo
 close(34)


     open(34,file='tilt2_1_8_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt2_1_8(i)
  enddo
 close(34)

     open(34,file='tilt2_8_14_'//trim(trString)//'.dat',form='unformatted')
  do i=1,istr
    read(34)tilt2_8_14(i)
  enddo
 close(34)


     open(34,file='tilt2_14_20_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt2_14_20(i)
  enddo
 close(34)


     open(34,file='tilt2_1_11_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt2_1_11(i)
  enddo
 close(34)



     open(34,file='tilt2_11_21_'//trim(trString)//'.dat',form='unformatted')

  do i=1,istr
    read(34)tilt2_11_21(i)
  enddo
 close(34)


 do i=1,istr
  heppw(1,i)=en(1,i)+en(2,i)
 enddo

 open(36,file='ffinal'//trim(trString)//'.dat')
 read(36,*)itermax,fe
 close(36)

!     open(43,file='th_cm'//trim(trString)//'.dat')
     open(44,file='th_eppw'//trim(trString)//'.dat')
     open(45,file='th_st'//trim(trString)//'.dat')
     open(46,file='th_str'//trim(trString)//'.dat')
     open(47,file='th_hel'//trim(trString)//'.dat')
     open(48,file='th_ran'//trim(trString)//'.dat')
!     open(49,file='th_len_hel'//trim(trString)//'.dat')
!     open(50,file='th_rg_hist'//trim(trString)//'.dat')
!     open(51,file='th_len_hist'//trim(trString)//'.dat')
!     open(52,file='th_len_beta'//trim(trString)//'.dat')
!     open(53,file='th_rg_beta'//trim(trString)//'.dat')
!     open(54,file='th_rg_cont'//trim(trString)//'.dat')
!     open(55,file='th_len_cont'//trim(trString)//'.dat')
!     open(56,file='th_cont_hist'//trim(trString)//'.dat')
!     open(57,file='th_ramach'//trim(trString)//'.dat')
!     open(58,file='th_psi'//trim(trString)//'.dat')
!     open(59,file='th_phi'//trim(trString)//'.dat')
!     open(60,file='th_cont_beta'//trim(trString)//'.dat')
     open(61,file='th_p_energy_hist'//trim(trString)//'.dat')
     open(62,file='th_sasa'//trim(trString)//'.dat')
     open(63,file='th_emm'//trim(trString)//'.dat')
     open(64,file='th_bond'//trim(trString)//'.dat')
     open(65,file='th_angle'//trim(trString)//'.dat')
     open(66,file='th_dihed'//trim(trString)//'.dat')
     open(67,file='th_imprp'//trim(trString)//'.dat')
     open(68,file='th_elect'//trim(trString)//'.dat')
     open(69,file='th_vdw'//trim(trString)//'.dat')
     open(70,file='th_kinetic'//trim(trString)//'.dat')
     open(71,file='th_bondP'//trim(trString)//'.dat')
     open(72,file='th_angleP'//trim(trString)//'.dat')
     open(73,file='th_dihedP'//trim(trString)//'.dat')
     open(74,file='th_imprpP'//trim(trString)//'.dat')
     open(75,file='th_electP'//trim(trString)//'.dat')
     open(76,file='th_vdwP'//trim(trString)//'.dat')
     open(77,file='th_emmP'//trim(trString)//'.dat')
     open(78,file='th_vol'//trim(trString)//'.dat')
     open(79,file='th_rg'//trim(trString)//'.dat')
     open(80,file='th_rg_INTERM'//trim(trString)//'.dat')
     open(81,file='th_cm_pc_up'//trim(trString)//'.dat')
     open(82,file='th_cm_pg_up'//trim(trString)//'.dat')
     open(83,file='th_zsc'//trim(trString)//'.dat')
     open(84,file='th_zsc_2d'//trim(trString)//'.dat')
     open(85,file='th_cm_pc_low'//trim(trString)//'.dat')
     open(86,file='th_cm_pg_low'//trim(trString)//'.dat')
     open(87,file='th_hel_sur'//trim(trString)//'.dat')
     open(88,file='th_hel_ins'//trim(trString)//'.dat')
     open(89,file='th_cm_pc_sur'//trim(trString)//'.dat')
     open(90,file='th_cm_pc_ins'//trim(trString)//'.dat')
     open(91,file='th_cm_pg_sur'//trim(trString)//'.dat')
     open(92,file='th_cm_pg_ins'//trim(trString)//'.dat')
     open(93,file='th_cm_pep_up'//trim(trString)//'.dat')
     open(94,file='th_cm_pep_low'//trim(trString)//'.dat')
     open(95,file='th_tilt'//trim(trString)//'.dat')
     open(96,file='th_zsc_ins'//trim(trString)//'.dat')
     open(97,file='th_zsc_sur'//trim(trString)//'.dat')
     open(98,file='th_zsc_2d_ins'//trim(trString)//'.dat')
     open(99,file='th_zsc_2d_sur'//trim(trString)//'.dat')
     open(100,file='th_utilt'//trim(trString)//'.dat')
     open(101,file='th_ro'//trim(trString)//'.dat')
     open(102,file='th_term2_Z_hist'//trim(trString)//'.dat')
     open(103,file='th_term1_Z_hist'//trim(trString)//'.dat')
     open(104,file='th_ins_sur'//trim(trString)//'.dat')
     open(105,file='th_aa_insert1_'//trim(trString)//'.dat')
     open(106,file='th_aa_insert2_'//trim(trString)//'.dat')
     open(107,file='th_z_ro15_'//trim(trString)//'.dat')
     open(108,file='th_z_ro19_'//trim(trString)//'.dat')
     open(109,file='th_lrc_hel_'//trim(trString)//'.dat')
     open(110,file='th_z_ro15_hel_'//trim(trString)//'.dat')
     open(111,file='th_z_ro19_hel_'//trim(trString)//'.dat')
     open(112,file='th_ro_hel_'//trim(trString)//'.dat')
     open(113,file='th_clust_helC_'//trim(trString)//'.dat')
     open(114,file='th_clust_ro_'//trim(trString)//'.dat')
     open(115,file='th_clust_tilt_'//trim(trString)//'.dat')

     open(116,file='th_zsc_2d_ins_'//trim(trString)//'.dat')
     open(117,file='th_zsc_2d_m-ins_'//trim(trString)//'.dat')
     open(118,file='th_zsc_2d_sur_'//trim(trString)//'.dat')
     open(119,file='th_zsc_2d_m-sur1_'//trim(trString)//'.dat')

     open(120,file='th_cm_pc_clust_ins_'//trim(trString)//'.dat')
     open(121,file='th_cm_pc_clust_m-ins_'//trim(trString)//'.dat')
     open(122,file='th_cm_pc_clust_sur_'//trim(trString)//'.dat')
     open(123,file='th_cm_pc_clust_m-sur1_'//trim(trString)//'.dat')

     open(124,file='th_cm_pg_clust_ins_'//trim(trString)//'.dat')
     open(125,file='th_cm_pg_clust_m-ins_'//trim(trString)//'.dat')
     open(126,file='th_cm_pg_clust_sur_'//trim(trString)//'.dat')
     open(127,file='th_cm_pg_clust_m-sur1_'//trim(trString)//'.dat')

     open(128,file='th_clust_helN_'//trim(trString)//'.dat')

     open(129,file='th_zsc_2d_m-sur2_'//trim(trString)//'.dat')
     open(130,file='th_cm_pc_clust_m-sur2_'//trim(trString)//'.dat')
     open(131,file='th_cm_pg_clust_m-sur2_'//trim(trString)//'.dat')
     open(132,file='th_scd_reg_'//trim(trString)//'.dat')
     open(133,file='th_lip_den_reg_'//trim(trString)//'.dat')
     open(134,file='th_nTerm_zTerm_'//trim(trString)//'.dat')
     open(135,file='th_cTerm_zMid_'//trim(trString)//'.dat')
     open(136,file='th_lipid_coord_'//trim(trString)//'.dat')

     open(137,file='th_fe_ion'//trim(trString)//'.dat')


!  do nt=nt1,nt2,dnt (do temperatures)
   nt=nt1

   pf=0.d0
   pf_ins=0.d0
   pf_sur=0.d0
   pf_ins_hel=0.d0
   pf_sur_hel=0.d0
   isHel_pf= 0.d0
   pf_ulr=0.d0
   pf_ulr_ins=0.d0
   pf_ulr_sur=0.d0
   zscAvg =  0.d0
   zscAvg2 =  0.d0
   zscAvg_sur =  0.d0
   zscAvg_ins =  0.d0

   avgIns1 = 0.d0
   avgSur1 = 0.d0
   avgUnb1 = 0.d0
   avgIns2 = 0.d0
   avgSur2 = 0.d0
   avgUnb2 = 0.d0
   avgInsR1_1 = 0.d0
   avgSurR1_1 = 0.d0
   avgUnbR1_1 = 0.d0
   avgInsR1_2 = 0.d0
   avgSurR1_2 = 0.d0
   avgUnbR1_2 = 0.d0
   avgInsR2_1 = 0.d0
   avgSurR2_1 = 0.d0
   avgUnbR2_1 = 0.d0
   avgInsR2_2 = 0.d0
   avgSurR2_2 = 0.d0
   avgUnbR2_2 = 0.d0

   ro_12_sur=0.d0
   ro_12_sur_sin=0.d0
   ro_12_sur_cos=0.d0
   ro_12_ins=0.d0
   ro_12_ins_sin=0.d0
   ro_12_ins_cos=0.d0

   ro_12_sur_avg=0.d0
   ro_12_ins_avg=0.d0
   ro_12_1=0.d0
   ro_12_2=0.d0
   ro_12_1_sin=0.d0
   ro_12_2_sin=0.d0
   ro_12_1_cos=0.d0
   ro_12_2_cos=0.d0

   ro_12_1_avg=0.d0
   ro_12_2_avg=0.d0

   ro_15_sur=0.d0
   ro_15_sur_sin=0.d0
   ro_15_sur_cos=0.d0
   ro_15_ins=0.d0
   ro_15_ins_sin=0.d0
   ro_15_ins_cos=0.d0

   ro_15_sur_hel=0.d0
   ro_15_ins_hel=0.d0
   ro_15_sur_avg=0.d0
   ro_15_ins_avg=0.d0
   ro_15_1=0.d0
   ro_15_1_sin=0.d0
   ro_15_1_cos=0.d0
   ro_15_2=0.d0
   ro_15_2_sin=0.d0
   ro_15_2_cos=0.d0
   ro_15_1_avg=0.d0
   ro_15_2_avg=0.d0

   ro_19_sur=0.d0
   ro_19_sur_sin=0.d0
   ro_19_sur_cos=0.d0
   ro_19_ins=0.d0
   ro_19_ins_sin=0.d0
   ro_19_ins_cos=0.d0
   ro_19_sur_hel=0.d0
   ro_19_ins_hel=0.d0
   ro_19_sur_avg=0.d0
   ro_19_ins_avg=0.d0
   ro_19_1=0.d0
   ro_19_1_sin=0.d0
   ro_19_1_cos=0.d0
   ro_19_2=0.d0
   ro_19_2_sin=0.d0
   ro_19_2_cos=0.d0
   ro_19_1_avg=0.d0
   ro_19_2_avg=0.d0

   cm_pg_sur=0.d0
   cm_pg_ins=0.d0
   cm_pc_sur=0.d0
   cm_pc_ins=0.d0
   hel_sur = 0.d0
   hel_ins = 0.d0
   cm_pep1=0.d0
   cm_pep2=0.d0
   cm_pc1=0.d0
   cm_pc2=0.d0
   cm_pg1=0.d0
   cm_pg2=0.d0
   eppw=0.d0
   st=0.d0
   st2=0.d0
   str = 0.d0
   str2 = 0.d0
   hel = 0.d0
   hel2 = 0.d0
   myRan = 0.d0
   myRan2 = 0.d0
   helLen = 0.d0
   betaCont = 0.d0
   betaLen = 0.d0
   betaRg = 0.d0
   contLen =0.d0
   contRg = 0.d0
   rgHist = 0.d0
   lenHist = 0.d0
   contHist = 0.d0
   ramach = 0.d0
   zVal_2d = 0.d0
   clustZVal_2d = 0.d0
   clustCmPC = 0.d0
   clustCmPG = 0.d0
   zVal_2d_sur = 0.d0
   zVal_2d_ins = 0.d0
   z_ro15 = 0.d0
   nTerm_cTerm = 0.d0
   mid_cTerm = 0.d0
   z_ro19 = 0.d0
   z_ro15_hel = 0.d0
   z_ro19_hel = 0.d0
   lrc_hel = 0.d0
   psiHist = 0.d0
   phiHist = 0.d0
   sasaAvg = 0.d0
   emmAvg = 0.d0
   bondAvg = 0.d0
   angleAvg = 0.d0
   dihedAvg = 0.d0
   imprpAvg = 0.d0
   electAvg = 0.d0
   vdwAvg = 0.d0

   emmAvgIon = 0.d0
   bondAvgIon = 0.d0
   angleAvgIon = 0.d0
   dihedAvgIon = 0.d0
   imprpAvgIon = 0.d0
   electAvgIon = 0.d0
   vdwAvgIon = 0.d0

   volAvg = 0.d0
   kineticAvg = 0.d0
   pEnergyHist = 0.d0
   cTermHist1 = 0.d0
   nTermHist1 = 0.d0
   cTermHist2 = 0.d0
   nTermHist2 = 0.d0
   rgAvg = 0.d0
   rgInterAvg = 0.d0

   emmPAvg =0.d0
   bondPAvg =0.d0
   anglePAvg =0.d0
   dihedPAvg =0.d0
   imprpPAvg =0.d0
   electPAvg =0.d0
   vdwPAvg =  0.d0

   emmPAvgIon =0.d0
   bondPAvgIon =0.d0
   anglePAvgIon =0.d0
   dihedPAvgIon =0.d0
   imprpPAvgIon =0.d0
   electPAvgIon =0.d0
   vdwPAvgIon =  0.d0

   do m=1,3
       scdPgRegCntAvg(m) = 0.d0
       scdPgRegAvg(m) = 0.d0
       scdPcRegCntAvg(m) = 0.d0
       scdPcRegAvg(m) = 0.d0
       pcRegLipDenArea(m) = 0.d0
       pgRegLipDenArea(m) = 0.d0
       pcRegLipDenVol(m) = 0.d0
       pgRegLipDenVol(m) = 0.d0
   enddo
!   totPAvg = 0.d0

   tempi = real(nt,long)*RT/temp0


!....computing partition function...........
   currRep = repmin
   currRepCnt = 0
   do i = 1,istr

     f1 = 0.d0
     hi = en(1,i) + en(2,i) + en(3,i) + p*en(5,i) + edelta
     do m=repmin,repmax
      hr = en(1,i) + &
          sqrt(restScale(m))*en(2,i) + &
          (restScale(m))*en(3,i) + &
          p*en(5,i) + edelta
      f1 = f1 + real(nk(m),long)*exp( fe(m) - hr/temp(m)+hi/tempi)
     enddo
     if(f1 > huge(real(100,long)))then
       print*,'f1=inf, f1, i and m=',f1,i,m !; stop
     endif
     pf = pf + 1.d0/f1


     do m=1,3
         scdPgRegCntAvg(m) = scdPgRegCntAvg(m) + ((1.0/f1) * scdPgRegCnt(m,i))
         scdPgRegAvg(m) = scdPgRegAvg(m)  + ((1.0/f1) * scdPgRegCnt(m,i) * scdPgReg(m,i))
         scdPcRegCntAvg(m) = scdPcRegCntAvg(m) + ((1.0/f1) * scdPcRegCnt(m,i))
         scdPcRegAvg(m) = scdPcRegAvg(m)  + ((1.0/f1) * scdPcRegCnt(m,i) * scdPcReg(m,i))
     enddo
     lipid_coord(:) = lipid_coord(:) + (real(lip_coord(:,i),long) / f1)
     PcRegLipDenVol(1) =  PcRegLipDenVol(1) + ((1.0/f1) *  regLipDen(1,i))
     PgRegLipDenVol(1) =  PgRegLipDenVol(1) + ((1.0/f1) *  regLipDen(2,i))
     PcRegLipDenArea(1) = PcRegLipDenArea(1) + ((1.0/f1) * regLipDen(3,i))
     PgRegLipDenArea(1) = PgRegLipDenArea(1) + ((1.0/f1) * regLipDen(4,i))

     PcRegLipDenVol(2) =  PcRegLipDenVol(2) + ((1.0/f1) *  regLipDen(5,i))
     PgRegLipDenVol(2) =  PgRegLipDenVol(2) + ((1.0/f1) *  regLipDen(6,i))
     PcRegLipDenArea(2) = PcRegLipDenArea(2) + ((1.0/f1) * regLipDen(7,i))
     PgRegLipDenArea(2) = PgRegLipDenArea(2) + ((1.0/f1) * regLipDen(8,i))

     PcRegLipDenVol(3) =  PcRegLipDenVol(3) + ((1.0/f1) *  regLipDen(9,i))
     PgRegLipDenVol(3) =  PgRegLipDenVol(3) + ((1.0/f1) *  regLipDen(10,i))
     PcRegLipDenArea(3) = PcRegLipDenArea(3) + ((1.0/f1) * regLipDen(11,i))
     PgRegLipDenArea(3) = PgRegLipDenArea(3) + ((1.0/f1) * regLipDen(12,i))
     pEneVal = emmP(i)
     pEneInt = int(pEneVal)  - 7000
     if (pEneInt < 0) then
         pEneInt = 0
     endif
     if (pEneInt >= 3000) then
         pEneInt = 2999
     endif
     pEnergyHist(pEneInt) = pEnergyHist(pEneInt) + (1.0/f1)

     termVal = termC_Z1(i)
     termInt = int(termVal) + 35
     if (termInt < 0) then
         termInt = 0
     endif
     if (termInt >= 70) then
         termInt = 69
     endif
     cTermHist1(TermInt) = cTermHist1(TermInt) + (1.0/f1)

     termVal = termC_Z2(i)
     termInt = int(termVal) + 35
     if (termInt < 0) then
         termInt = 0
     endif
     if (termInt >= 70) then
         termInt = 69
     endif
     cTermHist2(termInt) = cTermHist2(termInt) + (1.0/f1)

     termVal = termN_Z1(i)
     termInt = int(termVal) + 35
     if (termInt < 0) then
         termInt = 0
     endif
     if (termInt >= 70) then
         termInt = 69
     endif
     nTermHist1(TermInt) = nTermHist1(TermInt) + (1.0/f1)

     termVal = termN_Z2(i)
     termInt = int(termVal) + 35
     if (termInt < 0) then
         termInt = 0
     endif
     if (termInt >= 70) then
         termInt = 69
     endif
     nTermHist2(termInt) = nTermHist2(termInt) + (1.0/f1)

     totVal(:) = hst(:,i) + hstr(:,i) + hran(:,i) + hhel(:,i)
     totVal2(:) = hst2(:,i) + hstr2(:,i) + hran2(:,i) + hhel2(:,i)

!
!  Hist of Ro and Z
!  Upper leaflet
     nBucket = int(termN_Z1(i))
     if (nBucket < 0) then
       nBucket = 0
       print*, "ERROR - nZ is < 0 UP ", i, termN_Z1(i)
     endif
     if (nBucket >= 35) then
       nBucket = 34
       print*, "ERROR - nZ is > 35 UP ", i, termN_Z1(i)
     endif

     cBucket = int(termC_Z1(i))
     if (cBucket < 0) then
       cBucket = 0
       print*, "ERROR - cZ is < 0 UP ", i, termC_Z1(i)
     endif
     if (cBucket >= 35) then
       cBucket = 34
       print*, "ERROR - cZ is > 35 UP ", i, termC_Z1(i)
     endif


     mBucket = int(termM_Z1(i))
     if (mBucket < 0) then
       mBucket = 0
       print*, "ERROR - mZ is < 0 UP ", i, termM_Z1(i)
     endif
     if (mBucket >= 35) then
       mBucket = 34
       print*, "ERROR - mZ is > 35 UP ", i, termM_Z1(i)
     endif


     mid_cTerm(mBucket, cBucket) =  mid_cTerm(mBucket, cBucket)  + (1.0/f1)
     nTerm_cTerm(nBucket, cBucket) =  nTerm_cTerm(nBucket, cBucket)  + (1.0/f1)


     zBucket = int(termC_Z1(i))
     if (zBucket < 0) then
       zBucket = 0
       print*, "ERROR - Z is < 0 UP ", i, termC_Z1(i)
     endif
     if (zBucket >= 35) then
       zBucket = 34
       print*, "ERROR - Z is > 35 UP ", i, termC_Z1(i)
     endif



     ro_bucket = int(ro_15(1,i)/12)
     if (ro_bucket < 0) then
       ro_bucket = 0
       print*, "ERROR - ro15 is < 0 UP ", i, ro_15(1,i)
     endif
     if (ro_bucket >= 30) then
       ro_bucket = 29
       print*, "ERROR - ro15 is > 360 UP ", i, ro_15(1,i)
     endif
     z_ro15(ro_bucket, zBucket) = z_ro15(ro_bucket, zBucket)  + (1.0/f1)
     if (isHel_up(i) == 1) then 
         z_ro15_hel(ro_bucket, zBucket) = z_ro15(ro_bucket, zBucket)  + (1.0/f1)
         isHel_pf = isHel_pf + (1.0/f1)
     endif

     ro_bucket = int(ro_19(1,i)/12)
     if (ro_bucket < 0) then
       ro_bucket = 0
       print*, "ERROR - ro19 is < 0 UP ", i, ro_19(1,i)
     endif
     if (ro_bucket >= 30) then
       ro_bucket = 29
       print*, "ERROR - ro19 is > 360 UP ", i, ro_19(1,i)
     endif
     z_ro19(ro_bucket, zBucket) = z_ro19(ro_bucket, zBucket)  + (1.0/f1)
     if (isHel_up(i) == 1) then 
         z_ro19_hel(ro_bucket, zBucket) = z_ro19(ro_bucket, zBucket)  + (1.0/f1)
     endif
     currClust = clust_all(zBucket, ro_bucket)
     currClustSaddle = clust_all_saddle(zBucket, ro_bucket)
     clustCntSaddle(currClustSaddle) = clustCntSaddle(currClustSaddle) + (1.0/f1)
     currClust1 = currClust
     if (currClust > 0) then
       clustHelCnt = 0
       do j=1,11
         if (hhel(j,i) == 1) then
           clustHelCnt = clustHelCnt + 1
         endif
       enddo
       clustHelN(currClust) = clustHelN(currClust) + (clustHelCnt/f1)
       clustHelCnt = 0
       do j=12,21
         if (hhel(j,i) == 1) then
           clustHelCnt = clustHelCnt + 1
         endif
       enddo
       clustHelC(currClust) = clustHelC(currClust) + (clustHelCnt/f1)

       radVal = ro_12(1,i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustRo12_sin(currClust) = clustRo12_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustRo12_cos(currClust) = clustRo12_cos(currClust) + real(cosVal,long)/f1

       radVal = ro_15(1,i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustRo15_sin(currClust) = clustRo15_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustRo15_cos(currClust) = clustRo15_cos(currClust) + real(cosVal,long)/f1

       radVal = ro_19(1,i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustRo19_sin(currClust) = clustRo19_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustRo19_cos(currClust) = clustRo19_cos(currClust) + real(cosVal,long)/f1

       radVal = utilt1_6_14(i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustTiltN_sin(currClust) = clustTiltN_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustTiltN_cos(currClust) = clustTiltN_cos(currClust) + real(cosVal,long)/f1


       radVal = utilt1_15_20(i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustTiltC_sin(currClust) = clustTiltC_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustTiltC_cos(currClust) = clustTiltC_cos(currClust) + real(cosVal,long)/f1
     endif
     clustCnt(currClust) = clustCnt(currClust) + (1.0/f1)
!
!  Now Z ro for lower leaflet
!

     nBucket = -int(termN_Z2(i))
     if (nBucket < 0) then
       nBucket = 0
       print*, "ERROR - nZ is < -0 LO ", i, termN_Z2(i)
     endif
     if (nBucket >= 35) then
       nBucket = 34
       print*, "ERROR - nZ is > 35 LO ", i, termN_Z2(i)
     endif

     cBucket = -int(termC_Z2(i))
     if (cBucket < 0) then
       cBucket = 0
       print*, "ERROR - cZ is < 0 LO ", i, termC_Z2(i)
     endif
     if (cBucket >= 35) then
       cBucket = 34
       print*, "ERROR - cZ is > 35 LO ", i, termC_Z2(i)
     endif


     mBucket = -int(termM_Z2(i))
     if (mBucket < 0) then
       mBucket = 0
       print*, "ERROR - mZ is < 0 UP ", i, termM_Z2(i)
     endif
     if (mBucket >= 35) then
       mBucket = 34
       print*, "ERROR - mZ is > 35 UP ", i, termM_Z2(i)
     endif


     mid_cTerm(mBucket, cBucket) =  mid_cTerm(mBucket, cBucket)  + (1.0/f1)
     nTerm_cTerm(nBucket, cBucket) =  nTerm_cTerm(nBucket, cBucket)  + (1.0/f1)

     zBucket = -int(termC_Z2(i))
     if (zBucket < 0) then
       zBucket = 0
       print*, "ERROR - Z is < 0 LO ", i, termC_Z2(i)
     endif
     if (zBucket >= 35) then
       zBucket = 34
       print*, "ERROR - Z is > 35 LO ", i, termC_Z2(i)
     endif

     ro_bucket = int(ro2_15(1,i)/12)
     if (ro_bucket < 0) then
       ro_bucket = 0
       print*, "ERROR - ro15 is < 0 LO ", i, ro_15(2,i)
     endif
     if (ro_bucket >= 30) then
       ro_bucket = 29
       print*, "ERROR - ro15 is > 360 LO ", i, ro_15(2,i)
     endif
     z_ro15(ro_bucket, zBucket) = z_ro15(ro_bucket, zBucket)  + (1.0/f1)
     if (isHel_lo(i) == 1) then 
         z_ro15_hel(ro_bucket, zBucket) = z_ro15(ro_bucket, zBucket)  + (1.0/f1)
         isHel_pf = isHel_pf + (1.0/f1)
     endif

     ro_bucket = int(ro2_19(1,i)/12)
     if (ro_bucket < 0) then
       ro_bucket = 0
       print*, "ERROR - ro19 is < 0 LO ", i, ro_19(2,i)
     endif
     if (ro_bucket >= 30) then
       ro_bucket = 29
       print*, "ERROR - ro19 is > 360 LO ", i, ro_19(2,i)
     endif
     z_ro19(ro_bucket, zBucket) = z_ro19(ro_bucket, zBucket)  + (1.0/f1)
     if (isHel_lo(i) == 1) then 
         z_ro19_hel(ro_bucket, zBucket) = z_ro19(ro_bucket, zBucket)  + (1.0/f1)
     endif
     currClust = clust_all(zBucket, ro_bucket)
     currClust2 = currClust
     currClustSaddle = clust_all_saddle(zBucket, ro_bucket)
     clustCntSaddle(currClustSaddle) = clustCntSaddle(currClustSaddle) + (1.0/f1)

     if (currClust > 0) then
       clustHelCnt = 0
       do j=1,11
         if (hhel2(j,i) == 1) then
           clustHelCnt = clustHelCnt + 1
         endif
       enddo
       clustHelN(currClust) = clustHelN(currClust) + (clustHelCnt/f1)

       clustHelCnt = 0
       do j=12,21
         if (hhel2(j,i) == 1) then
           clustHelCnt = clustHelCnt + 1
         endif
       enddo
       clustHelC(currClust) = clustHelC(currClust) + (clustHelCnt/f1)

       radVal = ro2_12(1,i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustRo12_sin(currClust) = clustRo12_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustRo12_cos(currClust) = clustRo12_cos(currClust) + real(cosVal,long)/f1

       radVal = ro2_15(1,i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustRo15_sin(currClust) = clustRo15_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustRo15_cos(currClust) = clustRo15_cos(currClust) + real(cosVal,long)/f1

       radVal = ro2_19(1,i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustRo19_sin(currClust) = clustRo19_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustRo19_cos(currClust) = clustRo19_cos(currClust) + real(cosVal,long)/f1

       radVal = utilt2_6_14(i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustTiltN_sin(currClust) = clustTiltN_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustTiltN_cos(currClust) = clustTiltN_cos(currClust) + real(cosVal,long)/f1

       radVal = utilt2_15_20(i) / (180.0 / M_PI)
       sinVal = sin(radVal)
       clustTiltC_sin(currClust) = clustTiltC_sin(currClust) + real(sinVal,long)/f1
       cosVal = cos(radVal)
       clustTiltC_cos(currClust) = clustTiltC_cos(currClust) + real(cosVal,long)/f1

     endif
     clustCnt(currClust) = clustCnt(currClust) + (1.0/f1)


!
!  First structure
     if (i == 1) then
! 
!  Upper leaf surface
         if (termC_Z1(i) > z21_clust) then
             pf_sur = 1.d0/f1
             pf_ins = 0.d0
             cm_pc_sur(:) = real(hcm_pc_up(:,i), long)/f1
             cm_pc_ins = 0.d0
             cm_pg_sur(:) = real(hcm_pg_up(:,i), long)/f1
             cm_pg_ins = 0.d0
             hel_sur(:) = real(hhel(:,i), long)/f1
             hel_ins = 0.d0
             tilt1_8_sur = real(tilt1_1_8(i), long)/f1
             tilt8_14_sur = real(tilt1_8_14(i), long)/f1
             tilt14_20_sur = real(tilt1_14_20(i), long)/f1
             tilt1_11_sur = real(tilt1_1_11(i), long)/f1
             tilt11_21_sur = real(tilt1_11_21(i), long)/f1
             tilt1_8_ins = 0.d0
             tilt8_14_ins = 0.d0
             tilt14_20_ins = 0.d0
             tilt1_11_ins = 0.d0
             tilt11_21_ins = 0.d0
             ro_12_sur = real(ro_12(1,i), long)/f1
             radVal = ro_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_sur_sin = real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_sur_cos = real(cosVal,long)/f1

             ro_12_sur_avg = real(ro_12(2,i), long)/f1
             ro_15_sur = real(ro_15(1,i), long)/f1
             radVal = ro_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_sur_sin = real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_sur_cos = real(cosVal,long)/f1

             ro_15_sur_avg = real(ro_15(2,i), long)/f1
             ro_19_sur = real(ro_19(1,i), long)/f1
             radVal = ro_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_sur_sin = real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_sur_cos = real(cosVal,long)/f1
             ro_19_sur_avg = real(ro_19(2,i), long)/f1
             ro_12_ins = 0.d0
             ro_12_ins_sin = 0.d0
             ro_12_ins_cos = 0.d0
             ro_12_ins_avg = 0.d0
             ro_15_ins = 0.d0
             ro_15_ins_sin = 0.d0
             ro_15_ins_cos = 0.d0
             ro_15_ins_avg = 0.d0
             ro_19_ins = 0.d0
             ro_19_ins_sin = 0.d0
             ro_19_ins_cos = 0.d0
             ro_19_ins_avg = 0.d0
             if (isHel_up(i) == 1) then
                 ro_15_sur_hel = real(ro_15(1,i), long)/f1
                 ro_19_sur_hel = real(ro_19(1,i), long)/f1
                 pf_sur_hel = 1.d0/f1
                 pf_ins_hel = 0.d0
             else
                 ro_15_sur_hel = 0.d0
                 ro_19_sur_hel = 0.d0
                 ro_15_ins_hel = 0.d0
                 ro_19_ins_hel = 0.d0
                 pf_sur_hel = 0.d0
                 pf_ins_hel = 0.d0
             endif
                 

             utilt15_20_sur = real(utilt1_15_20(i), long)/f1
             utilt15_20_ins = 0.d0
             utilt6_14_sur = real(utilt1_6_14(i), long)/f1
             utilt6_14_ins = 0.d0
             if (utiltH1_15_20(i) > 0.0) then
                 utiltH15_20_sur = real(utiltH1_15_20(i), long)/f1
                 utiltH15_20 = real(utiltH1_15_20(i), long)/f1
                 utiltH15_20_ins = 0.d0
                 pf_ulr_sur = 1.d0/f1
                 pf_ulr = 1.d0/f1
                 pf_ulr_ins = 0.d0
             else
                 utiltH15_20_sur = 0.d0
                 utiltH15_20_ins = 0.d0
                 pf_ulr_sur = 0.d0
                 pf_ulr = 0.d0
                 pf_ulr_ins = 0.d0
             endif

             zscAvg_sur = real(zsc(:,i), long)/f1
             zscAvg_ins = 0.d0
! 
!  First struct
!  Upper leaf inserted
         else
             pf_ins = 1.d0/f1
             pf_sur = 0.d0
             cm_pc_ins(:) = real(hcm_pc_up(:,i), long)/f1
             cm_pc_sur = 0.d0
             cm_pg_ins(:) = real(hcm_pg_up(:,i), long)/f1
             cm_pg_sur = 0.d0
             hel_ins(:) = real(hhel(:,i), long)/f1
             hel_sur = 0.d0
             tilt1_8_ins = real(tilt1_1_8(i), long)/f1
             tilt8_14_ins = real(tilt1_8_14(i), long)/f1
             tilt14_20_ins = real(tilt1_14_20(i), long)/f1
             tilt1_11_ins = real(tilt1_1_11(i), long)/f1
             tilt11_21_ins = real(tilt1_11_21(i), long)/f1
             tilt1_8_sur = 0.d0
             tilt8_14_sur = 0.d0
             tilt14_20_sur = 0.d0
             tilt1_11_sur = 0.d0
             tilt11_21_sur = 0.d0

             ro_12_ins = real(ro_12(1,i), long)/f1
             radVal = ro_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_ins_sin = real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_ins_cos = real(cosVal,long)/f1

             ro_12_ins_avg = real(ro_12(2,i), long)/f1
             ro_15_ins = real(ro_15(1,i), long)/f1
             radVal = ro_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_ins_sin = real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_ins_cos = real(cosVal,long)/f1

             ro_15_ins_avg = real(ro_15(2,i), long)/f1
             ro_19_ins = real(ro_19(1,i), long)/f1
             radVal = ro_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_ins_sin = real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_ins_cos = real(cosVal,long)/f1
             ro_19_ins_avg = real(ro_19(2,i), long)/f1
             ro_12_sur = 0.d0
             ro_12_sur_sin = 0.d0
             ro_12_sur_cos = 0.d0
             ro_12_sur_avg = 0.d0
             ro_15_sur = 0.d0
             ro_15_sur_sin = 0.d0
             ro_15_sur_cos = 0.d0
             ro_15_sur_avg = 0.d0
             ro_19_sur = 0.d0
             ro_19_sur_sin = 0.d0
             ro_19_sur_cos = 0.d0
             ro_19_sur_avg = 0.d0
             if (isHel_up(i) == 1) then
                 ro_15_ins_hel = real(ro_15(1,i), long)/f1
                 ro_19_ins_hel = real(ro_19(1,i), long)/f1
                 pf_ins_hel = 1.d0/f1
                 pf_sur_hel = 0.d0
             else
                 ro_15_sur_hel = 0.d0
                 ro_19_sur_hel = 0.d0
                 ro_15_ins_hel = 0.d0
                 ro_19_ins_hel = 0.d0
             endif

             utilt15_20_ins = real(utilt1_15_20(i), long)/f1
             utilt15_20_sur = 0.d0
             utilt6_14_ins = real(utilt1_6_14(i), long)/f1
             utilt6_14_sur = 0.d0
             if (utiltH1_15_20(i) > 0.0) then
                 utiltH15_20_ins = real(utiltH1_15_20(i), long)/f1
                 utiltH15_20 = real(utiltH1_15_20(i), long)/f1
                 utiltH15_20_sur = 0.d0
                 pf_ulr_ins = 1.d0/f1
                 pf_ulr = 1.d0/f1
                 pf_ulr_sur = 0.d0
             else
                 utiltH15_20_sur = 0.d0
                 utiltH15_20_ins = 0.d0
                 pf_ulr_sur = 0.d0
                 pf_ulr = 0.d0
                 pf_ulr_ins = 0.d0
             endif

             zscAvg_ins = real(zsc(:,i), long)/f1
             zscAvg_sur = 0.d0
         endif

         do aa=1,11
             if (zsc(aa,i) < zSurface) then
                 avgIns1 = avgIns1 + (1.0/f1)
                 avgInsR1_1 = avgInsR1_1 + (1.0/f1)
             else if (zsc(aa,i) < zUnbound) then
                 avgSur1 = avgSur1 + (1.0/f1)
                 avgSurR1_1 = avgSurR1_1 + (1.0/f1)
             else
                 avgUnb1 = avgUnb1 + (1.0/f1)
                 avgUnbR1_1 = avgUnbR1_1 + (1.0/f1)
             endif
 
            if (zsc2(aa,i) < zSurface) then
                 avgIns2 = avgIns2 + (1.0/f1)
                 avgInsR1_2 = avgInsR1_2 + (1.0/f1)
             else if (zsc2(aa,i) < zUnbound) then
                 avgSur2 = avgSur2 + (1.0/f1)
                 avgSurR1_2 = avgSurR1_2 + (1.0/f1)
             else
                 avgUnb2 = avgUnb2 + (1.0/f1)
                 avgUnbR1_2 = avgUnbR1_2 + (1.0/f1)
             endif
         enddo

        do aa=12,21
             if (zsc(aa,i) < zSurface) then
                 avgIns1 = avgIns1 + (1.0/f1)
                 avgInsR2_1 = avgInsR2_1 + (1.0/f1)
             else if (zsc(aa,i) < zUnbound) then
                 avgSur1 = avgSur1 + (1.0/f1)
                 avgSurR2_1 = avgSurR2_1 + (1.0/f1)
             else
                 avgUnb1 = avgUnb1 + (1.0/f1)
                 avgUnbR2_1 = avgUnbR2_1 + (1.0/f1)
             endif

            if (zsc2(aa,i) < zSurface) then
                 avgIns2 = avgIns2 + (1.0/f1)
                 avgInsR2_2 = avgInsR2_2 + (1.0/f1)
             else if (zsc2(aa,i) < zUnbound) then
                 avgSur2 = avgSur2 + (1.0/f1)
                 avgSurR2_2 = avgSurR2_2 + (1.0/f1)
             else
                 avgUnb2 = avgUnb2 + (1.0/f1)
                 avgUnbR2_2 = avgUnbR2_2 + (1.0/f1)
             endif
         enddo

!
!  First struct
!  lower leaf
!  Surface
         if (-termC_Z2(i) > z21_clust) then
             pf_sur = pf_sur + 1.d0/f1
             cm_pc_sur(:) = cm_pc_sur(:) + real(hcm_pc_low(:,i), long)/f1
             cm_pg_sur(:) = cm_pg_sur(:) + real(hcm_pg_low(:,i), long)/f1
             hel_sur(:) = hel_sur(:) + real(hhel2(:,i), long)/f1
             tilt1_8_sur = tilt1_8_sur + real(tilt2_1_8(i), long)/f1
             tilt8_14_sur = tilt8_14_sur + real(tilt2_8_14(i), long)/f1
             tilt14_20_sur = tilt14_20_sur + real(tilt2_14_20(i), long)/f1
             tilt1_11_sur = tilt1_11_sur + real(tilt2_1_11(i), long)/f1
             tilt11_21_sur = tilt11_21_sur + real(tilt2_11_21(i), long)/f1
             zscAvg_sur = zscAvg_sur + real(zsc2(:,i), long)/f1

             roVal = ro2_12(1,i)
             ro_12_sur = ro_12_sur + real(roVal, long)/f1
             radVal = ro2_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_sur_sin = ro_12_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_sur_cos = ro_12_sur_cos + real(cosVal,long)/f1

             roVal = ro2_15(1,i)
             ro_15_sur = ro_15_sur + real(roVal, long)/f1
             radVal = ro2_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_sur_sin = ro_15_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_sur_cos = ro_15_sur_cos + real(cosVal,long)/f1

             roVal = ro2_19(1,i)
             radVal = ro2_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_sur_sin = ro_19_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_sur_cos = ro_19_sur_cos + real(cosVal,long)/f1
             ro_19_sur = ro_19_sur + real(roVal, long)/f1

             roVal = ro2_12(2,i)
             ro_12_sur_avg = ro_12_sur_avg + real(roVal, long)/f1

             roVal = ro2_15(2,i)
             ro_15_sur_avg = ro_15_sur_avg + real(roVal, long)/f1

             roVal = ro2_19(2,i)
             ro_19_sur_avg = ro_19_sur_avg + real(roVal, long)/f1

             if (isHel_lo(i) == 1) then
                 ro_15_sur_hel = ro_15_sur_hel + real(ro2_15(1,i), long)/f1
                 ro_19_sur_hel = ro_19_sur_hel + real(ro2_19(1,i), long)/f1
                 pf_sur_hel = pf_sur_hel + 1.d0/f1
             endif
!
!  First struct
!  lower leaf
!  inserted
         else
             pf_ins = pf_ins + 1.d0/f1
             cm_pc_ins(:) = cm_pc_ins(:) + real(hcm_pc_low(:,i), long)/f1
             cm_pg_ins(:) = cm_pg_ins(:) + real(hcm_pg_low(:,i), long)/f1
             hel_ins(:) = hel_ins(:) + real(hhel2(:,i), long)/f1
             tilt1_8_ins = tilt1_8_ins + real(tilt2_1_8(i), long)/f1
             tilt8_14_ins = tilt8_14_ins + real(tilt2_8_14(i), long)/f1
             tilt14_20_ins = tilt14_20_ins + real(tilt2_14_20(i), long)/f1
             tilt1_11_ins = tilt1_11_ins + real(tilt2_1_11(i), long)/f1
             tilt11_21_ins = tilt11_21_ins + real(tilt2_11_21(i), long)/f1
             zscAvg_ins = zscAvg_ins + real(zsc2(:,i), long)/f1

             roVal = ro2_12(1,i)
             ro_12_ins = ro_12_ins + real(roVal, long)/f1
             radVal = ro2_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_ins_sin = ro_12_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_ins_cos = ro_12_ins_cos + real(cosVal,long)/f1


             roVal = ro2_15(1,i)
             ro_15_ins = ro_15_ins + real(roVal, long)/f1
             radVal = ro2_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_ins_sin = ro_15_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_ins_cos = ro_15_ins_cos + real(cosVal,long)/f1

             roVal = ro2_19(1,i)
             ro_19_ins = ro_19_ins + real(roVal, long)/f1
             radVal = ro2_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_ins_sin = ro_19_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_ins_cos = ro_19_ins_cos + real(cosVal,long)/f1

             roVal = ro2_12(2,i)
             ro_12_ins_avg = ro_12_ins_avg + real(roVal, long)/f1
             roVal = ro2_15(2,i)
             ro_15_ins_avg = ro_15_ins_avg + real(roVal, long)/f1
             roVal = ro2_19(2,i)
             ro_19_ins_avg = ro_19_ins_avg + real(roVal, long)/f1

             if (isHel_lo(i) == 1) then
                 ro_15_ins_hel = ro_15_ins_hel + real(ro2_15(1,i), long)/f1
                 ro_19_ins_hel = ro_19_ins_hel + real(ro2_19(1,i), long)/f1
                 pf_ins_hel = pf_ins_hel + 1.d0/f1
             endif

         endif
         ro_12_1 = real(ro_12(1,i), long)/f1
         radVal = ro_12(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_12_1_sin = real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_12_1_cos = real(cosVal,long)/f1

         ro_12_1_avg = real(ro_12(2,i), long)/f1

         ro_15_1 = real(ro_15(1,i), long)/f1
         radVal = ro_15(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_15_1_sin = real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_15_1_cos = real(cosVal,long)/f1
         ro_15_1_avg = real(ro_15(2,i), long)/f1

         ro_19_1 = real(ro_19(1,i), long)/f1
         ro_19_1_avg = real(ro_19(2,i), long)/f1
         radVal = ro_19(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_19_1_sin = real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_19_1_cos = real(cosVal,long)/f1

         ro_12_2 = real(ro2_12(1,i), long)/f1
         radVal = ro2_12(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_12_2_sin = real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_12_2_cos = real(cosVal,long)/f1
         ro_12_2_avg = real(ro2_12(2,i), long)/f1

         ro_15_2 = real(ro2_15(1,i), long)/f1
         radVal = ro2_15(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_15_2_sin = real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_15_2_cos = real(cosVal,long)/f1
         ro_15_2_avg = real(ro2_15(2,i), long)/f1

         ro_19_2 = real(ro2_19(1,i), long)/f1
         radVal = ro2_19(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_19_2_sin = real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_19_2_cos = real(cosVal,long)/f1
         ro_19_2_avg = real(ro2_19(2,i), long)/f1

         zscAvg(:) = real(zsc(:,i), long)/f1
         zscAvg2(:) = real(zsc2(:,i), long)/f1
         cm_pc1(:) = real(hcm_pc_up(:,i), long)/f1
         cm_pc2(:) = real(hcm_pc_low(:,i), long)/f1
         cm_pg1(:) = real(hcm_pg_up(:,i), long)/f1
         cm_pg2(:) = real(hcm_pg_low(:,i), long)/f1
         cm_pep1(:) = real(hcm_pep1(:,i), long)/f1
         cm_pep2(:) = real(hcm_pep2(:,i), long)/f1
         eppw = real(heppw(1,i), long)/f1
         st(:) = real(hst(:,i), long)/f1
         str(:) = real(hstr(:,i), long)/f1
         hel(:) = real(hhel(:,i), long)/f1
         myRan(:) = real(hran(:,i), long)/f1
         st2(:) = real(hst2(:,i), long)/f1
         str2(:) = real(hstr2(:,i), long)/f1
         hel2(:) = real(hhel2(:,i), long)/f1
         myRan2(:) = real(hran2(:,i), long)/f1
         myPsi(:) = real(psiVal(:,i), long)/f1
         myPhi(:) = real(phiVal(:,i), long)/f1
         sasaAvg = real(sasa(i), long)/f1
         emmAvg = real(emm(i), long)/f1
         bondAvg = real(bond(i), long)/f1
         angleAvg = real(angle(i), long)/f1
         dihedAvg = real(dihed(i), long)/f1
         imprpAvg = real(imprp(i), long)/f1
         electAvg = real(elect(i), long)/f1
         vdwAvg = real(vdw(i), long)/f1

         emmPAvg = real(emm(i), long)/f1
         bondPAvg = real(bondP(i), long)/f1
         anglePAvg = real(angleP(i), long)/f1
         dihedPAvg = real(dihedP(i), long)/f1
         imprpPAvg = real(imprpP(i), long)/f1
         electPAvg = real(electP(i), long)/f1
         vdwPAvg = real(vdwP(i), long)/f1

         emmAvgIon = real(emm(i), long)/f1
         bondAvgIon = real(bond(i), long)/f1
         angleAvgIon = real(angle(i), long)/f1
         dihedAvgIon = real(dihed(i), long)/f1
         imprpAvgIon = real(imprp(i), long)/f1
         electAvgIon = real(elect(i), long)/f1
         vdwAvgIon = real(vdw(i), long)/f1

         emmPAvgIon = real(emm(i), long)/f1
         bondPAvgIon = real(bondP(i), long)/f1
         anglePAvgIon = real(angleP(i), long)/f1
         dihedPAvgIon = real(dihedP(i), long)/f1
         imprpPAvgIon = real(imprpP(i), long)/f1
         electPAvgIon = real(electP(i), long)/f1
         vdwPAvgIon = real(vdwP(i), long)/f1

         kineticAvg = real(kinetic(i), long)/f1
         volAvg = en(5,i)
         rgAvg = real((rg1(i) + rg2(i)), long)/f1
         rgInterAvg = real((rg1Inter(i) + rg2Inter(i)), long)/f1
         tilt1_8 = real((tilt1_1_8(i) + tilt2_1_8(i)), long)/f1
         tilt8_14 = real((tilt1_8_14(i) + tilt2_8_14(i)), long)/f1
         tilt14_20 = real((tilt1_14_20(i) + tilt2_14_20(i)), long)/f1
         tilt1_11 = real((tilt1_1_11(i) + tilt2_1_11(i)), long)/f1
         tilt11_21 = real((tilt1_11_21(i) + tilt2_11_21(i)), long)/f1
!         print*, 'TILT ', i, tilt1_8, tilt8_14, tilt14_20, tilt1_11, tilt11_21
     else
!
!  Not first struct
!  Upper leaf
!  Surface
        if (termC_Z1(i) > z21_clust) then
             pf_sur = pf_sur + 1.d0/f1
             cm_pc_sur(:) = cm_pc_sur(:) + real(hcm_pc_up(:,i), long)/f1
             cm_pg_sur(:) = cm_pg_sur(:) + real(hcm_pg_up(:,i), long)/f1
             hel_sur(:) = hel_sur(:) + real(hhel(:,i), long)/f1
             tilt1_8_sur = tilt1_8_sur + real(tilt1_1_8(i), long)/f1
             tilt8_14_sur = tilt8_14_sur + real(tilt1_8_14(i), long)/f1
             tilt14_20_sur = tilt14_20_sur + real(tilt1_14_20(i), long)/f1
             tilt1_11_sur = tilt1_11_sur + real(tilt1_1_11(i), long)/f1
             tilt11_21_sur = tilt11_21_sur + real(tilt1_11_21(i), long)/f1

             roVal = ro_12(1,i)
             ro_12_sur = ro_12_sur + real(roVal, long)/f1
             radVal = ro_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_sur_sin = ro_12_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_sur_cos = ro_12_sur_cos + real(cosVal,long)/f1

             roVal = ro_15(1,i)
             ro_15_sur = ro_15_sur + real(roVal, long)/f1
             radVal = ro_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_sur_sin = ro_15_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_sur_cos = ro_15_sur_cos + real(cosVal,long)/f1

             roVal = ro_19(1,i)
             ro_19_sur = ro_19_sur + real(roVal, long)/f1
             radVal = ro_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_sur_sin = ro_19_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_sur_cos = ro_19_sur_cos + real(cosVal,long)/f1

             roVal = ro_12(2,i)
             ro_12_sur_avg = ro_12_sur_avg + real(roVal, long)/f1
             roVal = ro_15(2,i)
             ro_15_sur_avg = ro_15_sur_avg + real(roVal, long)/f1
             roVal = ro_19(2,i)
             ro_19_sur_avg = ro_19_sur_avg + real(roVal, long)/f1

             utilt15_20_sur = utilt15_20_sur +real(utilt1_15_20(i), long)/f1
             utilt6_14_sur = utilt6_14_sur + real(utilt1_6_14(i), long)/f1
             if (utiltH1_15_20(i) > 0.0) then
                 utiltH15_20_sur = utiltH15_20_sur + real(utiltH1_15_20(i), long)/f1
                 utiltH15_20 = utiltH15_20 + real(utiltH1_15_20(i), long)/f1
                 pf_ulr_sur = pf_ulr_sur + 1.d0/f1
                 pf_ulr = pf_ulr + 1.d0/f1
             endif

             zscAvg_sur = zscAvg_sur + real(zsc(:,i), long)/f1

             if (isHel_up(i) == 1) then
                 ro_15_sur_hel = ro_15_sur_hel + real(ro_15(1,i), long)/f1
                 ro_19_sur_hel = ro_19_sur_hel + real(ro_19(1,i), long)/f1
                 pf_sur_hel = pf_sur_hel + 1.d0/f1
             endif

!
!  Not first struct
!  Upper leaf
!  inserted
         else
             pf_ins = pf_ins + 1.d0/f1
             cm_pc_ins(:) = cm_pc_ins(:) + real(hcm_pc_up(:,i), long)/f1
             cm_pg_ins(:) = cm_pg_ins(:) + real(hcm_pg_up(:,i), long)/f1
             hel_ins(:) = hel_ins(:) + real(hhel(:,i), long)/f1
             tilt1_8_ins = tilt1_8_ins + real(tilt1_1_8(i), long)/f1
             tilt8_14_ins = tilt8_14_ins + real(tilt1_8_14(i), long)/f1
             tilt14_20_ins = tilt14_20_ins + real(tilt1_14_20(i), long)/f1
             tilt1_11_ins = tilt1_11_ins + real(tilt1_1_11(i), long)/f1
             tilt11_21_ins = tilt11_21_ins + real(tilt1_11_21(i), long)/f1

             ro_12_ins = ro_12_ins + real(ro_12(1,i), long)/f1
             radVal = ro_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_ins_sin = ro_12_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_ins_cos = ro_12_ins_cos + real(cosVal,long)/f1

             ro_15_ins = ro_15_ins + real(ro_15(1,i), long)/f1
             radVal = ro_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_ins_sin = ro_15_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_ins_cos = ro_15_ins_cos + real(cosVal,long)/f1

             ro_19_ins = ro_19_ins + real(ro_19(1,i), long)/f1
             radVal = ro_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_ins_sin = ro_19_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_ins_cos = ro_19_ins_cos + real(cosVal,long)/f1

             ro_12_ins_avg = ro_12_ins_avg + real(ro_12(2,i), long)/f1
             ro_15_ins_avg = ro_15_ins_avg + real(ro_15(2,i), long)/f1
             ro_19_ins_avg = ro_19_ins_avg + real(ro_19(2,i), long)/f1

             utilt15_20_ins = utilt15_20_ins +real(utilt1_15_20(i), long)/f1
             utilt6_14_ins = utilt6_14_ins + real(utilt1_6_14(i), long)/f1
             if (utiltH1_15_20(i) > 0.0) then
                 utiltH15_20_ins = utiltH15_20_ins + real(utiltH1_15_20(i), long)/f1
                 utiltH15_20 = utiltH15_20 + real(utiltH1_15_20(i), long)/f1
                 pf_ulr_ins = pf_ulr_ins + 1.d0/f1
                 pf_ulr = pf_ulr + 1.d0/f1
             endif


             zscAvg_ins = zscAvg_ins + real(zsc(:,i), long)/f1
             if (isHel_up(i) == 1) then
                 ro_15_ins_hel = ro_15_ins_hel + real(ro_15(1,i), long)/f1
                 ro_19_ins_hel = ro_19_ins_hel + real(ro_19(1,i), long)/f1
                 pf_ins_hel = pf_ins_hel + 1.d0/f1
             endif

         endif
         if (-termC_Z2(i) > z21_clust) then
!
!  Not first struct
!  lower leaf
!  surface
             pf_sur = pf_sur + 1.d0/f1
             cm_pc_sur(:) = cm_pc_sur(:) + real(hcm_pc_low(:,i), long)/f1
             cm_pg_sur(:) = cm_pg_sur(:) + real(hcm_pg_low(:,i), long)/f1
             hel_sur(:) = hel_sur(:) + real(hhel2(:,i), long)/f1
             tilt1_8_sur = tilt1_8_sur + real(tilt2_1_8(i), long)/f1
             tilt8_14_sur = tilt8_14_sur + real(tilt2_8_14(i), long)/f1
             tilt14_20_sur = tilt14_20_sur + real(tilt2_14_20(i), long)/f1
             tilt1_11_sur = tilt1_11_sur + real(tilt2_1_11(i), long)/f1
             tilt11_21_sur = tilt11_21_sur + real(tilt2_11_21(i), long)/f1
             roVal = ro2_12(1,i)
             ro_12_sur = ro_12_sur + real(roVal, long)/f1
             radVal = ro2_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_sur_sin = ro_12_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_sur_cos = ro_12_sur_cos + real(cosVal,long)/f1

             roVal = ro2_15(1,i)
             ro_15_sur = ro_15_sur + real(roVal, long)/f1
             radVal = ro2_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_sur_sin = ro_15_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_sur_cos = ro_15_sur_cos + real(cosVal,long)/f1

             roVal = ro2_19(1,i)
             ro_19_sur = ro_19_sur + real(roVal, long)/f1
             radVal = ro2_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_sur_sin = ro_19_sur_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_sur_cos = ro_19_sur_cos + real(cosVal,long)/f1

             roVal = ro2_12(2,i)
             ro_12_sur_avg = ro_12_sur_avg + real(roVal, long)/f1
             roVal = ro2_15(2,i)
             ro_15_sur_avg = ro_15_sur_avg + real(roVal, long)/f1
             roVal = ro2_19(2,i)
             ro_19_sur_avg = ro_19_sur_avg + real(roVal, long)/f1

             utilt15_20_sur = utilt15_20_sur +real(utilt2_15_20(i), long)/f1
             utilt6_14_sur = utilt6_14_sur + real(utilt2_6_14(i), long)/f1
             if (utiltH2_15_20(i) > 0.0) then
                 utiltH15_20_sur = utiltH15_20_sur + real(utiltH2_15_20(i), long)/f1
                 utiltH15_20 = utiltH15_20 + real(utiltH2_15_20(i), long)/f1
                 pf_ulr_sur = 1.d0/f1
                 pf_ulr = 1.d0/f1
             endif
             zscAvg_sur = zscAvg_sur + real(zsc2(:,i), long)/f1
             if (isHel_lo(i) == 1) then
                 ro_15_sur_hel = ro_15_sur_hel + real(ro2_15(1,i), long)/f1
                 ro_19_sur_hel = ro_19_sur_hel + real(ro2_19(1,i), long)/f1
                 pf_sur_hel = pf_sur_hel + 1.d0/f1
             endif

!
!  Not first struct
!  lower leaf
!  inserted
         else
             pf_ins = pf_ins + 1.d0/f1
             cm_pc_ins(:) = cm_pc_ins(:) + real(hcm_pc_low(:,i), long)/f1
             cm_pg_ins(:) = cm_pg_ins(:) + real(hcm_pg_low(:,i), long)/f1
             hel_ins(:) = hel_ins(:) + real(hhel2(:,i), long)/f1
             tilt1_8_ins = tilt1_8_ins + real(tilt2_1_8(i), long)/f1
             tilt8_14_ins = tilt8_14_ins + real(tilt2_8_14(i), long)/f1
             tilt14_20_ins = tilt14_20_ins + real(tilt2_14_20(i), long)/f1
             tilt1_11_ins = tilt1_11_ins + real(tilt2_1_11(i), long)/f1
             tilt11_21_ins = tilt11_21_ins + real(tilt2_11_21(i), long)/f1

             roVal = ro2_12(1,i)
             ro_12_ins = ro_12_ins + real(roVal, long)/f1
             radVal = ro2_12(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_12_ins_sin = ro_12_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_12_ins_cos = ro_12_ins_cos + real(cosVal,long)/f1

             roVal = ro2_15(1,i)
             ro_15_ins = ro_15_ins + real(roVal, long)/f1
             radVal = ro2_15(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_15_ins_sin = ro_15_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_15_ins_cos = ro_15_ins_cos + real(cosVal,long)/f1

             roVal = ro2_19(1,i)
             ro_19_ins = ro_19_ins + real(roVal, long)/f1
             radVal = ro2_19(1,i) / (180.0 / M_PI)
             sinVal = sin(radVal)
             ro_19_ins_sin = ro_19_ins_sin + real(sinVal,long)/f1
             cosVal = cos(radVal)
             ro_19_ins_cos = ro_19_ins_cos + real(cosVal,long)/f1

             roVal = ro2_12(2,i)
             ro_12_ins_avg = ro_12_ins_avg + real(roVal, long)/f1
             roVal = ro2_15(2,i)
             ro_15_ins_avg = ro_15_ins_avg + real(roVal, long)/f1
             roVal = ro2_19(2,i)
             ro_19_ins_avg = ro_19_ins_avg + real(roVal, long)/f1

             utilt15_20_ins = utilt15_20_ins +real(utilt2_15_20(i), long)/f1
             utilt6_14_ins = utilt6_14_ins + real(utilt2_6_14(i), long)/f1
             if (utiltH2_15_20(i) > 0.0) then
                 utiltH15_20_ins = utiltH15_20_ins + real(utiltH2_15_20(i), long)/f1
                 utiltH15_20 = utiltH15_20 + real(utiltH2_15_20(i), long)/f1
                 pf_ulr_ins = 1.d0/f1
                 pf_ulr = 1.d0/f1
             endif
             zscAvg_ins = zscAvg_ins + real(zsc2(:,i), long)/f1
             if (isHel_up(i) == 1) then
                 ro_15_ins_hel = ro_15_ins_hel + real(ro2_15(1,i), long)/f1
                 ro_19_ins_hel = ro_19_ins_hel + real(ro2_19(1,i), long)/f1
                 pf_ins_hel = pf_ins_hel + 1.d0/f1
             endif

         endif
         ro_12_1 = ro_12_1 + real(ro_12(1,i), long)/f1
         radVal = ro_12(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_12_1_sin = ro_12_1_sin + real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_12_1_cos = ro_12_1_cos + real(cosVal,long)/f1

         ro_12_1_avg = ro_12_1_avg + real(ro_12(2,i), long)/f1
         ro_15_1 = ro_15_1 + real(ro_15(1,i), long)/f1
         radVal = ro_15(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_15_1_sin = ro_15_1_sin + real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_15_1_cos = ro_15_1_cos + real(cosVal,long)/f1
         ro_15_1_avg = ro_15_1_avg + real(ro_15(2,i), long)/f1
         ro_19_1 = ro_19_1 + real(ro_19(1,i), long)/f1
         radVal = ro_19(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_19_1_sin = ro_19_1_sin + real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_19_1_cos = ro_19_1_cos + real(cosVal,long)/f1
         ro_19_1_avg = ro_19_1_avg + real(ro_19(2,i), long)/f1

         ro_12_2 = ro_12_2 + real(ro2_12(1,i), long)/f1
         radVal = ro2_12(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_12_2_sin = ro_12_2_sin + real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_12_2_cos = ro_12_2_cos + real(cosVal,long)/f1
         ro_12_2_avg = ro_12_2_avg + real(ro2_12(2,i), long)/f1

         ro_15_2 = ro_15_2 + real(ro2_15(1,i), long)/f1
         radVal = ro2_15(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_15_2_sin = ro_15_2_sin + real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_15_2_cos = ro_15_2_cos + real(cosVal,long)/f1
         ro_15_2_avg = ro_15_2_avg + real(ro2_15(2,i), long)/f1


         ro_19_2 = ro_19_2 + real(ro2_19(1,i), long)/f1
         radVal = ro2_19(1,i) / (180.0 / M_PI)
         sinVal = sin(radVal)
         ro_19_2_sin = ro_19_2_sin + real(sinVal,long)/f1
         cosVal = cos(radVal)
         ro_19_2_cos = ro_19_2_cos + real(cosVal,long)/f1
         ro_19_2_avg = ro_19_2_avg + real(ro2_19(2,i), long)/f1

         zscAvg(:) = zscAvg(:) + real(zsc(:,i), long)/f1
         zscAvg2(:) = zscAvg2(:) + real(zsc2(:,i), long)/f1
         cm_pc1(:) = cm_pc1(:) + real(hcm_pc_up(:,i), long)/f1
         cm_pc2(:) = cm_pc2(:) + real(hcm_pc_low(:,i), long)/f1
         cm_pg1(:) = cm_pg1(:) + real(hcm_pg_up(:,i), long)/f1
         cm_pg2(:) = cm_pg2(:) + real(hcm_pg_low(:,i), long)/f1
         cm_pep1(:) = cm_pep1(:) + real(hcm_pep1(:,i), long)/f1
         cm_pep2(:) = cm_pep2(:) + real(hcm_pep2(:,i), long)/f1
         eppw = eppw + real(heppw(1,i), long)/f1
         st(:) = st(:) + real(hst(:,i), long)/f1
         str(:) = str(:) + real(hstr(:,i), long)/f1
         hel(:) = hel(:) + real(hhel(:,i), long)/f1
         myRan(:) = myRan(:) + real(hran(:,i), long)/f1
         st2(:) = st2(:) + real(hst2(:,i), long)/f1
         str2(:) = str2(:) + real(hstr2(:,i), long)/f1
         hel2(:) = hel2(:) + real(hhel2(:,i), long)/f1
         myRan2(:) = myRan2(:) + real(hran2(:,i), long)/f1
         myPsi(:) = real(psiVal(:,i), long)/f1
         myPhi(:) = real(phiVal(:,i), long)/f1
         sasaAvg = sasaAvg + real(sasa(i), long)/f1

         emmAvg = emmAvg + real(emm(i), long)/f1
         bondAvg = bondAvg + real(bond(i), long)/f1
         angleAvg = angleAvg + real(angle(i), long)/f1
         dihedAvg = dihedAvg + real(dihed(i), long)/f1
         imprpAvg = imprpAvg + real(imprp(i), long)/f1
         electAvg = electAvg + real(elect(i), long)/f1
         vdwAvg = vdwAvg + real(vdw(i), long)/f1
         emmPAvg = emmPAvg + real(emmP(i), long)/f1
         bondPAvg = bondPAvg + real(bondP(i), long)/f1
         anglePAvg = anglePAvg + real(angleP(i), long)/f1
         dihedPAvg = dihedPAvg + real(dihedP(i), long)/f1
         imprpPAvg = imprpPAvg + real(imprpP(i), long)/f1
         electPAvg = electPAvg + real(electP(i), long)/f1
         vdwPAvg = vdwPAvg + real(vdwP(i), long)/f1

         emmAvgIon = emmAvgIon + real(emmIon(i), long)/f1
         bondAvgIon = bondAvgIon + real(bondIon(i), long)/f1
         angleAvgIon = angleAvgIon + real(angleIon(i), long)/f1
         dihedAvgIon = dihedAvgIon + real(dihedIon(i), long)/f1
         imprpAvgIon = imprpAvgIon + real(imprpIon(i), long)/f1
         electAvgIon = electAvgIon + real(electIon(i), long)/f1
         vdwAvgIon = vdwAvgIon + real(vdwIon(i), long)/f1
         emmPAvgIon = emmPAvgIon + real(emmPIon(i), long)/f1
         bondPAvgIon = bondPAvgIon + real(bondPIon(i), long)/f1
         anglePAvgIon = anglePAvgIon + real(anglePIon(i), long)/f1
         dihedPAvgIon = dihedPAvgIon + real(dihedPIon(i), long)/f1
         imprpPAvgIon = imprpPAvgIon + real(imprpPIon(i), long)/f1
         electPAvgIon = electPAvgIon + real(electPIon(i), long)/f1
         vdwPAvgIon = vdwPAvgIon + real(vdwPIon(i), long)/f1

         kineticAvg = kineticAvg + real(kinetic(i), long)/f1
         volAvg = volAvg + (real(en(5,i), long)/f1)
         rgAvg = rgAvg + (real((rg1(i) + rg2(i)), long)/f1)
         rgInterAvg = rgInterAvg + (real((rg1Inter(i) + rg2Inter(i)), long)/f1)
         tilt1_8 = tilt1_8 + real((tilt1_1_8(i) + tilt2_1_8(i)), long)/f1
         tilt8_14 = tilt8_14 + real((tilt1_8_14(i) + tilt2_8_14(i)), long)/f1
         tilt14_20 = tilt14_20 + real((tilt1_14_20(i) + tilt2_14_20(i)), long)/f1
         tilt1_11 = tilt1_11 + real((tilt1_1_11(i) + tilt2_1_11(i)), long)/f1
         tilt11_21 = tilt11_21 + real((tilt1_11_21(i) + tilt2_11_21(i)), long)/f1
         utilt15_20 = utilt15_20 +real((utilt1_15_20(i) + utilt2_15_20(i)), long)/f1
         utilt6_14 = utilt6_14 + real((utilt1_6_14(i) + utilt2_6_14(i)), long)/f1
         do aa=1,11
             if (zsc(aa,i) < zSurface) then
                 avgIns1 = avgIns1 + (1.0/f1)
                 avgInsR1_1 = avgInsR1_1 + (1.0/f1)
             else if (zsc(aa,i) < zUnbound) then
                 avgSur1 = avgSur1 + (1.0/f1)
                 avgSurR1_1 = avgSurR1_1 + (1.0/f1)
             else
                 avgUnb1 = avgUnb1 + (1.0/f1)
                 avgUnbR1_1 = avgUnbR1_1 + (1.0/f1)
             endif

            if (zsc2(aa,i) < zSurface) then
                 avgIns2 = avgIns2 + (1.0/f1)
                 avgInsR1_2 = avgInsR1_2 + (1.0/f1)
             else if (zsc2(aa,i) < zUnbound) then
                 avgSur2 = avgSur2 + (1.0/f1)
                 avgSurR1_2 = avgSurR1_2 + (1.0/f1)
             else
                 avgUnb2 = avgUnb2 + (1.0/f1)
                 avgUnbR1_2 = avgUnbR1_2 + (1.0/f1)
             endif
         enddo

         do aa=12,21
             if (zsc(aa,i) < zSurface) then
                 avgIns1 = avgIns1 + (1.0/f1)
                 avgInsR2_1 = avgInsR2_1 + (1.0/f1)
             else if (zsc(aa,i) < zUnbound) then
                 avgSur1 = avgSur1 + (1.0/f1)
                 avgSurR2_1 = avgSurR2_1 + (1.0/f1)
             else
                 avgUnb1 = avgUnb1 + (1.0/f1)
                 avgUnbR2_1 = avgUnbR2_1 + (1.0/f1)
             endif

            if (zsc2(aa,i) < zSurface) then
                 avgIns2 = avgIns2 + (1.0/f1)
                 avgInsR2_2 = avgInsR2_2 + (1.0/f1)
             else if (zsc2(aa,i) < zUnbound) then
                 avgSur2 = avgSur2 + (1.0/f1)
                 avgSurR2_2 = avgSurR2_2 + (1.0/f1)
             else
                 avgUnb2 = avgUnb2 + (1.0/f1)
                 avgUnbR2_2 = avgUnbR2_2 + (1.0/f1)
             endif
         enddo

     endif
!     print*, 'Rep, f1', currRep, currRepCnt, (1.0/f1)
     totUsed(currRep) = totUsed(currRep) + (1.0/f1)
     currRepCnt = currRepCnt + 1
     if (currRepCnt >= nk(currRep)) then
         currRep = currRep + 1
         currRepCnt = 0
     endif
     helCnt = 0
     betaCnt = 0
     do j = 0, nmon
         betaCnt = betaCnt + hstr(j,i)
     enddo
     do j = 0, nmon
         helCnt = helCnt + hhel(j,i)
     enddo

     helCnt2 = 0
     betaCnt2 = 0
     do j = 0, nmon
         betaCnt2 = betaCnt2 + hstr2(j,i)
     enddo
     do j = 0, nmon
         helCnt2 = helCnt2 + hhel2(j,i)
     enddo

     lrc_hel(lrc_up(i),helCnt) = lrc_hel(lrc_up(i),helCnt) + (1.0/f1)
     lrc_hel(lrc_lo(i),helCnt2) = lrc_hel(lrc_lo(i),helCnt2) + (1.0/f1)

     lenVal = protLen(i) * 2.0
     lenValInt = int(lenVal)
     if (lenValInt >= 120) then
         lenValInt = 119
     endif

!     rgVal = rg(i) * 4.0
!     rgValInt = int(rgVal)
!     if (rgValInt >= 120) then
!         rgValInt = 119
!     endif

     contCntVal = contCnt(i)
     contCntValInt = int(contCntVal)
     if (contCntValInt >=21) then
         contCntValInt = 20
     endif

     if (helCnt < 0.0) then
         print*, "ERROR HELCNT ", helCnt
         helCnt = 0.0
     endif
     if (betaCnt < 0.0) then
         print*, "ERROR BETACNT ", betaCnt
         betaCnt = 0.0
     endif

     do n=1,20
         phi = phiVal(n,i)
         psi = psiVal(n,i)
         phiBuck = int((phi + 180)/5)
         psiBuck = int((psi + 180)/5)
         ramach(psiBuck,phiBuck) = ramach(psiBuck, phiBuck) + (1.0/f1)
         phiHist(n,phiBuck) = phiHist(n,phiBuck) + (1.0/f1)
         psiHist(n,psiBuck) = psiHist(n,psiBuck) + (1.0/f1)
     enddo

     do n=1,21
       zscVal = zsc(n,i)
       zBuck = int(zscVal*2)
       if (zBuck >= 69) then
           zBuck = 69
       endif
       if (zBuck < 0) then
           zBuck = 0
       endif
       zVal_2d(n,zBuck) = zVal_2d(n,zBuck) + (1.0/f1)
       if (zsc(nmon,i) > z21_clust) then
           zVal_2d_sur(n,zBuck) = zVal_2d_sur(n,zBuck) + (1.0/f1)
       else
           zVal_2d_ins(n,zBuck) = zVal_2d_ins(n,zBuck) + (1.0/f1)
       endif
       if (currClust1 .NE. 0) then
           clustZVal_2d(currClust1, n, zBuck) = clustZVal_2d(currClust1, n, zBuck) + (1.0/f1)
       endif


       zscVal = zsc2(n,i)
       zBuck = int(zscVal*2)
       if (zBuck >= 69) then
         zBuck = 69
       endif
       if (zBuck < 0) then
         zBuck = 0
       endif
       zVal_2d(n,zBuck) = zVal_2d(n,zBuck) + (1.0/f1)
       zVal_2d_sum = zVal_2d_sum + (2.0/f1)
       if (zsc2(nmon,i) > z21_clust) then
           zVal_2d_sur(n,zBuck) = zVal_2d_sur(n,zBuck) + (1.0/f1)
       else
           zVal_2d_ins(n,zBuck) = zVal_2d_ins(n,zBuck) + (1.0/f1)
       endif
!
       if (currClust2 .NE. 0) then
           clustZVal_2d(currClust2, n, zBuck) = clustZVal_2d(currClust2, n, zBuck) + (1.0/f1)
       endif
     enddo

     if (currClust1 .NE. 0) then
       clustCmPC(currClust1,:) = clustCmPC(currClust1,:) + (real(hcm_pc_up(:,i), long)/f1)
       clustCmPG(currClust1,:) = clustCmPG(currClust1,:) + (real(hcm_pg_up(:,i), long)/f1)
     endif

     if (currClust2 .NE. 0) then
       clustCmPC(currClust2,:) = clustCmPC(currClust2,:) + (real(hcm_pc_low(:,i), long)/f1)
       clustCmPG(currClust2,:) = clustCmPG(currClust2,:) + (real(hcm_pg_low(:,i), long)/f1)
     endif



     helLen(helCnt, lenValInt) = helLen(helCnt, lenValInt) + (1.0/f1)
     betaLen(betaCnt, lenValInt) = betaLen(betaCnt, lenValInt) + (1.0/f1)
     betaCont(betaCnt, contCntValInt) = betaCont(betaCnt, contCntValInt) + (1.0/f1)
!     betaRg(betaCnt, rgValInt) = betaRg(betaCnt, rgValInt) + (1.0/f1)

     contLen(contCntValInt, lenValInt) = contLen(contCntValInt, lenValInt) + (1.0/f1)
!     contRg(contCntValInt, rgValInt) = contRg(contCntValInt, rgValInt) + (1.0/f1)

!     rgHist(rgValInt) = rgHist(rgValInt) + (1.0/f1)

     lenHist(lenValInt) = lenHist(lenValInt) + (1.0/f1)
     contHist(contCntValInt) = contHist(contCntValInt) + (1.0/f1)

   enddo

   utilt15_20 = utilt15_20/pf
   utilt6_14 = utilt6_14/pf
   utilt15_20_sur = utilt15_20_sur/pf_sur
   utilt6_14_sur = utilt6_14_sur/pf_sur
   utilt15_20_ins = utilt15_20_ins/pf_ins
   utilt6_14_ins = utilt6_14_ins/pf_ins

   tilt1_8 = tilt1_8 /pf
   tilt8_14 = tilt8_14 /pf
   tilt14_20 = tilt14_20 /pf
   tilt1_11 = tilt1_11 /pf
   tilt11_21 = tilt11_21 /pf

   tilt1_8_sur = tilt1_8_sur /pf_sur
   tilt8_14_sur = tilt8_14_sur /pf_sur
   tilt14_20_sur = tilt14_20_sur /pf_sur
   tilt1_11_sur = tilt1_11_sur /pf_sur
   tilt11_21_sur = tilt11_21_sur /pf_sur

   tilt1_8_ins = tilt1_8_ins /pf_ins
   tilt8_14_ins = tilt8_14_ins /pf_ins
   tilt14_20_ins = tilt14_20_ins /pf_ins
   tilt1_11_ins = tilt1_11_ins /pf_ins
   tilt11_21_ins = tilt11_21_ins /pf_ins

   ro_12_1 = ro_12_1/pf
   ro_12_1_vec = atan(ro_12_1_sin/ro_12_1_cos) * (180.0/ M_PI)

   ro_12_1_avg = ro_12_1_avg/pf
   ro_12_2 = ro_12_2/pf
   ro_12_2_vec = atan(ro_12_2_sin/ro_12_2_cos) * (180.0/ M_PI)

   ro_12_2_avg = ro_12_2_avg/pf
   ro_12_ins = ro_12_ins/pf_ins
   ro_12_ins_vec = atan(ro_12_ins_sin/ro_12_ins_cos) * (180.0/ M_PI)
   ro_12_sur_vec = atan(ro_12_sur_sin/ro_12_sur_cos) * (180.0/ M_PI)

   ro_12_ins_avg = ro_12_ins_avg/pf_ins
   ro_12_sur = ro_12_sur/pf_sur
   ro_12_sur_avg = ro_12_sur_avg/pf_sur

   ro_15_1 = ro_15_1/pf
   ro_15_1_avg = ro_15_1_avg/pf
   ro_15_1_vec = atan(ro_15_1_sin/ro_15_1_cos) * (180.0/ M_PI)
   ro_15_2 = ro_15_2/pf
   ro_15_2_avg = ro_15_2_avg/pf
   ro_15_2_vec = atan(ro_15_2_sin/ro_15_2_cos) * (180.0/ M_PI)
   ro_15_ins = ro_15_ins/pf_ins
   ro_15_ins_vec = atan(ro_15_ins_sin/ro_15_ins_cos) * (180.0/ M_PI)
   ro_15_sur_vec = atan(ro_15_sur_sin/ro_15_sur_cos) * (180.0/ M_PI)
   ro_15_ins_avg = ro_15_ins_avg/pf_ins
   ro_15_sur = ro_15_sur/pf_sur
   ro_15_sur_avg = ro_15_sur_avg/pf_sur
   ro_15_sur_hel = ro_15_sur_hel / pf_sur_hel
   ro_15_ins_hel = ro_15_ins_hel / pf_ins_hel

   ro_19_1 = ro_19_1/pf
   ro_19_1_avg = ro_19_1_avg/pf
   ro_19_1_vec = atan(ro_19_1_sin/ro_19_1_cos) * (180.0/ M_PI)
   ro_19_2 = ro_19_2/pf
   ro_19_2_avg = ro_19_2_avg/pf
   ro_19_2_vec = atan(ro_19_2_sin/ro_19_2_cos) * (180.0/ M_PI)
   ro_19_ins = ro_19_ins/pf_ins
   ro_19_ins_vec = atan(ro_19_ins_sin/ro_19_ins_cos) * (180.0/ M_PI)
   ro_19_sur_vec = atan(ro_19_sur_sin/ro_19_sur_cos) * (180.0/ M_PI)
   ro_19_ins_avg = ro_19_ins_avg/pf_ins
   ro_19_sur = ro_19_sur/pf_sur
   ro_19_sur_avg = ro_19_sur_avg/pf_sur
   ro_19_sur_hel = ro_19_sur_hel / pf_sur_hel
   ro_19_ins_hel = ro_19_ins_hel / pf_ins_hel



   clustRo12_vec = atan(clustRo12_sin/ clustRo12_cos) * (180.0/ M_PI)
   clustRo15_vec = atan(clustRo15_sin/ clustRo15_cos) * (180.0/ M_PI)
   clustRo19_vec = atan(clustRo19_sin/ clustRo19_cos) * (180.0/ M_PI)


!   print*, "GET TILT VEC C SIN", clustTiltC_sin
!   print*, "GET TILT VEC C COS", clustTiltC_cos
!   print*, "GET TILT VEC N SIN", clustTiltN_sin
!   print*, "GET TILT VEC N COS", clustTiltN_cos
   clustTiltN_vec = atan(clustTiltN_sin/ clustTiltN_cos) * (180.0/ M_PI)
   clustTiltC_vec = atan(clustTiltC_sin/ clustTiltC_cos) * (180.0/ M_PI)
!   print*, "GET TILT VEC C VEC", clustTiltC_vec
!   print*, "GET TILT VEC N VEC", clustTiltN_vec

   do i=1,4
     clustZVal_2d(i,:,:) = clustZVal_2d(i,:,:)/clustCnt(i)
     clustCmPC(i,:) = clustCmPC(i,:) /clustCnt(i)
     clustCmPG(i,:) = clustCmPG(i,:) /clustCnt(i)
   enddo

   avgIns2 = avgIns2 / pf
   avgSur2 = avgSur2 / pf
   avgUnb2 = avgUnb2 / pf
   avgIns1 = avgIns1 / pf
   avgSur1 = avgSur1 / pf
   avgUnb1 = avgUnb1 / pf

   avgInsR1_2 = avgInsR1_2 / pf
   avgSurR1_2 = avgSurR1_2 / pf
   avgUnbR1_2 = avgUnbR1_2 / pf
   avgInsR1_1 = avgInsR1_1 / pf
   avgSurR1_1 = avgSurR1_1 / pf
   avgUnbR1_1 = avgUnbR1_1 / pf

   avgInsR2_2 = avgInsR2_2 / pf
   avgSurR2_2 = avgSurR2_2 / pf
   avgUnbR2_2 = avgUnbR2_2 / pf
   avgInsR2_1 = avgInsR2_1 / pf
   avgSurR2_1 = avgSurR2_1 / pf
   avgUnbR2_1 = avgUnbR2_1 / pf

   cm_pc1 = cm_pc1 / pf 
   cm_pc2 = cm_pc2 / pf 
   cm_pg1 = cm_pg1 / pf 
   cm_pg2 = cm_pg2 / pf 
   cm_pep1 = cm_pep1 / pf
   cm_pep2 = cm_pep2 / pf
   eppw = eppw / pf
   st = st / pf
   st2 = st2 / pf
   str = str /pf
   str2 = str2 /pf
   hel = hel /pf
   hel2 = hel2 /pf
   hel_sur = hel_sur/pf_sur
   hel_ins = hel_ins/pf_ins
   cm_pc_sur = cm_pc_sur/pf_sur
   cm_pc_ins = cm_pc_ins/pf_ins
   cm_pg_sur = cm_pg_sur/pf_sur
   cm_pg_ins = cm_pg_ins/pf_ins
   myRan = myRan /pf
   myRan2 = myRan2 /pf
   zscAvg = zscAvg /pf
   zscAvg_sur = zscAvg_sur /pf_sur
   zscAvg_ins = zscAvg_ins /pf_ins
   zscAvg2 = zscAvg2 /pf
   zVal_2d = zVal_2d /pf
   zVal_2d_sur = zVal_2d_sur /pf_sur
   zVal_2d_ins = zVal_2d_ins /pf_ins
   z_ro15 = z_ro15/(pf*2)
   z_ro19 = z_ro19/(pf*2)
   z_ro15_hel = z_ro15_hel/isHel_pf
   z_ro19_hel = z_ro19_hel/isHel_pf
   lrc_hel = lrc_hel/(pf*2)
!   helLen = helLen / pf
!   betaCont = betaCont / pf
!   betaLen = betaLen / pf
!   betaRg = betaRg / pf
!   contLen = contLen / pf
!   contRg = contRg / pf
!   rgHist = rgHist / pf
!   lenHist = lenHist / pf
!   contHist = contHist / pf
!   ramach = ramach / (pf*20)
!   psiHist = psiHist / (pf)
!   phiHist = phiHist / (pf)
   pEnergyHist = pEnergyHist / pf
   cTermHist1 = cTermHist1 /pf
   cTermHist2 = cTermHist2 /pf
   nTermHist1 = nTermHist1 /pf
   nTermHist2 = nTermHist2 /pf

   sasaAvg = sasaAvg / pf
   emmAvg = emmAvg / pf
   bondAvg = bondAvg / pf
   angleAvg = angleAvg / pf
   dihedAvg = dihedAvg / pf
   imprpAvg = imprpAvg / pf
   electAvg = electAvg / pf
   vdwAvg = vdwAvg / pf
   volAvg = volAvg / pf
   kineticAvg = kineticAvg / pf

   emmPAvg = emmPAvg / pf
   bondPAvg = bondPAvg /pf
   anglePAvg = anglePAvg /pf
   dihedPAvg = dihedPAvg /pf
   imprpPAvg = imprpPAvg /pf
   electPAvg = electPAvg /pf
   vdwPAvg =  vdwPAvg /pf

   emmAvgIon = emmAvgIon / pf
   bondAvgIon = bondAvgIon / pf
   angleAvgIon = angleAvgIon / pf
   dihedAvgIon = dihedAvgIon / pf
   imprpAvgIon = imprpAvgIon / pf
   electAvgIon = electAvgIon / pf
   vdwAvgIon = vdwAvgIon / pf

   emmPAvgIon = emmPAvgIon / pf
   bondPAvgIon = bondPAvgIon /pf
   anglePAvgIon = anglePAvgIon /pf
   dihedPAvgIon = dihedPAvgIon /pf
   imprpPAvgIon = imprpPAvgIon /pf
   electPAvgIon = electPAvgIon /pf
   vdwPAvgIon =  vdwPAvgIon /pf

   rgAvg = rgAvg / (pf *2)
   rgInterAvg = rgInterAvg / (pf *2)
!   totPAvg = totPAvg /pf



 ramachSum = 0.0
 do i=0,72
   do j=0,72
     ramachSum = ramachSum +  ramach(i,j)
   enddo
 enddo

   lenProb = 0

!   write(43,*)real(nt),real(cm,long)
   write(44,*)real(nt),real(eppw,long)
   write(45,*)real(nt),real(((st+st2)/2),long)
   write(46,*)real(nt),real(((str+str2)/2),long)
   write(47,*)real(nt),real(((hel+hel2)/2),long)
   write(48,*)real(nt),real(((myRan+myRan2)/2),long)
!   write(49,*)real(nt),real(helLen,long)
!   write(50,*)real(nt),real(rgHist,long)
!   write(51,*)real(nt),real(lenHist,long)
!   write(52,*)real(nt),real(betaLen,long)
!   write(53,*)real(nt),real(betaRg,long)
!   write(54,*)real(nt),real(contRg,long)
!   write(55,*)real(nt),real(contLen,long)
!   write(56,*)real(nt),real(contHist,long)
!   write(57,*)real(ramach,long)
!   write(58,*)real(psiHist,long)
!   write(59,*)real(phiHist,long)
!   write(60,*)real(betaCont,long)
   write(61,*)real(pEnergyHist,long)
   write(62,*)real(sasaAvg,long)
   write(63,*)real(emmAvg,long)
   write(64,*)real(bondAvg,long)
   write(65,*)real(angleAvg,long)
   write(66,*)real(dihedAvg,long)
   write(67,*)real(imprpAvg,long)
   write(68,*)real(electAvg,long)
   write(69,*)real(vdwAvg,long)
   write(70,*)real(kineticAvg,long)

   write(71,*)real(bondPAvg,long)
   write(72,*)real(anglePAvg,long)
   write(73,*)real(dihedPAvg,long)
   write(74,*)real(imprpPAvg,long)
   write(75,*)real(electPAvg,long)
   write(76,*)real(vdwPAvg,long)
   write(77,*)real(emmPAvg,long)
   write(78,*)real(volAvg,long)

   write(137,*)"EMM ION   ",   real(emmAvgIon,long),  real(emmPAvgIon,long)
   write(137,*)"BOND ION  ",  real(bondAvgIon,long),  real(bondPAvgIon,long)
   write(137,*)"ANGLE ION ", real(angleAvgIon,long),  real(anglePAvgIon,long)
   write(137,*)"Dihed ION ", real(dihedAvgIon,long),  real(dihedPAvgIon,long)
   write(137,*)"impro ION ", real(imprpAvgIon,long),  real(imprpPAvgIon,long)
   write(137,*)"elect ION ", real(electAvgIon,long),  real(electPAvgIon,long)
   write(137,*)"VDW ION   ",   real(vdwAvgIon,long),  real(vdwPAvgIon,long)

   write(79,*)real(rgAvg,long)
   write(80,*)real(rgInterAvg,long)
   write(81,*)real(cm_pc1,long)
   write(82,*)real(cm_pg1,long)
   write(83,*)real(nt),real(((zscAvg+zscAvg2)/2),long)
   do i=1,21
       write(84,*)i,real((zVal_2d(i,:)),long)
   enddo
   write(85,*)real(cm_pc2,long)
   write(86,*)real(cm_pg2,long)

   write(87,*)real(hel_sur,long)
   write(88,*)real(hel_ins,long)
   write(89,*)real(cm_pc_sur,long)
   write(90,*)real(cm_pc_ins,long)
   write(91,*)real(cm_pg_sur,long)
   write(92,*)real(cm_pg_ins,long)
   write(93,*)real(cm_pep1,long)
   write(94,*)real(cm_pep2,long)
   write(95,*)real(tilt1_8/2,long),real(tilt1_8_ins,long), real(tilt1_8_sur,long)
   write(95,*)real(tilt8_14/2,long), real(tilt8_14_ins,long), real(tilt8_14_sur,long)
   write(95,*)real(tilt14_20/2,long), real(tilt14_20_ins,long),  real(tilt14_20_sur,long)
   write(95,*)real(tilt1_11/2,long), real(tilt1_11_ins,long), real(tilt1_11_sur,long)
   write(95,*)real(tilt11_21/2,long), real(tilt11_21_ins,long), real(tilt11_21_sur,long)
   write(96,*)real(nt),real(zscAvg_ins,long)
   write(97,*)real(nt),real(zscAvg_sur,long)
   do i=1,21
       write(98,*)i,real((zVal_2d_ins(i,:)),long)
   enddo
   do i=1,21
       write(99,*)i,real((zVal_2d_sur(i,:)),long)
   enddo
   write(100,*)real(utilt15_20/2,long),real(utilt15_20_ins,long), real(utilt15_20_sur,long)
   write(100,*)real(utiltH15_20/2,long),real(utiltH15_20_ins,long), real(utiltH15_20_sur,long)
   write(100,*)real(utilt6_14/2,long),real(utilt6_14_ins,long), real(utilt6_14_sur,long)

   write(100,*)real(utilt15_20/2,long),real(utilt15_20_ins,long), real(utilt15_20_sur,long)

   if (ro_12_1_cos < 0.0) then
       if (ro_12_1_sin < 0.0) then
           ro_12_1_vec = ro_12_1_vec + 180
       else
           ro_12_1_vec = ro_12_1_vec + 180
       endif
   else
       if (ro_12_1_sin < 0.0) then
           ro_12_1_vec = ro_12_1_vec + 360
       endif
   endif

   if (ro_12_2_cos < 0.0) then
       if (ro_12_2_sin < 0.0) then
           ro_12_2_vec = ro_12_2_vec + 180
       else
           ro_12_2_vec = ro_12_2_vec + 180
       endif
   else
       if (ro_12_2_sin < 0.0) then
           ro_12_2_vec = ro_12_2_vec + 360
       endif
   endif



   if (ro_15_1_cos < 0.0) then
       if (ro_15_1_sin < 0.0) then
           ro_15_1_vec = ro_15_1_vec + 180
       else
           ro_15_1_vec = ro_15_1_vec + 180
       endif
   else
       if (ro_15_1_sin < 0.0) then
           ro_15_1_vec = ro_15_1_vec + 360
       endif
   endif

   if (ro_15_2_cos < 0.0) then
       if (ro_15_2_sin < 0.0) then
           ro_15_2_vec = ro_15_2_vec + 180
       else
           ro_15_2_vec = ro_15_2_vec + 180
       endif
   else
       if (ro_15_2_sin < 0.0) then
           ro_15_2_vec = ro_15_2_vec + 360
       endif
   endif

   if (ro_19_1_cos < 0.0) then
       if (ro_19_1_sin < 0.0) then
           ro_19_1_vec = ro_19_1_vec + 180
       else
           ro_19_1_vec = ro_19_1_vec + 180
       endif
   else
       if (ro_19_1_sin < 0.0) then
           ro_19_1_vec = ro_19_1_vec + 360
       endif
   endif

   if (ro_19_2_cos < 0.0) then
       if (ro_19_2_sin < 0.0) then
           ro_19_2_vec = ro_19_2_vec + 180
       else
           ro_19_2_vec = ro_19_2_vec + 180
       endif
   else
       if (ro_19_2_sin < 0.0) then
           ro_19_2_vec = ro_19_2_vec + 360
       endif
   endif


   if (ro_12_ins_cos < 0.0) then 
       if (ro_12_ins_sin < 0.0) then
           ro_12_ins_vec = ro_12_ins_vec + 180
       else
           ro_12_ins_vec = ro_12_ins_vec + 180
       endif
   else
       if (ro_12_ins_sin < 0.0) then
           ro_12_ins_vec = ro_12_ins_vec + 360
       endif
   endif

   if (ro_12_sur_cos < 0.0) then 
       if (ro_12_sur_sin < 0.0) then
           ro_12_sur_vec = ro_12_sur_vec + 180
       else
           ro_12_sur_vec = ro_12_sur_vec + 180
       endif
   else
       if (ro_12_sur_sin < 0.0) then
           ro_12_sur_vec = ro_12_sur_vec + 360
       endif
   endif

   if (ro_15_ins_cos < 0.0) then 
       if (ro_15_ins_sin < 0.0) then
           ro_15_ins_vec = ro_15_ins_vec + 180
       else
           ro_15_ins_vec = ro_15_ins_vec + 180
       endif
   else
       if (ro_15_ins_sin < 0.0) then
           ro_15_ins_vec = ro_15_ins_vec + 360
       endif
   endif

   if (ro_15_sur_cos < 0.0) then 
       if (ro_15_sur_sin < 0.0) then
           ro_15_sur_vec = ro_15_sur_vec + 180
       else
           ro_15_sur_vec = ro_15_sur_vec + 180
       endif
   else
       if (ro_15_sur_sin < 0.0) then
           ro_15_sur_vec = ro_15_sur_vec + 360
       endif
   endif


   if (ro_19_ins_cos < 0.0) then 
       if (ro_19_ins_sin < 0.0) then
           ro_19_ins_vec = ro_19_ins_vec + 180
       else
           ro_19_ins_vec = ro_19_ins_vec + 180
       endif
   else
       if (ro_19_ins_sin < 0.0) then
           ro_19_ins_vec = ro_19_ins_vec + 360
       endif
   endif

   if (ro_19_sur_cos < 0.0) then 
       if (ro_19_sur_sin < 0.0) then
           ro_19_sur_vec = ro_19_sur_vec + 180
       else
           ro_19_sur_vec = ro_19_sur_vec + 180
       endif
   else
       if (ro_19_sur_sin < 0.0) then
           ro_19_sur_vec = ro_19_sur_vec + 360
       endif
   endif

   do i=1,numStates
     if (clustRo12_cos(i) < 0.0) then
       clustRo12_vec(i) = clustRo12_vec(i) + 180
     else
       if (clustRo12_sin(i) < 0.0) then
         clustRo12_vec(i) = clustRo12_vec(i) + 360
       endif
     endif

     if (clustRo15_cos(i) < 0.0) then
       clustRo15_vec(i) = clustRo15_vec(i) + 180
     else
       if (clustRo15_sin(i) < 0.0) then
         clustRo15_vec(i) = clustRo15_vec(i) + 360
       endif
     endif

     if (clustRo19_cos(i) < 0.0) then
       clustRo19_vec(i) = clustRo19_vec(i) + 180
     else
       if (clustRo19_sin(i) < 0.0) then
         clustRo19_vec(i) = clustRo19_vec(i) + 360
       endif
     endif

     if (clustTiltN_cos(i) < 0.0) then
       clustTiltN_vec(i) = clustTiltN_vec(i) + 180
     else
       if (clustTiltN_sin(i) < 0.0) then
         clustTiltN_vec(i) = clustTiltN_vec(i) + 360
       endif
     endif

     if (clustTiltC_cos(i) < 0.0) then
       clustTiltC_vec(i) = clustTiltC_vec(i) + 180
     else
       if (clustTiltC_sin(i) < 0.0) then
         clustTiltC_vec(i) = clustTiltC_vec(i) + 360
       endif
     endif
   enddo

   write(101,*)real(ro_12_1,long), real(ro_12_1_avg,long), real(ro_12_1_vec, long), &
     real(ro_12_2,long), real(ro_12_2_avg,long), real(ro_12_2_vec, long), &
     real(ro_12_ins,long), real(ro_12_ins_avg,long), real(ro_12_ins_vec, long), &
     real(ro_12_sur,long), real(ro_12_sur_avg,long), real(ro_12_sur_vec, long)

   write(101,*)real(ro_15_1,long), real(ro_15_1_avg,long), real(ro_15_1_vec, long), &
     real(ro_15_2,long), real(ro_15_2_avg,long), real(ro_15_2_vec, long), &
     real(ro_15_ins,long), real(ro_15_ins_avg,long), real(ro_15_ins_vec, long), &
     real(ro_15_sur,long), real(ro_15_sur_avg,long), real(ro_15_sur_vec, long)

   write(101,*)real(ro_19_1,long), real(ro_19_1_avg,long), real(ro_19_1_vec, long), &
     real(ro_19_2,long), real(ro_19_2_avg,long), real(ro_19_2_vec, long), &
     real(ro_19_ins,long), real(ro_19_ins_avg,long), real(ro_19_ins_vec, long), &
     real(ro_19_sur,long), real(ro_19_sur_avg,long), real(ro_19_sur_vec, long)

!   do i=0,69
!      write(102,*)i-34.5,cTermHist1(i), cTermHist2(i), nTermHist1(i), nTermHist2(i)
!   enddo
   do i=0,34
      off = 34-i
      write(102,*)i+0.5,cTermHist2(off), nTermHist2(off)
   enddo
   do i=0,34
      off = 35+i
      write(103,*)i+0.5,cTermHist1(off), nTermHist1(off)
   enddo
!   enddo   (tempemperatures)

  write(104,*) "INS/SUR ", pf_ins/pf_sur

  write(105,*) avgInsR1_1, avgInsR2_1, avgIns1
  write(105,*) avgSurR1_1, avgSurR2_1, avgSur1
  write(105,*) avgUnbR1_1, avgUnbR2_1, avgUnb1

  write(106,*) avgInsR1_2, avgInsR2_2, avgIns2
  write(106,*) avgSurR1_2, avgSurR2_2, avgSur2
  write(106,*) avgUnbR1_2, avgUnbR2_2, avgUnb2

!   Each row is the same angle (0-6, 6-12 ... 354-360)
!     Each column is same depth (0-0.5, 0.5-1.0 .... 34.5-35.0)
  do i=0,29
    write(107,*)real((z_ro15(i,:)),long)
  enddo
  do i=0,29
    write(108,*)real((z_ro19(i,:)),long)
  enddo

  do i=0,29
    write(110,*)real((z_ro15_hel(i,:)),long)
  enddo
  do i=0,29
    write(111,*)real((z_ro19_hel(i,:)),long)
  enddo

  do i=0,20
    write(109,*)real((lrc_hel(i,:)),long)
  enddo

   write(112,*) real(ro_15_ins_hel,long), real(ro_15_sur_hel,long)
   write(112,*) real(ro_19_ins_hel,long), real(ro_19_sur_hel,long)

  do i=1,numStates
   clustHelC(i) = clustHelC(i)/clustCnt(i)
   write(113,*) real(clustHelC(i),long)
  enddo

  do i=1,numStates
   clustHelN(i) = clustHelN(i)/clustCnt(i)
   write(128,*) real(clustHelN(i),long)
  enddo


  write(114,*) clustRo12_vec
  write(114,*) clustRo15_vec
  write(114,*) clustRo19_vec

  write(115,*) clustTiltN_vec
  write(115,*) clustTiltC_vec

!
!  SB_vals.dat 1 
!  INS_vals.dat 2 
!  MI_vals.dat 3 
!  MSB1_vals.dat 4  
!  MSB2_vals.dat 5
!
!   INS=2
  do i=1,21
      write(116,*)i,real((clustZVal_2d(2,i,:)),long)
  enddo
  write(120,*)real(clustCmPC(2,:),long)
  write(124,*)real(clustCmPG(2,:),long)

! M-INS=3
  do i=1,21
      write(117,*)i,real((clustZVal_2d(3,i,:)),long)
  enddo
  write(121,*)real(clustCmPC(3,:),long)
  write(125,*)real(clustCmPG(3,:),long)

! SB=1
  do i=1,21
      write(118,*)i,real((clustZVal_2d(1,i,:)),long)
  enddo
  write(122,*)real(clustCmPC(1,:),long)
  write(126,*)real(clustCmPG(1,:),long)

! M-SUR1=4
  do i=1,21
      write(119,*)i,real((clustZVal_2d(4,i,:)),long)
  enddo
  write(123,*)real(clustCmPC(4,:),long)
  write(127,*)real(clustCmPG(4,:),long)

! M-SER2=5
  do i=1,21
      write(129,*)i,real((clustZVal_2d(5,i,:)),long)
  enddo
  write(130,*)real(clustCmPC(5,:),long)
  write(131,*)real(clustCmPG(5,:),long)

  do i=1,3
      scdPgRegAvg(i) = scdPgRegAvg(i) / scdPgRegCntAvg(i)
      scdPcRegAvg(i) = scdPcRegAvg(i) / scdPcRegCntAvg(i)
  enddo
  do i=1,3
      write(132,*)scdPcRegAvg(i), scdPgRegAvg(i)
  enddo

  do i=1,3
      pcRegLipDenArea(i) = pcRegLipDenArea(i) / pf
      pgRegLipDenArea(i) = pgRegLipDenArea(i) / pf
      pcRegLipDenVol(i) = pcRegLipDenVol(i) / pf
      pgRegLipDenVol(i) = pgRegLipDenVol(i) / pf
  enddo
  do i=1,3
      write(133,*)PcRegLipDenArea(i), PgRegLipDenArea(i), PcRegLipDenVol(i), PgRegLipDenVol(i)
  enddo

  do i=0,34
    nTerm_cTerm(i,:) =  nTerm_cTerm(i,:) / (pf * 2)
    write(134,*)real((nTerm_cTerm(i,:)),long)
  enddo

  do i=0,34
    mid_cTerm(i,:) =  mid_cTerm(i,:) / (pf * 2)
    write(135,*)real((mid_cTerm(i,:)),long)
  enddo

  do i=1,4
      lipid_coord(i) = lipid_coord(i) / pf
      write(136,*)real((lipid_coord(i)),long)
  enddo

  close(43)
  close(44)
  close(45)
  close(46)
  close(47)
  close(48)
  close(49)
  close(50)
  close(51)
  close(52)
  close(53)
  close(54)
  close(55)
  close(56)
  close(57)
  close(58)
  close(59)
  close(60)
  close(61)
  close(62)
  close(63)
  close(64)
  close(65)
  close(66)
  close(67)
  close(68)
  close(69)
  close(70)
  close(71)
  close(72)
  close(73)
  close(75)
  close(75)
  close(76)
  close(77)
  close(78)
  close(79)
  close(80)
  close(81)
  close(82)
  close(83)
  close(84)
  close(85)
  close(85)
  close(86)
  close(87)
  close(88)
  close(89)
  close(90)
  close(91)
  close(92)
  close(93)
  close(94)
  close(95)
  close(96)
  close(97)
  close(98)
  close(99)
  close(100)
  close(101)
  close(102)
  close(103)
  close(104)
  close(105)
  close(106)
  close(107)
  close(108)
  close(109)
  close(110)
  close(111)
  close(112)
  close(113)
  close(114)
  close(115)
  close(116)
  close(117)
  close(118)
  close(119)
  close(120)
  close(121)
  close(122)
  close(123)
  close(124)
  close(125)
  close(126)
  close(127)
  close(128)
  close(129)
  close(130)
  close(131)
  close(132)
  close(133)
  close(134)
  close(135)
  close(136)

  usedTotal = sum(totUsed)


  open(49,file='th_rep_used_'//trim(trString)//'.dat')

  do currRep = repmin,repmax
      write(49,*)currRep,real(totUsed(currRep)/usedTotal, long)
  enddo
  close(49)

  open(49,file='th_clust_prob'//trim(trString)//'.dat')

  do currClust = 0, numStates
      write(49,*)currClust, clustCntSaddle(currClust)/(pf*2), clustCnt(currClust)/(pf*2)
  enddo
  close(49)




end program hist7


