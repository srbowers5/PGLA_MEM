!    TO", hcm_pc(1:ncont,istr)
!   Input
!      trmin, trmax
!      nrep, repmin, repmax
!      tsim = number of structs
!      teq = number needed to equilibrate
!      e1, e2, edist not used
!
!   Output
!      th_cm.dat - average contact map for each trajectory
!      th_eppw.dat - average energy (from p-->p and p-->w contacts, no w-->w)
!      th_st.dat - average turns
!      histstr.dat strand
!      histst.dat turn
!      histran.dat random coil
!      histhel.dat helix
!
!   Params to modify
!      trmin, trmax
!      repmin, repmax
!      nmon - number of AA
!      nmol - number of peptide molecules
!      interval - 
!      nbin - 
!      nbinz
!      na - number of atoms
!      ncont - contacts = ((nmon * (nmon-1))/2)
!      tsim - steps in simulation
!      teq - steps to equilibration?
program hist1z

implicit none
integer, parameter :: in_tr=1
integer, parameter :: long = selected_real_kind(20,1400)
!integer, parameter :: long = selected_real_kind(15,307)
integer, parameter :: ilong = selected_int_kind(10)
integer, parameter :: trmin=1
integer, parameter :: trmax=1
integer, parameter :: repmin=1,repmax=20, interval=1,nbin=80,&
                      nbinz=80,nmon=21,nmol=2,na=308,ncont_leaf=105,ncont=210, ncont_pep=210, trFirst=1, trLast=1
!integer, dimension(trmin:trmax), parameter :: tsim=(/9999,9999,9999,9999,9999,9999,9999,9999,9999,9999/),&
!                                              teq=(/0,0,0,0,0,0,0,0,0,0/)
integer, dimension(trmin:trmax), parameter :: tsim=(/50000/),&
                                              teq=(/0/)

! Not used  p.
real (long), parameter :: eps=1.d-7, e1 = -4000.d0, e2 = 10000.d0, edist=2.d0, &
                          edelta=0.d0, &
                          p=1.458397d-05

integer :: ic,it,icount,tr,rep,k,i,j,l, m,n,tprod, istr, cnt
integer :: ehp, ehpw, ehw, ehtots
real (long) ::  emin,emax,ep,epw,ew,etot,etots,eminb,emaxb,emins,emaxs,emean, SasaVal

real (long) ::  bondPVal, anglePVal, dihedPVal, imprpPVal, electPVal, vdwPVal, EmmPVal
real (long) :: bondVal, angleVal, dihedVal, imprpVal, electVal, &
        vdwVal, kineticVal, EmmVal, volVal
real (long) ::  bondPValIon, anglePValIon, dihedPValIon, imprpPValIon, electPValIon, vdwPValIon, EmmPValIon
real (long) :: bondValIon, angleValIon, dihedValIon, imprpValIon, electValIon, vdwValIon, EmmValIon
real (long) :: rg1IntVal, rg2IntVal, rg1Val, rg2Val 
real (long), dimension(3) :: box
integer, dimension(:,:), allocatable :: hss,heppw,hst,hstr,hhel,hran
integer, dimension(:,:), allocatable :: hss2,heppw2,hst2,hstr2,hhel2,hran2
integer, dimension(:,:), allocatable :: rg, rgIng
integer, dimension(:,:), allocatable :: hcm_pep1, hcm_pep2
integer, dimension(:,:), allocatable :: hcm_pc_up, hcm_pc_low, hcm_pg_up, hcm_pg_low
integer, dimension(:), allocatable :: totVal, totVal2
real (long), dimension(:), allocatable :: protLen,short_cont,long_cont
real (long), dimension(:), allocatable :: sasa, emm, &
    bond, angle, dihed, imprp, elect, vdw, kinetic, vol, &
    bondP, angleP, dihedP, imprpP, electP, vdwP, emmP
real (long), dimension(:), allocatable :: emmIon, bondIon, angleIon, dihedIon, imprpIon, electIon, vdwIon, &
    bondPIon, anglePIon, dihedPIon, imprpPIon, electPIon, vdwPIon, emmPIon

real (long), dimension(:), allocatable :: tilt1_1_8, tilt1_8_14, tilt1_14_20, tilt1_1_11, tilt1_11_21
real (long), dimension(:), allocatable :: tilt2_1_8, tilt2_8_14, tilt2_14_20, tilt2_1_11, tilt2_11_21
real (long), dimension(:), allocatable :: utilt1_15_20, utiltH1_15_20, utilt1_6_14
real (long), dimension(:), allocatable :: utilt2_15_20, utiltH2_15_20, utilt2_6_14
real (long), dimension(:), allocatable :: tilt1_15_20, tilt1_6_14, tilt2_15_20, tilt2_6_14

real (long), dimension(:), allocatable :: rg1Int, rg2Int, rg1, rg2
real (long), dimension(:), allocatable :: termC_Z1, termN_Z1, termC_Z2, termN_Z2, termM_Z1, termM_Z2, termAll_Z1, termAll_Z2

integer, dimension(:), allocatable :: LRC1, LRC2
integer, dimension(:), allocatable :: isHel1, isHel2


real (long), dimension(:,:), allocatable :: regionDen
real (long), dimension(:,:), allocatable :: lip_coord
real (long), dimension(:,:), allocatable :: en,psiVal1,psiVal2,phiVal1,phiVal2,zsc,zsc2,scdPg,scdPc,scdPgReg,scdpCreg
integer, dimension(:,:), allocatable :: scdPgRegCnt,scdpCregCnt
real (long), dimension(:,:), allocatable :: ro_12, ro2_12, ro_15, ro2_15, ro_19, ro2_19
real (long), dimension(311) :: dummy
integer, dimension(210) :: dummycm
integer, dimension(420) :: dummycm_pep
character(len=nmon) :: dchar
character(len=80), dimension(1:4) :: path,locpath
character(len=16):: stepString
character(len=16):: repString
character(len=80), dimension(1:4) :: intermpath, datapath, utiltpath, scdpath, xscpath
character(len=80), dimension(1:4) :: fedir, sasapath, eneAApath, enepath, eneionpath, eneAAionpath
character(len=80) :: myFileName
real                  :: myLen

type trajectory
  integer               :: step
  real                  :: enp
  real                  :: enpw
  real                  :: enw
  real                  :: entot
  real                  :: ent,entb,ents
  real                  :: v
  real                  :: eppw
  real                  :: protLen
  real                  :: rg1
  real                  :: rg2
  real                  :: rg1Int
  real                  :: rg2Int
  real                  :: ion_cont
  real                  :: short_cont
  real                  :: long_cont
  real                  :: sasa
  real                  :: emm
  real                  :: bond
  real                  :: angle
  real                  :: dihed
  real                  :: imprp
  real                  :: elect
  real                  :: vdw
  real                  :: kinetic
  real                  :: bondP
  real                  :: angleP
  real                  :: dihedP
  real                  :: imprpP
  real                  :: electP
  real                  :: emmP
  real                  :: vdwP
  real                  :: vol

  real                  :: emmIon
  real                  :: bondIon
  real                  :: angleIon
  real                  :: dihedIon
  real                  :: imprpIon
  real                  :: electIon
  real                  :: vdwIon
  real                  :: bondPIon
  real                  :: anglePIon
  real                  :: dihedPIon
  real                  :: imprpPIon
  real                  :: electPIon
  real                  :: emmPIon
  real                  :: vdwPIon

  real, dimension(nmon-2) :: psiVal1
  real, dimension(nmon-2) :: psiVal2
  real, dimension(nmon-2) :: phiVal1
  real, dimension(nmon-2) :: phiVal2
  real, dimension(nmon) :: zsc
  real, dimension(nmon) :: zsc2
  real, dimension(2) :: ro_12
  real, dimension(2) :: ro2_12
  real, dimension(2) :: ro_15
  real, dimension(2) :: ro2_15
  real, dimension(2) :: ro_19
  real, dimension(2) :: ro2_19
  real, dimension(2) :: termC_Z
  real, dimension(2) :: termN_Z
  real, dimension(2) :: termM_Z
  real, dimension(2) :: termAll_Z
  real :: tilt1_1_8
  real :: tilt1_8_14
  real :: tilt1_14_20
  real :: tilt1_1_11
  real :: tilt1_11_21
  real :: tilt2_1_8
  real :: tilt2_8_14
  real :: tilt2_14_20
  real :: tilt2_1_11
  real :: tilt2_11_21
  real :: utilt1_15_20
  real :: utiltH1_15_20
  real :: utilt1_6_14
  real :: utilt2_15_20
  real :: utiltH2_15_20
  real :: utilt2_6_14
  real :: tilt1_6_14
  real :: tilt2_6_14
  real :: tilt1_15_20
  real :: tilt2_15_20
  real, dimension(13) :: scdPg
  real, dimension(13) :: scdPc
  real, dimension(3) :: scdPgReg
  real, dimension(3) :: scdPcReg
  real, dimension(12) :: regionDen
  real, dimension(4) :: lip_coord
  integer, dimension(3) :: scdPgRegCnt
  integer, dimension(3) :: scdPcRegCnt
  integer :: LRC1
  integer :: LRC2
  integer :: isHel1
  integer :: isHel2
  integer, dimension(ncont_leaf) :: cm_pg_up, cm_pg_low, cm_pc_up, cm_pc_low
  integer, dimension(ncont_pep) :: cm_pep1
  integer, dimension(ncont_pep) :: cm_pep2
  integer, dimension(nmon) :: ss,st,sstr,shel,sran
  integer, dimension(nmon) :: ss2,st2,sstr2,shel2,sran2
  character(len=1), dimension(1:nmon) :: sschar1, sschar2
end type trajectory
type (trajectory), dimension(:), allocatable  :: trj

integer :: vals, tens


 path(1) = './'
 path(2) = './'
 path(3) = './'
 path(4) = './'

 locpath(1) = './'
 locpath(2) = './'
 locpath(3) = './'
 locpath(4) = './'

 fedir(1) = '/home/dad/MM-GBSA_HOT_1/'
 fedir(2) = '/home/dad/MM-GBSA_HOT_2/'
 fedir(3) = '/home/dad/MM-GBSA_HOT_3/'

 datapath(1) = '/SCHOOL/DO_REST_PCPG_HOT_H2O_1_CONT/'
 datapath(2) = '/SCHOOL/DO_REST_PCPG_HOT_H2O_2/'
 datapath(3) = '/SCHOOL/DO_REST_PCPG_HOT_H2O_CONT3_3/'

 utiltpath(1) = '/home/dad/GET_TILT_ULRICH_HOT_1/DATA/'
 utiltpath(2) = '/home/dad/GET_TILT_ULRICH_HOT_2/DATA/'
 utiltpath(3) = '/home/dad/GET_TILT_ULRICH_HOT_3/DATA/'

 enepath(1) = trim(fedir(1))//'ENERGY_DIR/'
 enepath(2) = trim(fedir(2))//'ENERGY_DIR/'
 enepath(3) = trim(fedir(3))//'ENERGY_DIR/'

 eneAApath(1) = trim(fedir(1))//'ENERGY_DIR_AA/'
 eneAApath(2) = trim(fedir(2))//'ENERGY_DIR_AA/'
 eneAApath(3) = trim(fedir(3))//'ENERGY_DIR_AA/'

 eneionpath(1) = trim(fedir(1))//'ENERGY_DIR_IONS/'
 eneionpath(2) = trim(fedir(2))//'ENERGY_DIR_IONS/'
 eneionpath(3) = trim(fedir(3))//'ENERGY_DIR_IONS/'

 eneAAionpath(1) = trim(fedir(1))//'ENERGY_DIR_AA_IONS/'
 eneAAionpath(2) = trim(fedir(2))//'ENERGY_DIR_AA_IONS/'
 eneAAionpath(3) = trim(fedir(3))//'ENERGY_DIR_AA_IONS/'




 sasapath(1) = trim(fedir(1))//'SASA_DIR/'
 sasapath(2) = trim(fedir(2))//'SASA_DIR/'
 sasapath(3) = trim(fedir(3))//'SASA_DIR/'

 scdpath(1) = trim(datapath(1))//'SCD_DIR/'
 scdpath(2) = trim(datapath(2))//'SCD_DIR/'
 scdpath(3) = trim(datapath(3))//'SCD_DIR/'

 xscpath(1) = trim(datapath(1))
 xscpath(2) = trim(datapath(2))
 xscpath(3) = trim(datapath(3))

 intermpath(1) = trim(datapath(1))//'interm_out/'
 intermpath(2) = trim(datapath(2))//'interm_out/'
 intermpath(3) = trim(datapath(3))//'interm_out/'


 emin = 100000.d0 ; emax = -100000.d0
 emin = 100000.d0 ; emax = -100000.d0
 eminb = 100000.d0 ; emaxb = -100000.d0
 emins = 100000.d0 ; emaxs = -100000.d0
 icount=0; istr=0
 allocate(hcm_pg_up(ncont_leaf,sum(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pg_low(ncont_leaf,sum(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pc_up(ncont_leaf,sum(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pc_low(ncont_leaf,sum(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pep1(ncont_pep,sum(tsim-teq)*(repmax-repmin+1)),   &
          hcm_pep2(ncont_pep,sum(tsim-teq)*(repmax-repmin+1)),   &
          hss(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hst(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hstr(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hhel(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hran(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hss2(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hst2(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hstr2(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hhel2(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          hran2(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          heppw(1,sum(tsim-teq)*(repmax-repmin+1)),     &
          sasa(sum(tsim-teq)*(repmax-repmin+1)),     &
          emm(sum(tsim-teq)*(repmax-repmin+1)),     &
          emmP(sum(tsim-teq)*(repmax-repmin+1)),     &
          bond(sum(tsim-teq)*(repmax-repmin+1)),     &
          angle(sum(tsim-teq)*(repmax-repmin+1)),     &
          dihed(sum(tsim-teq)*(repmax-repmin+1)),     &
          imprp(sum(tsim-teq)*(repmax-repmin+1)),     &
          elect(sum(tsim-teq)*(repmax-repmin+1)),     &
          vdw(sum(tsim-teq)*(repmax-repmin+1)),     &
          bondP(sum(tsim-teq)*(repmax-repmin+1)),     &
          angleP(sum(tsim-teq)*(repmax-repmin+1)),     &
          dihedP(sum(tsim-teq)*(repmax-repmin+1)),     &
          imprpP(sum(tsim-teq)*(repmax-repmin+1)),     &
          electP(sum(tsim-teq)*(repmax-repmin+1)),     &
          vdwP(sum(tsim-teq)*(repmax-repmin+1)),     &
          kinetic(sum(tsim-teq)*(repmax-repmin+1)),     &
          emmIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          emmPIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          bondIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          angleIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          dihedIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          imprpIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          electIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          vdwIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          bondPIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          anglePIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          dihedPIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          imprpPIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          electPIon(sum(tsim-teq)*(repmax-repmin+1)),     &
          vdwPIon(sum(tsim-teq)*(repmax-repmin+1)),     &

          psiVal1(nmon-2,sum(tsim-teq)*(repmax-repmin+1)), &
          psiVal2(nmon-2,sum(tsim-teq)*(repmax-repmin+1)), &
          phiVal1(nmon-2,sum(tsim-teq)*(repmax-repmin+1)), &
          phiVal2(nmon-2,sum(tsim-teq)*(repmax-repmin+1)), &
          utilt1_15_20(sum(tsim-teq)*(repmax-repmin+1)), &
          utiltH1_15_20(sum(tsim-teq)*(repmax-repmin+1)), &
          utilt1_6_14(sum(tsim-teq)*(repmax-repmin+1)), &
          utilt2_15_20(sum(tsim-teq)*(repmax-repmin+1)), &
          utiltH2_15_20(sum(tsim-teq)*(repmax-repmin+1)), &
          utilt2_6_14(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_1_8(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_8_14(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_14_20(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_1_11(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_11_21(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_1_8(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_8_14(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_14_20(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_1_11(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_11_21(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_6_14(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_6_14(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt1_15_20(sum(tsim-teq)*(repmax-repmin+1)), &
          tilt2_15_20(sum(tsim-teq)*(repmax-repmin+1)), &
          vol(sum(tsim-teq)*(repmax-repmin+1)),     &
          en(5,sum(tsim-teq)*(repmax-repmin+1)), &
          rg1Int(sum(tsim-teq)*(repmax-repmin+1)),     &
          rg2Int(sum(tsim-teq)*(repmax-repmin+1)),     &
          rg1(sum(tsim-teq)*(repmax-repmin+1)),     &
          totVal(nmon), &
          totVal2(nmon), &
          ro_12(2,sum(tsim-teq)*(repmax-repmin+1)), &
          ro2_12(2,sum(tsim-teq)*(repmax-repmin+1)), &
          ro_15(2,sum(tsim-teq)*(repmax-repmin+1)), &
          ro2_15(2,sum(tsim-teq)*(repmax-repmin+1)), &
          ro_19(2,sum(tsim-teq)*(repmax-repmin+1)), &
          ro2_19(2,sum(tsim-teq)*(repmax-repmin+1)), &
          scdPg(13,sum(tsim-teq)*(repmax-repmin+1)), &
          scdPc(13,sum(tsim-teq)*(repmax-repmin+1)), &
          scdPgReg(3,sum(tsim-teq)*(repmax-repmin+1)), &
          scdPcReg(3,sum(tsim-teq)*(repmax-repmin+1)), &
          scdPgRegCnt(3,sum(tsim-teq)*(repmax-repmin+1)), &
          scdPcRegCnt(3,sum(tsim-teq)*(repmax-repmin+1)), &
          regionDen(12,sum(tsim-teq)*(repmax-repmin+1)), &
          lip_coord(4,sum(tsim-teq)*(repmax-repmin+1)), &
          termC_Z1(sum(tsim-teq)*(repmax-repmin+1)), &
          termN_Z1(sum(tsim-teq)*(repmax-repmin+1)), &
          termC_Z2(sum(tsim-teq)*(repmax-repmin+1)), &
          termN_Z2(sum(tsim-teq)*(repmax-repmin+1)), &
          termM_Z1(sum(tsim-teq)*(repmax-repmin+1)), &
          termM_Z2(sum(tsim-teq)*(repmax-repmin+1)), &
          termAll_Z1(sum(tsim-teq)*(repmax-repmin+1)), &
          termAll_Z2(sum(tsim-teq)*(repmax-repmin+1)), &
          zsc(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          zsc2(nmon,sum(tsim-teq)*(repmax-repmin+1)),    &
          LRC1(sum(tsim-teq)*(repmax-repmin+1)),    &
          LRC2(sum(tsim-teq)*(repmax-repmin+1)),    &
          isHel1(sum(tsim-teq)*(repmax-repmin+1)),    &
          isHel2(sum(tsim-teq)*(repmax-repmin+1)),    &
          rg2(sum(tsim-teq)*(repmax-repmin+1)))

! read*,in_tr
 print*, "IN TR", in_tr
 do tr=in_tr,in_tr
   if(tsim(tr)==0)cycle
   tprod = tsim(tr)-teq(tr)
   allocate(trj(tprod))
   do rep=repmin,repmax
     if (rep < 10) then
       repString = char(48+rep)
     else
       tens = rep/10
       vals = rep-(tens*10)
       repString = char(48+tens)//char(48+vals)
     endif


    print*,'tr, rep =',tr,rep,trim(path(tr))
    open(21,file=trim(locpath(tr))//'ENERGY/Eprest2namd'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(22,file=trim(locpath(tr))//'ENERGY/Epwrest2namd'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(23,file=trim(locpath(tr))//'ENERGY/Ewrest2namd'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(24,file=trim(locpath(tr))//'ENERGY/Etotrest2namd'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(27,file=trim(xscpath(tr))//'cryst_$RANGE_'//char(48+tr)//'_'//trim(repString)//'.cry')
    open(28,file=trim(datapath(tr))//'SSU'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
    open(29,file=trim(datapath(tr))//'SSL'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
     open(35,file=trim(enepath(tr))//'energy'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
!     open(36,file=trim(sasapath(tr))//'sasa'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
     open(37,file=trim(eneAApath(tr))//'energyAA'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
    open(38,file=trim(intermpath(tr))//'r'//char(48+tr)//'_'//trim(repString)//'b_$RANGE.dat')
!    open(39,file=trim(datapath(tr))//'rg_'//char(48+tr)//'_'//trim(repString)//'_b_all.dat')
    open(40,file=trim(intermpath(tr))//'contlist'//char(48+tr)//'_'//trim(repString)//'_lpPCb_$RANGE.dat')
    open(41,file=trim(intermpath(tr))//'contlist'//char(48+tr)//'_'//trim(repString)//'_lpPGb_$RANGE.dat')
    open(42,file=trim(intermpath(tr))//'zsc'//char(48+tr)//'_'//trim(repString)//'b_$RANGE.dat')
    open(43,file=trim(intermpath(tr))//'cm'//char(48+tr)//'_'//trim(repString)//'b_$RANGE.dat')
    open(46,file=trim(intermpath(tr))//'dhWT'//char(48+tr)//'_'//trim(repString)//'b_$RANGE.dat')
    open(47,file=trim(utiltpath(tr))//'tilt'//char(48+tr)//'_$RANGE_15-20_'//trim(repString)//'.dat')
    open(49,file=trim(utiltpath(tr))//'tilt'//char(48+tr)//'_$RANGE_6-14_'//trim(repString)//'.dat')
    open(50,file=trim(utiltpath(tr))//'roUp'//char(48+tr)//'_$RANGE_6-14_'//trim(repString)//'.dat')
    open(51,file=trim(utiltpath(tr))//'roUp'//char(48+tr)//'_$RANGE_15-20_'//trim(repString)//'.dat')
    open(52,file=trim(utiltpath(tr))//'roLo'//char(48+tr)//'_$RANGE_6-14_'//trim(repString)//'.dat')
    open(53,file=trim(utiltpath(tr))//'roLo'//char(48+tr)//'_$RANGE_15-20_'//trim(repString)//'.dat')
    open(54,file=trim(intermpath(tr))//'term_com'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(55,file=trim(intermpath(tr))//'LRC_up'//char(48+tr)//'_'//trim(repString)//'.dat')
    open(56,file=trim(intermpath(tr))//'LRC_lo'//char(48+tr)//'_'//trim(repString)//'.dat')
    open(57,file=trim(utiltpath(tr))//'isHelUp'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
    open(58,file=trim(utiltpath(tr))//'isHelLo'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
    open(59,file=trim(scdpath(tr))//'pg_scd_step_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(60,file=trim(scdpath(tr))//'pc_scd_step_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(61,file=trim(scdpath(tr))//'pg_scd_reg_step_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(62,file=trim(scdpath(tr))//'pc_scd_reg_step_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(63,file=trim(scdpath(tr))//'pg_scd_reg_step_cnt_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(64,file=trim(scdpath(tr))//'pc_scd_reg_step_cnt_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(65,file=trim(datapath(tr))//'reg_den_step_'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(66,file=trim(intermpath(tr))//'lip_coord_num'//char(48+tr)//'_'//trim(repString)//'_$RANGE.dat')
    open(67,file=trim(eneionpath(tr))//'energy'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')
    open(68,file=trim(eneAAionpath(tr))//'energyAA'//char(48+tr)//'_$RANGE_'//trim(repString)//'.dat')

    print*, 'For traj ', tr, teq(tr), tsim(tr)
    print*, 'Done open '
!...read in file headers...
!    read(21,*)dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar
!    print*, 'Read Epr '
!    read(24,*)dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar
!    print*, 'Read Etot '
!    read(21,*)dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar
!    read(24,*)dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar

!    read(35,*)dchar, dummy(1:20)
!    read(36,*)dummy(1:2)


    do k=1,teq(tr)
      read(21,*)dummy(1:11)
      read(22,*)dchar,dummy(2:6) 
      read(23,*)dummy(1:11)
      read(24,*)dummy(1:11)
      read(27,*)dchar,dchar,dchar,dchar,dchar,dchar,dchar,dchar, dchar
      read(28,'(a21)')dchar
      read(29,'(a21)')dchar
      read(35,*)dchar, dummy(1:20)
!      read(36,*)dummy(1:2)
      read(37,*)dchar, dummy(1:20)
      read(38,*)dummy(1:5)
!      read(39,*)dummy(1:3)
      read(40,*)it,dummycm
      read(41,*)it,dummycm
      read(42,*)it,dummy(1:42)
      read(43,*)it,dummycm_pep
!      read(44,*)dummy(1:3)
!      read(45,*)dummy(1:3)
      read(46,*)it,dummy(1:38)
      read(50,*)dummy(1:10)
      read(51,*)dummy(1:8)
      read(52,*)dummy(1:10)
      read(53,*)dummy(1:8)
      read(54,*)dummy(1:8)
      read(55,*)it
      read(56,*)it
      read(59,*)dummy(1:26)
      read(60,*)dummy(1:26)
      read(61,*)dummy(1:6)
      read(62,*)dummy(1:6)
      read(63,*)dummy(1:6)
      read(64,*)dummy(1:6)
      read(65,*)dummy(1:12)
      read(66,*)dummy(1:4)
      read(67,*)dchar, dummy(1:20)
      read(68,*)dchar, dummy(1:20)
    enddo   

    do k=teq(tr)+1,tsim(tr)
      icount = icount + 1
      i=k-teq(tr)
      read(21,*)dummy(1:10),ep
      trj(i)%enp = ep
      read(22,*)dchar,dummy(2:5),epw
      trj(i)%enpw = epw
      read(23,*)dummy(1:10),ew
      trj(i)%enw = ew
      read(24,*)dummy(1:10),etot
      trj(i)%entot = etot
      trj(i)%eppw = ep + epw
      read(27,*)dchar,box(1),box(2), box(3), dchar, dchar,dchar,dchar,dchar
      trj(i)%v = box(1)*box(2)*box(3)
      trj(i)%ent = ep + epw + ew + p*trj(i)%v
      read(28,'(21a1)')trj(i)%sschar1
      read(29,'(21a1)')trj(i)%sschar2

      read(35,*)dchar, dummy(2:2), bondVal, angleVal, dihedVal, imprpVal, electVal, &
        vdwVal, dummy(9:10), kineticVal, dummy(12:13), EmmVal, dummy(15:18), volVal, dummy(20:21)
      read(37,*)dchar, dummy(2:2), bondPVal, anglePVal, dihedPVal, imprpPVal, electPVal, &
        vdwPVal, dummy(9:13), EmmPVal, dummy(15:21)

      read(67,*)dchar, dummy(2:2), bondValIon, angleValIon, dihedValIon, imprpValIon, electValIon, &
        vdwValIon, dummy(9:10), dummy(11:11), dummy(12:13), EmmValIon, dummy(15:18), dummy(19:19), dummy(20:21)
      read(68,*)dchar, dummy(2:2), bondPValIon, anglePValIon, dihedPValIon, imprpPValIon, electPValIon, &
        vdwPValIon, dummy(9:13), EmmPValIon, dummy(15:21)

      read(38,*)dchar, rg1IntVal, rg2IntVal, dummy(4:5)
!      read(39,*)dchar, rg1Val, rg2Val
!      read(40,'(i6,210(i1))')it,trj(i)%cm_pc(:)
      read(40,*)it,trj(i)%cm_pc_up(:),trj(i)%cm_pc_low
      read(41,*)it,trj(i)%cm_pg_up(:),trj(i)%cm_pg_low
!      read(43,*)it,trj(i)%cm_pep(:)
!      read(43,'(i6,420(i1))')it,trj(i)%cm_pep(:)
      read(42,*)it,trj(i)%zsc(:),trj(i)%zsc2(:)
      read(43,'(i6,210(i1),210(i1))')it,trj(i)%cm_pep1(:),trj(i)%cm_pep2(:)
!      read(44,*)it,trj(i)%tilt1_6_14, trj(i)%tilt2_6_14
!      read(45,*)it,trj(i)%tilt1_15_20, trj(i)%tilt2_15_20

      read(46,*)it,trj(i)%phiVal1(:),trj(i)%phiVal2(:),trj(i)%psiVal1(:),trj(i)%psiVal2(:)
      read(47,*)it,trj(i)%utilt1_15_20,trj(i)%utilt2_15_20
!      read(48,*)it,trj(i)%utiltH1_15_20,trj(i)%utiltH2_15_20
      read(49,*)it,trj(i)%utilt1_6_14,trj(i)%utilt2_6_14

      read(50,*)dummy(1:6),trj(i)%ro_12(1),dummy(8:9),trj(i)%ro_12(2)
      read(51,*)trj(i)%ro_15(1),dummy(2:4),trj(i)%ro_19(1),dummy(5:5),trj(i)%ro_15(2),trj(i)%ro_19(2)
      read(52,*)dummy(1:6),trj(i)%ro2_12(1),dummy(8:9),trj(i)%ro2_12(2)
!      print*, "READLO HI 12 ", i, trj(i)%ro2_12(1), trj(i)%ro2_12(2)
      read(53,*)trj(i)%ro2_15(1),dummy(2:4),trj(i)%ro2_19(1),dummy(5:5),trj(i)%ro2_15(2),trj(i)%ro2_19(2)
!      print*, "READLO HI 15, 19 ", i, trj(i)%ro2_15(1), trj(i)%ro2_19(1), trj(i)%ro2_15(2),trj(i)%ro2_19(2)
      read(54,*)it,trj(i)%termN_Z(1), trj(i)%termC_Z(1), trj(i)%termN_Z(2), trj(i)%termC_Z(2), trj(i)%termM_Z(1), trj(i)%termM_Z(2), trj(i)%termAll_Z(1), trj(i)%termAll_Z(2)
      read(55,*)trj(i)%LRC1
      read(56,*)trj(i)%LRC2
      read(57,*)trj(i)%isHel1
      read(58,*)trj(i)%isHel2
      read(59,*)trj(i)%scdPg(:),dummy(1:13)
      read(60,*)trj(i)%scdPc(:),dummy(14:26)
      read(61,*)trj(i)%scdPgReg(:),dummy(4:6)
      read(62,*)trj(i)%scdPcReg(:),dummy(4:6)
      read(63,*)trj(i)%scdPgRegCnt(:),dummy(4:6)
      read(64,*)trj(i)%scdPcRegCnt(:),dummy(4:6)
      read(65,*)trj(i)%regionDen(:)
      read(66,*)trj(i)%lip_coord(:)

      trj(i)%bond = bondVal
      trj(i)%angle = angleVal
      trj(i)%dihed = dihedVal
      trj(i)%imprp = imprpVal
      trj(i)%elect = electVal
      trj(i)%vdw = vdwVal
      trj(i)%bondP = bondPVal
      trj(i)%angleP = anglePVal
      trj(i)%dihedP = dihedPVal
      trj(i)%imprpP = imprpPVal
      trj(i)%electP = electPVal
      trj(i)%vdwP = vdwPVal

      trj(i)%bondIon = bondValIon
      trj(i)%angleIon = angleValIon
      trj(i)%dihedIon = dihedValIon
      trj(i)%imprpIon = imprpValIon
      trj(i)%electIon = electValIon
      trj(i)%vdwIon = vdwValIon
      trj(i)%bondPIon = bondPValIon
      trj(i)%anglePIon = anglePValIon
      trj(i)%dihedPIon = dihedPValIon
      trj(i)%imprpPIon = imprpPValIon
      trj(i)%electPIon = electPValIon
      trj(i)%vdwPIon = vdwPValIon


      trj(i)%rg1Int = rg1IntVal
      trj(i)%rg2Int = rg2IntVal
      trj(i)%rg1 = rg1Val
      trj(i)%rg2 = rg2Val

      trj(i)%vol = volVal
      trj(i)%kinetic = kineticVal
      trj(i)%emm = EmmVal
      trj(i)%emmP = EmmPVal
      trj(i)%emmIon = EmmValIon
      trj(i)%emmPIon = EmmPValIon

!      read(36,*) dummy(1:1), SasaVal
      trj(i)%sasa = SasaVal
    enddo
!    print*, "DONE TILT "
!    flush(6)
    close(21)
    close(22)
    close(23)
    close(24)
    close(27)
    close(28)
    close(29)
    close(35)
!    close(36)
    close(37)
    close(38)
!    close(39)
    close(40)
    close(41)
    close(42)
    close(43)
!    close(44)
!    close(45)
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
!    close(59)
!    close(60)
    close(61)
    close(62)
    close(63)
    close(64)
    close(65)
    close(66)
    close(67)
    close(68)

    if(minval(trj%ent) < emin) emin = minval(trj%ent)
    if(maxval(trj%ent) > emax) emax = maxval(trj%ent)
    print*,'Enthalpy min/max:',emin,emax

!    if(minval(trj%entb) < eminb) eminb = minval(trj%entb)
!    if(maxval(trj%entb) > emaxb) emaxb = maxval(trj%entb)
!    print*,'Biased Enthalpy min/max:',eminb,emaxb

!    if(minval(trj%ents) < emins) emins = minval(trj%ents)
!    if(maxval(trj%ents) > emaxs) emaxs = maxval(trj%ents)
!    print*,'Scaled Enthalpy min/max:',emins,emaxs

    print*, "MIN ENT ", minval(trj%ent)
    print*, "MAX ENT ", maxval(trj%ent)
    print*, "Write to structs ", istr
    flush(6)
    do k=1, tprod
      istr = istr + 1

      en(1,istr) = trj(k)%enp
      en(2,istr) = trj(k)%enpw
      en(3,istr) = trj(k)%enw
      en(4,istr) = trj(k)%entot
      en(5,istr) = trj(k)%v

!      print*, "EN ", istr, en(:,istr)
      heppw(1,istr) = trj(k)%eppw

      hcm_pc_up(1:ncont_leaf,istr) = trj(k)%cm_pc_up
      hcm_pc_low(1:ncont_leaf,istr) = trj(k)%cm_pc_low
      hcm_pg_up(1:ncont_leaf,istr) = trj(k)%cm_pg_up
      hcm_pg_low(1:ncont_leaf,istr) = trj(k)%cm_pg_low
      hcm_pep1(1:ncont_pep,istr) = trj(k)%cm_pep1
      hcm_pep2(1:ncont_pep,istr) = trj(k)%cm_pep2
  
      trj(k)%ss=0
      trj(k)%st=0
      trj(k)%shel=0
      trj(k)%sstr=0
      trj(k)%sran=0
      trj(k)%ss2=0
      trj(k)%st2=0
      trj(k)%shel2=0
      trj(k)%sstr2=0
      trj(k)%sran2=0
      do m=1,nmon
       trj(k)%ss(m)=0
       trj(k)%st(m)=0
       trj(k)%shel(m)=0
       trj(k)%sstr(m)=0
       trj(k)%sran(m)=0
       trj(k)%ss2(m)=0
       trj(k)%st2(m)=0
       trj(k)%shel2(m)=0
       trj(k)%sstr2(m)=0
       trj(k)%sran2(m)=0
       j=0
       select case(trj(k)%sschar1(m))
        case('G','H','I')
         j=1
        case('T')
         j=2
        case('C')
         j=3
        case('S')
         j=4
        case('E','B','b')
         j=5
       end select
       if(j.eq.0) then
        print*, 'Error: SS type not selected.'
        stop
       endif
       trj(k)%ss(m)=j
       if(j.eq.2) then       ! turn
        trj(k)%st(m)=1
       else if(j.eq.1) then  ! helix
        trj(k)%shel(m)=1
       else if(j.eq.5) then  ! strand
        trj(k)%sstr(m)=1
       else
        trj(k)%sran(m)=1
       endif

       select case(trj(k)%sschar2(m))
        case('G','H','I')
         j=1
        case('T')
         j=2
        case('C')
         j=3
        case('S')
         j=4
        case('E','B','b')
         j=5
       end select
       if(j.eq.0) then
        print*, 'Error: SS type not selected.'
        stop
       endif
       trj(k)%ss2(m)=j
       if(j.eq.2) then       ! turn
        trj(k)%st2(m)=1
       else if(j.eq.1) then  ! helix
        trj(k)%shel2(m)=1
       else if(j.eq.5) then  ! strand
        trj(k)%sstr2(m)=1
       else
        trj(k)%sran2(m)=1
       endif

      enddo
      termN_Z1(istr) = trj(k)%termN_Z(1)
      termN_Z2(istr) = trj(k)%termN_Z(2)
      termC_Z1(istr) = trj(k)%termC_Z(1)
      termC_Z2(istr) = trj(k)%termC_Z(2)
      termM_Z1(istr) = trj(k)%termM_Z(1)
      termM_Z2(istr) = trj(k)%termM_Z(2)
      termAll_Z1(istr) = trj(k)%termAll_Z(1)
      termAll_Z2(istr) = trj(k)%termAll_Z(2)



      tilt1_1_8(istr) = trj(k)%tilt1_1_8
      tilt1_8_14(istr) = trj(k)%tilt1_8_14
      tilt1_14_20(istr) = trj(k)%tilt1_14_20
      tilt1_1_11(istr) = trj(k)%tilt1_1_11
      tilt1_11_21(istr) = trj(k)%tilt1_11_21
      tilt2_1_8(istr) = trj(k)%tilt2_1_8
      tilt2_8_14(istr) = trj(k)%tilt2_8_14
      tilt2_14_20(istr) = trj(k)%tilt2_14_20
      tilt2_1_11(istr) = trj(k)%tilt2_1_11
      tilt2_11_21(istr) = trj(k)%tilt2_11_21

      utilt1_15_20(istr) = trj(k)%utilt1_15_20
      utiltH1_15_20(istr) = trj(k)%utiltH1_15_20
      utilt1_6_14(istr) = trj(k)%utilt1_6_14
      utilt2_15_20(istr) = trj(k)%utilt2_15_20
      utiltH2_15_20(istr) = trj(k)%utiltH2_15_20
      utilt2_6_14(istr) = trj(k)%utilt2_6_14

      hst(1:nmon,istr) = trj(k)%st
      hstr(1:nmon,istr) = trj(k)%sstr
      hhel(1:nmon,istr) = trj(k)%shel
      hran(1:nmon,istr) = trj(k)%sran
      totVal = hst(1:nmon,istr) + hstr(1:nmon,istr) + hhel(1:nmon,istr) + hran(1:nmon,istr)
      hst2(1:nmon,istr) = trj(k)%st2
      hstr2(1:nmon,istr) = trj(k)%sstr2
      hhel2(1:nmon,istr) = trj(k)%shel2
      hran2(1:nmon,istr) = trj(k)%sran2
      totVal2 = hst2(1:nmon,istr) + hstr2(1:nmon,istr) + hhel2(1:nmon,istr) + hran2(1:nmon,istr)
      sasa(istr) = trj(k)%sasa

      emm(istr) = trj(k)%emm
      emmP(istr) = trj(k)%emmP
      bond(istr) = trj(k)%bond
      angle(istr) = trj(k)%angle
      dihed(istr) = trj(k)%dihed
      imprp(istr) = trj(k)%imprp
      elect(istr) = trj(k)%elect
      vdw(istr) = trj(k)%vdw
      vol(istr) = trj(k)%vol
      kinetic(istr) = trj(k)%kinetic
      bondP(istr) = trj(k)%bondP
      angleP(istr) = trj(k)%angleP
      dihedP(istr) = trj(k)%dihedP
      imprpP(istr) = trj(k)%imprpP
      electP(istr) = trj(k)%electP
      vdwP(istr) = trj(k)%vdwP

      emmIon(istr) = trj(k)%emmIon
      emmPIon(istr) = trj(k)%emmPIon
      bondIon(istr) = trj(k)%bondIon
      angleIon(istr) = trj(k)%angleIon
      dihedIon(istr) = trj(k)%dihedIon
      imprpIon(istr) = trj(k)%imprpIon
      electIon(istr) = trj(k)%electIon
      vdwIon(istr) = trj(k)%vdwIon
      bondPIon(istr) = trj(k)%bondPIon
      anglePIon(istr) = trj(k)%anglePIon
      dihedPIon(istr) = trj(k)%dihedPIon
      imprpPIon(istr) = trj(k)%imprpPIon
      electPIon(istr) = trj(k)%electPIon
      vdwPIon(istr) = trj(k)%vdwPIon

      rg1(istr) = trj(k)%rg1
      rg2(istr) = trj(k)%rg2
      rg1Int(istr) = trj(k)%rg1Int
      rg2Int(istr) = trj(k)%rg2Int
      LRC1(istr) = trj(k)%LRC1
      LRC2(istr) = trj(k)%LRC2
      isHel1(istr) = trj(k)%isHel1
      isHel2(istr) = trj(k)%isHel2
      zsc(1:nmon,istr) = trj(k)%zsc
      zsc2(1:nmon,istr) = trj(k)%zsc2
      scdPg(1:13,istr)=trj(k)%scdPg
      scdPc(1:13,istr)=trj(k)%scdPc
      scdPgReg(1:3,istr)=trj(k)%scdPgReg
      print*, "MOVE PG ", k, istr, scdPgReg(1,istr), scdPgReg(2,istr), scdPgReg(3,istr)
      scdPcReg(1:3,istr)=trj(k)%scdPcReg
      scdPgRegCnt(1:3,istr)=trj(k)%scdPgRegCnt
      scdPcRegCnt(1:3,istr)=trj(k)%scdPcRegCnt
      regionDen(1:12,istr)=trj(k)%regionDen
      lip_coord(1:4,istr)=trj(k)%lip_coord
      ro_12(1:2,istr) = trj(k)%ro_12
      ro2_12(1:2,istr) = trj(k)%ro2_12
      ro_15(1:2,istr) = trj(k)%ro_15
      ro2_15(1:2,istr) = trj(k)%ro2_15
      ro_19(1:2,istr) = trj(k)%ro_19
      ro2_19(1:2,istr) = trj(k)%ro2_19
      phiVal1(1:nmon-2,istr) = trj(k)%phiVal1
      phiVal2(1:nmon-2,istr) = trj(k)%phiVal1
      psiVal1(1:nmon-2,istr) = trj(k)%psiVal1
      psiVal2(1:nmon-2,istr) = trj(k)%psiVal1
!      print*, "DONE ", istr, rep
      flush(6)

    enddo  !ending structures for rep
!    emean=0.0
!    do k=1,istr
!     emean = emean+en(1,istr)+en(2,istr)
!    enddo
!    emean=emean/istr
!    print*,emean,istr
   enddo
   deallocate(trj)
 enddo
!    prev was enddo for each traj



 do i=1,istr
  !Sanity checks?

 enddo

! emean=0.0
! do k=1,istr
!  emean = emean+en(1,k)+en(2,k)
! enddo
! emean=emean/istr
! print*,emean,istr


   if (in_tr < 10) then 
      open(34,file='histcm_pc1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histcm_pc1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hcm_pc_up(:,i)
   enddo
  close(34)
   if (in_tr < 10) then 
      open(34,file='histcm_pc2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histcm_pc2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hcm_pc_low(:,i)
   enddo
  close(34)

   if (in_tr < 10) then
      open(34,file='histcm_pg1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histcm_pg1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hcm_pg_up(:,i)
   enddo
  close(34)
   if (in_tr < 10) then
      open(34,file='histcm_pg2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histcm_pg2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hcm_pg_low(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='histcm_pep1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histcm_pep1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hcm_pep1(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='histcm_pep2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histcm_pep2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hcm_pep2(:,i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='zsc1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='zsc1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)zsc(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='zsc2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='zsc2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)zsc2(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='scdPg_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='scdPg_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)scdPg(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='scdPc_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='scdPc_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)scdPc(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='scdPgReg_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='scdPgReg_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     print*, "WRITE PG ", i, scdPgReg(1,i),scdPgReg(2,i),scdPgReg(3,i)
     write(34)scdPgReg(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='scdPcReg_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='scdPcReg_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)scdPcReg(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='scdPgRegCnt_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='scdPgRegCnt_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)scdPgRegCnt(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='scdPcRegCnt_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='scdPcRegCnt_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)scdPcRegCnt(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='regionDen_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='regionDen_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)regionDen(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='lipid_coord_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='lipid_coord_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)lip_coord(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='nTerm1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='nTerm1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termN_Z1(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='nTerm2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='nTerm2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termN_Z2(i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='cTerm1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='cTerm1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termC_Z1(i)
   enddo
   close(34)

  if (in_tr < 10) then
      open(34,file='cTerm2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='cTerm2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termC_Z2(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='midCOM1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='midCOM1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termM_Z1(i)
   enddo
   close(34)

  if (in_tr < 10) then
      open(34,file='midCOM2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='midCOM2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termM_Z2(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='allCOM1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='allCOM1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termAll_Z1(i)
   enddo
   close(34)

  if (in_tr < 10) then
      open(34,file='allCOM2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='allCOM2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)termAll_Z2(i)
   enddo
  close(34)





  if (in_tr < 10) then
      open(34,file='ro12_1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='ro12_1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)ro_12(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='ro12_2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='ro12_2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)ro2_12(:,i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='ro15_1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='ro15_1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)ro_15(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='ro15_2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='ro15_2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)ro2_15(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='ro19_1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='ro19_1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)ro_19(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='ro19_2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='ro19_2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)ro2_19(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='histst'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histst1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hst(:,i)
   enddo
  close(34)
  
  if (in_tr < 10) then
      open(34,file='histst2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histst2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hst2(:,i)
   enddo
  close(34)
  
  if (in_tr < 10) then
      open(34,file='energyi'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='energyi1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)en(:,i)
   enddo
  close(34)

!
!    Save all state info
  if (in_tr < 10) then
      open(34,file='histstr'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histstr1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hstr(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='histstr2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histstr2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hstr2(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='histhel'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histhel1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hhel(:,i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='histhel2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histhel2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hhel2(:,i)
   enddo

  close(34)

   if (in_tr < 10) then
      open(34,file='histran'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histran1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hran(:,i)
   enddo
  close(34)

   if (in_tr < 10) then
      open(34,file='histran2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='histran2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)hran2(:,i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='sasa'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='sasa1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)sasa(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='emm'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='emm1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)emm(i)
   enddo
  close(34)
 
  if (in_tr < 10) then
      open(34,file='emmP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='emmP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)emmP(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='emmIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='emmIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)emmIon(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='emmPIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='emmPIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)emmPIon(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='rg1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='rg1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)rg1(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='rg2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='rg2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)rg2(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='rg1Int_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='rg1Int_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)rg1Int(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='rg2Int_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='rg2Int_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)rg2Int(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='LRC1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='LRC1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)LRC1(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='LRC2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='LRC2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)LRC2(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='isHel1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='isHel1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)isHel1(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='isHel2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='isHel2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)isHel2(i)
   enddo
  close(34)





  if (in_tr < 10) then
      open(34,file='bond'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='bond1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)bond(i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='bondP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='bondP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)bondP(i)
   enddo
  close(34)


if (in_tr < 10) then
      open(34,file='angle'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='angle1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)angle(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='angleP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='angleP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)angleP(i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='dihed'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='dihed1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)dihed(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='dihedP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='dihedP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)dihedP(i)
   enddo
  close(34)
 
if (in_tr < 10) then
      open(34,file='imprp'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='imprp1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)imprp(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='imprpP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='imprpP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)imprpP(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='elect'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='elect1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)elect(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='electP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='electP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)electP(i)
   enddo
  close(34)


if (in_tr < 10) then
      open(34,file='vdw'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='vdw1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)vdw(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='vdwP'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='vdwP1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)vdwP(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='bondIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='bondIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)bondIon(i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='bondPIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='bondPIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)bondPIon(i)
   enddo
  close(34)


if (in_tr < 10) then
      open(34,file='angleIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='angleIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)angleIon(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='anglePIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='anglePIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)anglePIon(i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='dihedIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='dihedIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)dihedIon(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='dihedPIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='dihedPIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)dihedPIon(i)
   enddo
  close(34)
 
if (in_tr < 10) then
      open(34,file='imprpIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='imprpIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)imprpIon(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='imprpPIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='imprpPIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)imprpPIon(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='electIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='electIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)electIon(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='electPIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='electPIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)electPIon(i)
   enddo
  close(34)


if (in_tr < 10) then
      open(34,file='vdwIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='vdwIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)vdwIon(i)
   enddo
  close(34)
  if (in_tr < 10) then
      open(34,file='vdwPIon'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='vdwPIon1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)vdwPIon(i)
   enddo
  close(34)





  if (in_tr < 10) then
      open(34,file='vol'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='vol1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)vol(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='kinetic'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='kinetic1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)kinetic(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='utilt1_15_20_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='utilt1_15_20_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)utilt1_15_20(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='utilt2_15_20_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='utilt2_15_20_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)utilt2_15_20(i)
   enddo
  close(34)

!  if (in_tr < 10) then
!      open(34,file='utiltH1_15_20_'//char(48+in_tr)//'.dat',form='unformatted')
!   else
!      open(34,file='utiltH1_15_20_1'//char(38+in_tr)//'.dat',form='unformatted')
!   endif
!   do i=1,istr
!     write(34)utiltH1_15_20(i)
!   enddo
!  close(34)

!  if (in_tr < 10) then
!      open(34,file='utiltH2_15_20_'//char(48+in_tr)//'.dat',form='unformatted')
!   else
!      open(34,file='utiltH2_15_20_1'//char(38+in_tr)//'.dat',form='unformatted')
!   endif
!   do i=1,istr
!     write(34)utiltH2_15_20(i)
!   enddo
!  close(34)

  if (in_tr < 10) then
      open(34,file='utilt1_6_14_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='utilt1_6_14_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)utilt1_6_14(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='utilt2_6_14_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='utilt2_6_14_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)utilt2_6_14(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt1_1_8_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt1_1_8_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_1_8(i)
   enddo
  close(34)





  if (in_tr < 10) then
      open(34,file='tilt1_8_14_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt1_8_14_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_8_14(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt1_14_20_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt1_14_20_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_14_20(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt1_1_11_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt1_1_11_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_1_11(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt1_11_21_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt1_11_21_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_11_21(i)
   enddo
  close(34)



  if (in_tr < 10) then
      open(34,file='tilt2_1_8_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_1_8_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_1_8(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt2_8_14_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_8_14_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_8_14(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt2_14_20_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_14_20_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_14_20(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt2_1_11_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_1_11_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_1_11(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt2_11_21_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_11_21_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_11_21(i)
   enddo
  close(34)


  if (in_tr < 10) then
      open(34,file='tilt2_6_14_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_6_14_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_6_14(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt_6_14_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt_6_14_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_6_14(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt2_15_20_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt2_15_20_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt2_15_20(i)
   enddo
  close(34)

  if (in_tr < 10) then
      open(34,file='tilt1_15_20_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='tilt1_15_20_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)tilt1_15_20(i)
   enddo
  close(34)

 if (in_tr < 10) then
      open(34,file='phi1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='phi1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)phiVal1(:,i)
   enddo
  close(34)

 if (in_tr < 10) then
      open(34,file='phi2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='phi2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)phiVal2(:,i)
   enddo
  close(34)

 if (in_tr < 10) then
      open(34,file='psi1_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='psi1_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)psiVal1(:,i)
   enddo
  close(34)

 if (in_tr < 10) then
      open(34,file='psi2_'//char(48+in_tr)//'.dat',form='unformatted')
   else
      open(34,file='psi2_1'//char(38+in_tr)//'.dat',form='unformatted')
   endif
   do i=1,istr
     write(34)psiVal2(:,i)
   enddo
  close(34)







end program hist1z


