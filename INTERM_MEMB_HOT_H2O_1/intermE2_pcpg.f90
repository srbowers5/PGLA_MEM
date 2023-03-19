! Set na (number atoms)
!     na_tot (number atoms)
!     nmon (aa/pep)
!     nmol (num pep)
!     nlip (num lipids)
!     nNa, nCl, nWater
!     nla  (lipid atoms? should match count(.not.protein))
!     nha  (num hydrogen atoms)

module parameters
 implicit none
 integer, parameter :: not_now=1
 integer, parameter :: short = selected_int_kind(2) 
 integer, parameter :: na=25184,nmon=21,nmol=2,npmon=nmon*nmol, nlip=98, first_pep=11299, &
                       nNa=38, nCl=10, nwater=4410,   &
                       nmontot=npmon+nlip+nNa+nCl+nwater, nha=15982,nla=24576,na_tot=25184, &
                       itot_rg=3*(nmon**2-nmon),                           &
                       runmin=3,runmax=3,repmin=1,repmax=1,&
                       run1=3,run2=3,rep1=1,rep2=1, nrun=1,nrep=1,numstr=$SIZE_RUN
 integer, parameter :: NTERM_START=1,NTERM_END=11, CTERM_START=12,CTERM_END=21,MID_START=6
 integer, dimension(runmin:runmax,repmin:repmax), parameter:: &
          strnum=reshape([numstr],[nrun,nrep])

!...nha - # of any H atoms in all peptides
 real, parameter :: pi = 3.1415926535898, eps=1.e-6
 real, parameter :: rc=4.6, rc_cm = 6.5, rc_NO=3.5, alpha=2.*pi/3.,rDAmin=2.9,&
                    Ehb0=-0.5
 character(len=9), parameter :: HBdefinition = 'geometric'
 character(len=256), parameter:: &
                  path='$PATH'

 character(len=16):: stepString
 character(len=16):: repString
 character(len=16):: runString


 type ProteinStructure
  integer            :: atomnum
  real, dimension(3) :: pos
  type (ProteinStructure), pointer :: next
 end type ProteinStructure
 type (ProteinStructure), pointer :: stpi
 type (ProteinStructure), dimension(:,:), pointer  :: stp

 real, dimension(numstr,3) :: box 
 logical, dimension(nlip) :: isPC
 logical,  dimension(na) :: backbone, terminal,hydro,protein,& 
                            choline,phg,glycerol,fa1,fa2,ionatom,wateratom
 logical,  dimension(nmontot) :: chainterminal, peptide, lipid,water,ionNa,ionCl
 character(len=4),  dimension(na) :: atomname
 character(len=4),  dimension(nmontot) :: sequence
 integer, dimension(na) :: resid
 integer, dimension(nmontot) :: chain, charge, hp
 real, dimension(numstr,npmon,3) :: posCa, posN, posC, posO, posH,poss
 real, dimension(numstr,npmon,3) :: term_com(numstr,6)
 real, dimension(numstr,nwater,3) :: posW
 real, dimension(numstr,nNa,3) :: posNai
 real, dimension(numstr,nCl,3) :: posCli
 real, dimension(numstr,nlip,5,3) :: posLP
 real, dimension(numstr) :: z0
 integer :: aanum
 integer :: chain1, chain2
 real :: nterm1_com, nterm2_com, cterm1_com, cterm2_com
 integer :: nterm1_cnt, nterm2_cnt, cterm1_cnt, cterm2_cnt

 real :: mid1_com, mid2_com
 integer :: mid1_cnt, mid2_cnt
 character(len=80) :: dcdfile
 character(len=80) :: pdbfile
 integer :: run,ipim, ippi
 real :: plim, fcbl

!  Input parameters
integer :: runnum, first_struct, tr, rep

end module parameters



program interm
 use parameters
 implicit none
 integer :: i,j,k,l, n, str
 integer :: vals, tens
 type (ProteinStructure), pointer :: stpii

 read*, tr, rep, runnum, first_struct
 print*, "TR ", tr, "REP", rep, "RUN ", runnum, "FIRST_ST ", first_struct

! do rep=rep1,rep2

  if (rep < 10) then
      repString = char(48+rep)
  else
      tens = rep/10
      vals = rep-(tens*10)
      repString = char(48+tens)//char(48+vals)
  endif

  if (runnum < 10) then
      runString = char(48+runnum)
  else
      tens = runnum/10
      vals = runnum-(tens*10)
      runString = char(48+tens)//char(48+vals)
  endif

   print*,"PATH (", path, ")"

   open(24,file=trim(path)//'interm_out/conthb'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(25,file=trim(path)//'interm_out/r'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(27,file=trim(path)//'interm_out/contlist'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(28,file=trim(path)//'interm_out/cm'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(31,file=trim(path)//'interm_out/dhWT'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(35,file=trim(path)//'interm_out/hb'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(36,file=trim(path)//'interm_out/hblist_'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(38,file=trim(path)//'interm_out/hbres_'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')

   open(43,file=trim(path)//'interm_out/contlist'//char(48+tr)//'_'//trim(repString)//'_lpb'//trim(runString)//'.dat')
   open(45,file=trim(path)//'interm_out/llist'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(49,file=trim(path)//'interm_out/aabd'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(50,file=trim(path)//'interm_out/cm'//char(48+tr)//'_'//trim(repString)//'_lpb'//trim(runString)//'.dat')
   open(51,file=trim(path)//'interm_out/zsc'//char(48+tr)//'_'//trim(repString)//'b'//trim(runString)//'.dat')
   open(52,file=trim(path)//'interm_out/contlist'//char(48+tr)//'_'//trim(repString)//'_lpPCb'//trim(runString)//'.dat')
   open(53,file=trim(path)//'interm_out/contlist'//char(48+tr)//'_'//trim(repString)//'_lpPGb'//trim(runString)//'.dat')
   open(54,file=trim(path)//'interm_out/term_com'//char(48+tr)//'_'//trim(repString)//'_'//trim(runString)//'.dat')


  ipim=0; plim=0.0; fcbl=0.0; ippi=0
  do run=run1, run2

  if(strnum(run,rep)==0)cycle
  allocate(stp(strnum(run,rep),nmontot))

!  print*, "calling ", run,  rep
  call readin
  call unitcell
  call coarsegrain
  call structure
  call WTDihedral
  call Hbonds

! releasing the memory associated with structures tagged by pointers ... 
  do str=1,strnum(run,rep)
    do n=1,nmontot
     stpi => stp(str,n); l=0
     do while(associated(stpi%next))
        stpii => stpi
        stpi => stpi%next
        if(l>0)deallocate(stpii)
        l = l + 1
     enddo
     deallocate(stpi)
    enddo
   enddo
  deallocate(stp)

  enddo


  close(24)
  close(25)
  close(27)
  close(28)
  close(31)
  close(35)
  close(36)
  close(38)
  close(43)
  close(45)
  close(49)
  close(50)
  close(51)
  close(52)
  close(53)
  close(54)
 
  print*,'Fraction of structures with peptide image interactions:', &
          real(ipim)/sum(strnum(run1:run2,rep))
  print*,'Lipid fraction crossbridging ab images:', plim/sum(strnum(run1:run2,rep))
  print*,'Lipid fraction crosslinked by ab1 and ab2:',fcbl/sum(strnum(run1:run2,rep))
  print*,'Fraction of structures with peptide-peptide interaction through BL:', &
          real(ippi)/sum(strnum(run1:run2,rep))


! enddo

end program interm

!...........................................................................

subroutine readin
 use parameters
 implicit none
 integer :: i, j, k, l,m,n, q,atomnum, resnum,str,istr,icount,iresnum,iresnum0,&
            nstr, ntitle, na0, ic, nacount, stri
 integer, dimension(9) :: dumi
 real, dimension(3) :: coord
 real(8), dimension(6) :: ibox
 real, dimension(na_tot,3) :: ipos
 real :: dumr
 character(len=4) :: residue,atom, segment,segment0,dcdhdr
 character(len=1) :: chn
 character(len=6)  :: lineStart
 character(len=3)  :: endRecord
 character(len=80), dimension(2) :: title
 integer :: hthou, tthou, thou, hun, tens, tmp_stri


!    print*, "readit ", run, rep
!...initializing.............
   backbone = .false.; terminal = .false.; hydro = .false. 
   chainterminal=.false.; protein = .false.; lipid= .false.
   peptide=.false.; choline=.false.; phg=.false.; glycerol=.false.
   fa1=.false.; fa2=.false.
   sequence = 'xxxx'; atomname = 'xxxx'
   chain=0; charge=0; hp=0; resid=0
   nullify(stpi)
   do str=1,strnum(run,rep)
    do n = 1, nmontot
       nullify(stp(str,n)%next)
    enddo
   enddo   

!.......read topology file 
100 format(a6,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3,18x,a4,3x)
   pdbfile=trim(path)//'/spt.pdb'
   open(21,file=pdbfile)

    read(21,'()')
    iresnum0=0; resnum = 0; segment0='xxxx'
    do n = 1, na
      read(21,100)lineStart,atomnum,atom,residue,chn,iresnum,coord,segment
      if(n/=atomnum)then
        print*,'Error: n /= atomnum'; stop
      endif 
      if( (iresnum/=iresnum0).or.(segment/=segment0) ) then
        iresnum0=iresnum; segment0=segment
        resnum = resnum + 1
        ic=0
      endif
      ic = ic + 1

      atomname(atomnum) = atom
      resid(atomnum) = resnum
!      if(( residue/='DMPG' ).and.( residue/='DMPE' ).and.( residue/='DMPC' ).and.( residue/='CLA ' ).and.( residue/='SOD ' ).and.( residue/='TIP3' )) then
      if(( residue/='DMPG' ).and.( residue/='DMPE' ).and.( residue/='DMPC' ).and.( residue/='CLA ' ).and. &
        ( residue/='SOD ' ).and.( residue/='TIP3' )) then
         protein(atomnum) = .true.
         if(atom==' N  ')then
              backbone(atomnum) = .true.
         elseif(atom==' CA ')then
              backbone(atomnum) = .true.
         elseif(atom==' C  ')then
              backbone(atomnum) = .true.
         elseif(atom==' O  ')then
              backbone(atomnum) = .true.
         elseif(atom==' HN ')then
              backbone(atomnum) = .true.
         elseif(atom==' HA ')then
              backbone(atomnum) = .true.
         elseif(atom==' HA1')then
              backbone(atomnum) = .true.
         elseif(atom==' HA2')then
              backbone(atomnum) = .true.
         endif

         if(atom==' CAY')then
              terminal(atomnum) = .true.
         elseif(atom==' HY1')then
              terminal(atomnum) = .true.
         elseif(atom==' HY2')then
              terminal(atomnum) = .true.
         elseif(atom==' HY3')then
              terminal(atomnum) = .true.
         elseif(atom==' CY ')then
              terminal(atomnum) = .true.
         elseif(atom==' OY ')then
              terminal(atomnum) = .true.
         elseif(atom==' NT ')then
              terminal(atomnum) = .true.
         elseif(atom==' HT1')then
              terminal(atomnum) = .true.
         elseif(atom==' HT2')then
              terminal(atomnum) = .true.
         elseif(atom==' HT3')then
              terminal(atomnum) = .true.
         endif
         if(terminal(atomnum)) chainterminal(resnum)=.true. 

         if( protein(atomnum).and.(atom == ' CA ') ) then
          peptide(resnum) = .true.
          sequence(resnum) = residue
          chain(resnum)   = ((resnum-(nlip+1))/nmon) + 1
          if( (residue=='LYS ').or.(residue=='ARG ') )charge(resnum)= 1
          if( (residue=='GLU ').or.(residue=='ASP ') )charge(resnum)=-1
          if( (residue=='VAL ').or.(residue=='ILE ').or.(residue=='LEU ').or.&
              (residue=='PHE ').or.(residue=='TRP ').or.(residue=='CYS ').or.&
              (residue=='ALA ').or.(residue=='MET ').or.(residue=='PRO ') )&
              hp(resnum)=1
          endif
      endif
      if((index(atom,'H')==1).or.(index(atom,'H')==2))hydro(atomnum) = .true.

      if( residue=='DMPC')then
        if(atom == ' P  ') then
          sequence(resnum) = residue
          chain(resnum) = resnum
          lipid(resnum)= .true.
          isPC(iresnum)= .true.
        endif
        if( (ic <= 19).or.(ic==24) )choline(atomnum)=.true.
        if( (ic>=20).and.(ic<=22) )phg(atomnum)=.true.
        if( (ic==23).or.((ic >=25).and.(ic<=30)).or.((ic >=36).and.(ic<=39)))&
              glycerol(atomnum)=.true.
        if(((ic >=31).and.(ic<=35)).or.((ic>=45).and.(ic<=81)))fa1(atomnum)=.true.
        if(((ic >=40).and.(ic<=44)).or.((ic>=82).and.(ic<=118)))fa2(atomnum)=.true.
      endif


      if(residue=='DMPE')then 
        if(atom == ' P  ') then
          sequence(resnum) = residue
          chain(resnum) = resnum+nmol
          lipid(resnum)= .true.
        endif
!
!  CHECK
        if( (ic <= 10).or.(ic==15) )choline(atomnum)=.true.
        if( (ic>=11).and.(ic<=13) )phg(atomnum)=.true.
        if( (ic==14).or.((ic >=16).and.(ic<=21)).or.((ic >=27).and.(ic<=30)))&
              glycerol(atomnum)=.true.
        if(((ic >=22).and.(ic<=26)).or.((ic>=36).and.(ic<=72)))fa1(atomnum)=.true.
        if(((ic >=31).and.(ic<=35)).or.((ic>=73).and.(ic<=109)))fa2(atomnum)=.true.
      endif


      if(residue=='DMPG')then
        if(atom == ' P  ') then
          sequence(resnum) = residue
          chain(resnum) = resnum+2
          lipid(resnum)= .true.
          isPC(iresnum)= .false.
        endif

        if( (ic <= 12).or.(ic==17) )choline(atomnum)=.true.
        if( (ic>=13).and.(ic<=15) )phg(atomnum)=.true.
        if( (ic==16).or.((ic >=18).and.(ic<=23)).or.((ic >=29).and.(ic<=32)))&
              glycerol(atomnum)=.true.
        if(((ic >=24).and.(ic<=28)).or.((ic>=38).and.(ic<=74)))fa1(atomnum)=.true.
        if(((ic >=33).and.(ic<=37)).or.((ic>=75).and.(ic<=111)))fa2(atomnum)=.true.
      endif
!
!  END CHECK


      if( residue=='TIP3')then
        if(atom == ' OH2') then
          sequence(resnum) = residue
          chain(resnum) = resnum - npmon + nmol
          water(resnum)= .true.
        endif
        wateratom(atomnum)=.true.
      endif

      if( residue=='SOD ')then
          sequence(resnum) = residue
          chain(resnum) = resnum - npmon + nmol
          ionNa(resnum)= .true.
          ionatom(atomnum)=.true.
      endif

      if( residue=='CLA ')then
          sequence(resnum) = residue
          chain(resnum) = resnum - npmon + nmol
          ionCl(resnum)= .true.
          ionatom(atomnum)=.true.
      endif

!        print*,resnum,sequence(resnum),chain(resnum),charge(resnum), hp(resnum)
    enddo
    read(21,'(a3)')endRecord
    if(endRecord/='END')then 
       print*,'Error in reading pdb file',endRecord; stop
    endif

    if( nmontot /= resnum ) then
      print*,'Error: nmontot and resnum do not match', nmontot, resnum; stop
    endif
    if( nlip /= count(lipid) ) then
      print*,'Error: nlip and count(lipid) do not match',nlip,count(lipid)
      stop
    endif
! Error: nla is incorrect       11346       24576

    if( nla /= count(.not.protein) ) then
      print*,'Error: nla is incorrect',nla,count(.not.protein); stop
    endif
   close(21)
   print*,'Topology file read'


!...reading dcd coordinates...........
   do str=1,numstr

    stri=str+first_struct-1
    tmp_stri = stri
    hthou = tmp_stri / 100000
    tmp_stri = tmp_stri - (hthou * 100000)
    tthou = tmp_stri / 10000
    tmp_stri = tmp_stri - (tthou * 10000)
    thou = tmp_stri / 1000
    tmp_stri = tmp_stri - (thou * 1000)
    hun = tmp_stri / 100
    tmp_stri = tmp_stri - (hun * 100)
    tens = tmp_stri / 10
    tmp_stri = tmp_stri - (tens * 10)

    if (hthou > 0) then
        stepString = char(hthou+48)//char(tthou+48)//char(thou+48)//char(hun+48)//char(tens+48)//char(tmp_stri+48)
    else
        stepString = char(tthou+48)//char(thou+48)//char(hun+48)//char(tens+48)//char(tmp_stri+48)
    endif

    dcdfile=trim(path)//'/output/pgl22_mem'//char(48+tr)//'_'//trim(stepString)//'_'//trim(repString)//'.dcd'


   open(21,file=trim(dcdfile),form='unformatted')
!   print*,"READ", dcdfile
!   flush(6)
   print*,'reading str = ', str, stri, rep
   flush(6)

   read(21) dcdhdr, nstr, dumi(1:8), dumr, dumi(1:9)
   if(nstr/=1) then
    print*, 'Error: nstr /= 1', nstr; stop
   endif
   read(21) ntitle, title(1:ntitle)
   if(ntitle/=2) then
       print*, 'Error: ntitle /= 2', ntitle; stop
   endif
   read(21)na0
   if(na0/=na_tot) then
     print*, 'Error: na0 /= na_tot', na0; stop
   endif

   read(21) ibox
!   print*,rep,str,ibox(1),ibox(3),ibox(6)
   read(21) ipos(:,1)
   read(21) ipos(:,2)
   read(21) ipos(:,3)
   box(str,:)=(/ibox(1),ibox(3),ibox(6)/)
   do i=1,na
    stpi => stp(str,resid(i))
    do while(associated(stpi%next))
      stpi => stpi%next
    enddo
    stpi%atomnum = i
    stpi%pos = ipos(i,:)
    allocate(stpi%next); stpi => stpi%next; nullify(stpi%next)
   enddo
   close(21)

  enddo
  print*,'Trajectory segment loaded....',tr,run,rep

  icount=count(hydro .eqv. .true.)
  if( icount /= nha)then
       print*,'Error: wrong number of H',icount,nha; stop
  endif

   if( (maxval(chain)/=nmol+nlip+nwater+nNa+nCl).or.(minval(chain)/=1) )then 
    print*,'Error:  wrong chain values', minval(chain), maxval(chain); stop
   endif



!
!  Checks not valid
!!!  do str=1,strnum(run,rep)
!!!   do n=1,nmontot
!!!    stpi => stp(str,n)
!!!    nacount = 0
!!!    do while(associated(stpi%next))
!!!      nacount = nacount + 1      
!!!      stpi => stpi%next
!!!    enddo
!!!    if(n>npmon.and.n<121) then
!!!      if(nacount/=118)then
!!!        print*,'Error: wrong na in lipid',n; stop
!!!      endif
!!!    endif
!!!    if(n==2.or.n==13) then
!!!      if(nacount/=11)then
!!!        print*,'Error: wrong na in ser', nacount; stop
!!!      endif 
!!!    endif
!!!   enddo
!!!  enddo

end subroutine readin

!...........................................................................

subroutine unitcell
 use parameters
 implicit none
 integer :: n, i, m, str, ix,iy,iz, l
 integer, dimension(nmol+nlip+nwater+nNa+nCl) :: ina
 real, dimension(3) :: shift, cmposi
 real, dimension(3,nmol+nlip+nwater+nNa+nCl) :: cmpos, shiftmin
 real :: d, rmin


  do str = 1,strnum(run,rep)

!   print*, "UNIT CELL ", str
    cmpos=0.0; ina=0
    do n=1,nmontot

     stpi => stp(str,n)
     do while(associated(stpi%next))
      if( .not.hydro(stpi%atomnum) ) then
                   cmpos(:,chain(n)) = cmpos(:,chain(n)) + stpi%pos
                   ina(chain(n)) = ina(chain(n)) + 1
      endif
      stpi => stpi%next
     enddo

    enddo
    do i=1,nmol+nlip+nwater+nNa+nCl
      cmpos(:,i) = cmpos(:,i) / ina(i)
    enddo
    if(sum(ina) /= na-nha)then
      print*,'Error in unitcell: sum(ina) /= na-nha ',sum(ina),na-nha; stop
    endif

!....translating back to the unit cell...

    do l=1,nmol+nlip+nwater+nNa+nCl
      do m=1,3
        shiftmin(m,l) = box(str,m)*nint(cmpos(m,l)/box(str,m))
      enddo
    enddo

    do i=1,nmol+nlip+nwater+nNa+nCl
      do m=1,3 
        if(abs(nint(shiftmin(m,i)/box(str,m))) > 1)then 
         print*,'Some chains are out of the box',nint(shiftmin(m,i)/box(str,m)); stop
        endif 
      enddo
      cmpos(:,i) = cmpos(:,i) - shiftmin(:,i)
    enddo
    do n=1,nmontot
     stpi => stp(str,n)
     do while(associated(stpi%next))
         stpi%pos = stpi%pos - shiftmin(:,chain(n))
         stpi => stpi%next
     enddo
    enddo
!.......................................

    do l=1,nmol+nlip+nwater+nNa+nCl
     do m=1,3
      if( abs(cmpos(m,l)) > box(str,m)/2.+ 5.)then 
        print*,'Error: peptide/lipid outside the unitcell',str,l,cmpos(:,l),box(str,:); stop
      endif
     enddo
    enddo

!.....placing ligands around the peptide.....
!
!    shiftmin=0.d0
!    do l=nmol+1,nmol+nlip    
!
!       rmin=100000.
!       do ix = -1,1
!       do iy = -1,1
!       do iz = -1,1
!        shift(1) = ix*box(1); shift(2) = iy*box(2); shift(3) = iz*box(3)
!        cmposi = cmpos(:,l) + shift
!        do n=1,npmon
!          stpi => stp(str,n)
!          do while(associated(stpi%next))
!           if( .not.hydro(stpi%atomnum) ) then
!             d = sum( (cmposi-stpi%pos)**2 )
!             if(d<rmin)then 
!               rmin=d 
!               shiftmin(:,l)=shift
!             endif
!           endif
!           stpi => stpi%next
!          enddo
!        enddo
!       enddo
!       enddo
!       enddo   
!    enddo

!    do n=npmon+1,nmontot
!     stpi => stp(str,n)
!     do while(associated(stpi%next))
!         stpi%pos = stpi%pos + shiftmin(:,chain(n))
!         stpi => stpi%next
!     enddo
!    enddo

  enddo

!...output translated structure....
!100 format(a6,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3,18x,a4,3x)
!  if(rep < 10) then
!   pdbfilei='./COL'//char(48+run)//'/abf_quench'//char(48+tr)//'_'//char(48+rep)//'_tr_translated.pdb'
!  elseif(rep >= 10)then
!   pdbfilei='./COL'//char(48+run)//'/abf_quench'//char(48+tr)//'_'//char(48+rep/10)// &
!                          char(48+rep-rep/10*10)//'_tr_translated.pdb'
!  endif

!  open(91,file=pdbfilei)
!  do str = 1,strnum(run,rep)
!    do n=1,nmontot
!     stpi => stp(str,n)
!     do while(associated(stpi%next))
!      write(91,100)'ATOM  ',stpi%atomnum,atomname(stpi%atomnum),sequence(n),' ',n,stpi%pos
!      stpi => stpi%next
!     enddo
!    enddo
!    write(91,'(a3)')'END'
!  enddo
!  close(91)

end subroutine unitcell

!............................................................

subroutine coarsegrain
 use parameters
 implicit none
 integer :: n, str, sgcount, k, m
 integer, dimension(5) :: lpcount
 real :: z0i

 posCa = 0.0; posN = 0.0; posH = 0.0; posC = 0.0; posO = 0.0; poss = 0.0
 posLP = 0.0; z0 = 0.0
!...assigning backbone CA and side chain centers of masses.....

  do str = 1,strnum(run,rep)
    z0i = 0.0
    nterm1_com = 0.0
    nterm2_com = 0.0
    cterm1_com = 0.0
    cterm2_com = 0.0
    mid1_com = 0.0
    mid2_com = 0.0
    nterm1_cnt = 0
    nterm2_cnt = 0
    cterm1_cnt = 0
    cterm2_cnt = 0
    mid1_cnt = 0
    mid2_cnt = 0
    do n=1,nmontot
       if(peptide(n))then
         m = n - nlip
         stpi => stp(str,n); sgcount = 0
         do while(associated(stpi%next))

            if( atomname(stpi%atomnum) == ' CA ') posCa(str,m,:) = stpi%pos
            if( atomname(stpi%atomnum) == ' N  ') posN(str,m,:) = stpi%pos
            if( atomname(stpi%atomnum) == ' HN ') posH(str,m,:) = stpi%pos
            if( atomname(stpi%atomnum) == ' C  ') posC(str,m,:) = stpi%pos
            if( atomname(stpi%atomnum) == ' O  ') posO(str,m,:) = stpi%pos

            if(hydro(stpi%atomnum)) then
             stpi => stpi%next; cycle
            endif
 
            if (m > nmon) then
                aanum = m - nmon
                if (aanum <= NTERM_END) then
                    nterm2_com = nterm2_com + stpi%pos(3);
                    nterm2_cnt = nterm2_cnt + 1
!                    print*,m, aanum, "nterm2", stpi%pos(3), nterm2_com, nterm2_cnt, nterm2_com/nterm2_cnt
                    if (aanum >= MID_START) then
                        mid2_com = mid2_com + stpi%pos(3);
                        mid2_cnt = mid2_cnt + 1
                    endif
                else
                    cterm2_com = cterm2_com + stpi%pos(3);
                    cterm2_cnt = cterm2_cnt + 1
!                    print*,m, aanum, "cterm2", stpi%pos(3), cterm2_com, cterm2_cnt, cterm2_com/cterm2_cnt
                endif
            else
                aanum = m
                if (aanum <= NTERM_END) then
                    nterm1_com = nterm1_com + stpi%pos(3);
                    nterm1_cnt = nterm1_cnt + 1
!                    print*,m, aanum, "nterm1", stpi%pos(3), nterm1_com, nterm1_cnt, nterm1_com/nterm1_cnt
                    if (aanum >= MID_START) then
                        mid1_com = mid1_com + stpi%pos(3);
                        mid1_cnt = mid1_cnt + 1
                    endif
                else
                    cterm1_com = cterm1_com + stpi%pos(3);
                    cterm1_cnt = cterm1_cnt + 1
!                    print*,m, aanum, "cterm1", stpi%pos(3), cterm1_com, cterm1_cnt, cterm1_com/cterm1_cnt
                endif
            endif

            if( (hydro(stpi%atomnum)).or.(backbone(stpi%atomnum)).or.&
                (terminal(stpi%atomnum)) )then
             stpi => stpi%next; cycle
            endif
            if(hydro(stpi%atomnum))then 
               print*,'Error: H considered for poss'; stop
            endif
            poss(str,m,:) = poss(str,m,:) + stpi%pos(:); sgcount = sgcount + 1
            stpi => stpi%next
         enddo
         term_com(str,1) = nterm1_com/nterm1_cnt
         term_com(str,2) = cterm1_com/cterm1_cnt
         term_com(str,3) = nterm2_com/nterm2_cnt
         term_com(str,4) = cterm2_com/cterm2_cnt
         term_com(str,5) = mid1_com/mid1_cnt
         term_com(str,6) = mid2_com/mid2_cnt
!         print*,'str ', str, nterm1_com, cterm1_com, nterm2_com, cterm2_com, nterm1_cnt,nterm2_cnt,cterm1_cnt,cterm2_cnt
         if(sgcount > 0)then
           poss(str,m,:) = poss(str,m,:)/sgcount
         else
           poss(str,m,:) = posCa(str,m,:)
         endif
       elseif(lipid(n))then
         m = n
         stpi => stp(str,n); lpcount = 0
         do while(associated(stpi%next))
            if(hydro(stpi%atomnum))then
               stpi => stpi%next; cycle
            endif
            if(choline(stpi%atomnum))then 
              posLP(str,m,1,:) = posLP(str,m,1,:) + stpi%pos
              lpcount(1) = lpcount(1) + 1
            elseif(phg(stpi%atomnum))then 
              if(atomname(stpi%atomnum) == ' P  ')then 
                 posLP(str,m,2,:) = stpi%pos
                 lpcount(2) = lpcount(2) + 1
                 z0i = z0i + posLP(str,m,2,3)
              endif
            elseif(glycerol(stpi%atomnum))then 
              posLP(str,m,3,:) = posLP(str,m,3,:) + stpi%pos
              lpcount(3) = lpcount(3) + 1
            elseif(fa1(stpi%atomnum))then 
              posLP(str,m,4,:) = posLP(str,m,4,:) + stpi%pos
              lpcount(4) = lpcount(4) + 1
            elseif(fa2(stpi%atomnum))then 
              posLP(str,m,5,:) = posLP(str,m,5,:) + stpi%pos
              lpcount(5) = lpcount(5) + 1
            else
              print*,'Error: unassigned ligand atom'; stop
            endif  
            stpi => stpi%next
         enddo
         do k=1,5
          posLP(str,m,k,:) = posLP(str,m,k,:)/lpcount(k)
         enddo
!         if( (lpcount(1)/=7).or.(lpcount(2)/=1).or.(lpcount(3)/=6).or. &
!         if( (lpcount(2)/=1).or.(lpcount(3)/=6).or. &
!             (lpcount(4)/=15).or.(lpcount(5)/=15) )then 
!!           print*,'Error: incorrect lpcount',lpcount; stop
!           print*,'Error: incorrect lpcount',lpcount(1), lpcount(2), lpcount(3), lpcount(4), lpcount(5)
!         endif
       elseif(water(n))then
         m = n-npmon-nlip-nNa
         stpi => stp(str,n)
         do while(associated(stpi%next))
            if( atomname(stpi%atomnum) == ' OH2') posW(str,m,:) = stpi%pos
            stpi => stpi%next
         enddo
       elseif(ionNa(n))then
         m = n-npmon-nlip
         stpi => stp(str,n)
         do while(associated(stpi%next))
            if( atomname(stpi%atomnum) == ' SOD') posNai(str,m,:) = stpi%pos
            stpi => stpi%next
         enddo
       elseif(ionCl(n))then
         m = n-npmon-nlip-nwater-nNa
         stpi => stp(str,n)
         do while(associated(stpi%next))
            if( atomname(stpi%atomnum) == ' CLA') posCli(str,m,:) = stpi%pos
            stpi => stpi%next
         enddo
       else
         print*,'Error: unassigned residue'; stop
       endif

    enddo
    z0(str) = z0i/nlip
  enddo

!  do str = 1,strnum(run,rep)
!     do m=1,nlig
!       do k=1,3
!          print*,str,m,posl(str,m,k,:)
!       enddo
!     enddo
!  enddo


end subroutine coarsegrain

!...........................................................................

subroutine structure
 use parameters
 implicit none
 integer, dimension(npmon,npmon) :: cm
 integer, dimension(npmon,nlip,5) :: cmlp
 integer, dimension(nmol,nmol) :: cntHH, cntEL, cntall
 integer, dimension(npmon) :: contlist
 integer, dimension(npmon,5) :: contlistlp, contlistPG, contlistPC
 integer, dimension(nlip) :: llist
 integer, dimension(nmol+nlip+nwater+nNa+nCl) :: ina
 integer, dimension(2,-1:1,-1:1,-1:1) :: crosslink
 real :: r, d, rmin, rb1
 real, dimension(nmol) :: rg, r1N
 real, dimension(3,nmol+nlip+nwater+nNa+nCl) :: cmpos
 real, dimension(npmon) :: rb, zsc
 real, dimension(3) :: posli, possi, shift, shiftmin, cmposi
 integer :: i,j, ic, n, str, q, k, l, m, i1, j1, j2,inorm, stri, ix,iy,iz, ilim, icbl,is
 type (ProteinStructure), pointer :: stpi1,stpi2
 logical :: CrossImage

  do str = 1,strnum(run,rep)
!...no coarse graining here .............
!..calculating centers of mass (geometrical) for individual peptides...
   cmpos = 0.0; ina=0
   do n=1,nmontot

    stpi => stp(str,n)
    do while(associated(stpi%next))
      if( .not.hydro(stpi%atomnum) ) then
                   cmpos(:,chain(n)) = cmpos(:,chain(n)) + stpi%pos
                   ina(chain(n)) = ina(chain(n)) + 1
      endif
      stpi => stpi%next
    enddo

   enddo
   do i=1,nmol+nlip+nwater+nNa+nCl
     cmpos(:,i) = cmpos(:,i) / ina(i)
   enddo
   if(sum(ina) /= na-nha)then
     print*,'Error: sum(ina) /= na-nha ',sum(ina),na-nha; stop
   endif

!...constructing contact map and Rg ................

   cm = 0; cmlp=0; rg = 0.0
   ic = 0
   do i=1,npmon-1
     do j=i+1,npmon
         r = sqrt( sum( (poss(str,i,:)-poss(str,j,:))**2 ) )
         if(chain(i+nlip)==chain(j+nlip))then 
           rg(chain(i)) = rg(chain(i)) + r*r; ic = ic + 1
         endif
         if( (j-i < 2).and.(chain(i+nlip)==chain(j+nlip)) )cycle
         if(r <= rc_cm) cm(i,j) = 1
     enddo
   enddo      
   
!...collecting Zcom for side chains
   do i = 1,npmon
     if(poss(str,i,3)>0)then
	   zsc(i) = poss(str,i,3) - z0(str)
	 elseif(poss(str,i,3)<0)then
	   zsc(i) = z0(str) - poss(str,i,3)
	 endif
   enddo
   

!...checking peptide self - interaction..............
   CrossImage=.false.
   do i=1,npmon
     do j=1,npmon
        if(chain(i+nlip)/=chain(j+nlip))cycle
        do ix = -1,1
        do iy = -1,1
        do iz = -1,1
          if( abs(ix)+abs(iy)+abs(iz)==0)cycle
          shift(1)=ix*box(str,1); shift(2)=iy*box(str,2); shift(3)=iz*box(str,3)
          possi =  poss(str,j,:) + shift
          r = sqrt( sum( (possi-poss(str,i,:))**2 ) )
          if(r <= rc_cm) CrossImage=.true.
        enddo
        enddo
        enddo
     enddo
   enddo
   if(CrossImage)ipim=ipim+1

!....protein-ligand interactions ...

     do i=1,npmon
      do j=1,nlip
       do ix = -1,1
       do iy = -1,1
       do iz = -1,1
         shift(1)=ix*box(str,1); shift(2)=iy*box(str,2); shift(3)=iz*box(str,3)
         do k=1,5
          r = sqrt( sum( (poss(str,i,:)+shift-posLP(str,j,k,:))**2 ) )
          if(r < rc_cm) cmlp(i,j,k) = 1
         enddo
       enddo
       enddo
       enddo
      enddo
     enddo

!..checking for lipids cross-bridging the peptides....
     ilim=0
     do j=1,nlip
       crosslink=0
       do i=1,npmon
        do ix = -1,1
        do iy = -1,1
        do iz = -1,1
          shift(1)=ix*box(str,1); shift(2)=iy*box(str,2); shift(3)=iz*box(str,3)
          do k=1,5
            r = sqrt( sum( (poss(str,i,:)+shift-posLP(str,j,k,:))**2 ) )
            if(r < rc_cm) crosslink(chain(i),ix,iy,iz) = 1
          enddo
        enddo
        enddo
        enddo
       enddo
       if( (sum(crosslink(1,:,:,:))>1).or.(sum(crosslink(2,:,:,:))>1) )ilim=ilim+1
     enddo

!......................................

     do i=1,npmon-1
      do j=i+1,npmon
        if(chain(i+nlip)==chain(j+nlip))then 
           rg(chain(i+nlip)) = rg(chain(i+nlip)) + sum((posCa(str,i,:)-posCa(str,j,:))**2)
           ic = ic + 1
        endif
      enddo
     enddo

     rg = 2.d0*rg
     ic = 2*ic
     do i=1,npmon
      do j=1,npmon
        if(i == j) cycle
        if(chain(i+nlip)==chain(j+nlip))then
          rg(chain(i+nlip)) = rg(chain(i+nlip)) + sum( (poss(str,i,:)-posCa(str,j,:))**2)
          ic = ic + 1
        endif
      enddo
     enddo
     rg = sqrt( rg / real( itot_rg ) )
     
!
!  FIX
!     if(ic /=  itot_rg*nmol ) then
!       print*,'Error: ic /=  itot_rg', ic,itot_rg; stop
!     endif
!
!  END FIX

!.....computing binding distances for aa....
     rb=0.0
     do i = 1, npmon
      rmin=1000000.0
      do ix = -1,1
      do iy = -1,1
      do iz = -1,1
        shift(1)=ix*box(str,1); shift(2)=iy*box(str,2); shift(3)=iz*box(str,3)
        do j=1,nlip
          do k=1,5
           r = sqrt( sum( (poss(str,i,:)-posLP(str,j,k,:)-shift)**2 ) )
           if(r<rmin) rmin=r
          enddo
        enddo
      enddo
      enddo
      enddo
      rb(i-nla) = rmin
     enddo
 if (not_now == 1) then

   contlist = 0; cntHH = 0; cntEL = 0; cntall = 0
   do i=1,npmon-1
     do j=i+1,npmon
        chain1 = chain(i+nlip)
        chain2 = chain(j+nlip)
        if(chain1 > chain2)then
           print*,'Error: chain(i+nlip)>chain(j+nlip)', i, j, chain(i+nlip), chain(j+nlip); stop
        endif
        if (chain1 > 2) print*, "CHAIN1 ", chain1
        if (chain2 > 2) print*, "CHAIN2 ", chain2
        if(cm(i,j)>1)then 
          print*,'Error: cm > 1'; stop
        endif
        if( hp(i)*hp(j)==1 ) then 
          if(chain1==chain2)then 
           contlist(i) = contlist(i) + cm(i,j)
           contlist(j) = contlist(j) + cm(i,j)
          endif
          cntHH(chain1,chain2) = cntHH(chain1,chain2) + cm(i,j)
        endif
        if( charge(i)*charge(j) == -1 ) &
            cntEL(chain1,chain2) = cntEL(chain1,chain2) + cm(i,j)
        if (cm(i,j) == 1) then
            cntall(chain1,chain2) = cntall(chain1,chain2) +1
        endif
     enddo
   enddo
   if(cntall(1,2) > 0) ippi = ippi + 1

   contlistlp = 0; llist=0
   contlistPC = 0;
   contlistPG = 0;
   do i=1,npmon
     do j=1,nlip 
        do k=1,5
         if( cmlp(i,j,k)>0 ) then 
           contlistlp(i,k) = contlistlp(i,k) + 1
           llist(j)=1
           if  (isPC(j) .eqv. .true.) then
               contlistPC(i,k) = contlistPC(i,k) + 1
           else
               contlistPG(i,k) = contlistPG(i,k) + 1
           endif
         endif
        enddo
     enddo
   enddo
   if(sum(llist)>0)plim=plim + ilim/real(sum(llist))
   icbl=0
   do j=1,nlip
    if((sum(cmlp(1:nmon,j,:)) > 0).and.(sum(cmlp(1+nmon:npmon,j,:)) > 0))icbl=icbl+1
   enddo
   if(sum(llist)>0)fcbl = fcbl + icbl/real(sum(llist))


   if( sum(cntall) /= sum(cm(1:npmon,1:npmon)) ) then
     print*,'Error: cntall and cm do not match',cntall,sum(cm(1:npmon,1:npmon));stop
   endif

!   distKD(1) = sqrt( sum( (poss(str,14,:)-poss(str,19,:))**2 ) )
!   distKD(2) = sqrt( sum( (poss(str,45,:)-poss(str,50,:))**2 ) )

   do m = 1, nmol
    r1N(m) = sqrt( sum((posCa(str,m*nmon,:)-posCa(str,(m-1)*nmon+1,:))**2) )
    if(r1N(m) > 100.)then 
      print*,'Error: r1N is too large',str,r1N; stop
    endif
   enddo
 endif

   stri = sum(strnum(runmin:run-1,rep)) + str
   write(24,'(i6,6(1x,i3))')stri,(cntHH(i,i),i=1,2),(cntEL(i,i),i=1,2),&
                                 (cntall(i,i),i=1,2)
   write(25,'(i6,4(1x,f8.4))')stri,rg, r1N
   write(27,'(i6,62(i2))')stri,contlist
   write(28,'(i6,1922(i1))')stri,((cm(i,j),j=i+1,nmon),i=1,nmon-1),&
                                 ((cm(i,j),j=i+1,npmon),i=nmon+1,npmon-1)

   write(43,'(i6,210(i2))')stri,(contlistlp(i,:),i=1,npmon)
   write(52,'(i6,210(i2))')stri,(contlistPC(i,:),i=1,npmon)
   write(53,'(i6,210(i2))')stri,(contlistPG(i,:),i=1,npmon)
   write(45,'(i6,98(i2))')stri,llist
   write(49,'(i6,62(f6.2,1x))')stri,rb
   write(50,'(i6,5704(i1))')stri,((sum(cmlp(i,j,:)),j=1,nlip),i=1,npmon)
   write(51,'(i6,22(1x,f7.3))')stri,zsc(:)
   write(54,'(i6,6(1x,f7.3))')stri,term_com(str,1:6)-z0(str)
   
 enddo 

end subroutine structure


!.................................................................

subroutine WTDihedral
 use parameters
 implicit none
 integer :: i,j, m, k, str, stri
 real, dimension(npmon-2*nmol) :: dh_phi, dh_ksi
 real, dimension(4,3) :: posi

 interface 
  function  DihedralAngle(pos)
    implicit none
    real, dimension(4,3), intent(in) :: pos
    real :: DihedralAngle
  end function  DihedralAngle
 end interface 

 do str = 1,strnum(run,rep)
 
  dh_phi = 0.0; dh_ksi = 0.0;  k=0
  do m=1,nmol
   do i=(m-1)*nmon+2,m*nmon-1
     k = k + 1

     posi(1,:) = posC(str,i-1,:); posi(2,:) = posN(str,i,:)
     posi(3,:) = posCa(str,i,:);  posi(4,:) = posC(str,i,:)
     dh_phi(k) = DihedralAngle(posi)
     dh_phi(k) = 180.*dh_phi(k)/pi

     posi(1,:) = posN(str,i,:); posi(2,:) = posCa(str,i,:)
     posi(3,:) = posC(str,i,:); posi(4,:) = posN(str,i+1,:)
     dh_ksi(k) = DihedralAngle(posi)
     dh_ksi(k) = 180.*dh_ksi(k)/pi

   enddo
  enddo

  stri = sum(strnum(runmin:run-1,rep)) + str
  write(31,*)stri,dh_phi,dh_ksi

 enddo


end subroutine WTDihedral


!.........................................................................

function DihedralAngle(pos)
 implicit none
 real :: x21,y21,z21,x23,y23,z23,x34,y34,z34,r34dotr23,r21dotr23,r23sq,&
        onex,oney,onez,twox,twoy,twoz,rone,rtwo,rprod,cosphi,arg,sinphi,sign
 real, dimension(4,3), intent(in) :: pos
 real :: DihedralAngle

     x21 = pos(1,1)-pos(2,1)
     y21 = pos(1,2)-pos(2,2)
     z21 = pos(1,3)-pos(2,3)

     x23 = pos(3,1)-pos(2,1)
     y23 = pos(3,2)-pos(2,2)
     z23 = pos(3,3)-pos(2,3)

     x34 = pos(4,1)-pos(3,1)
     y34 = pos(4,2)-pos(3,2)
     z34 = pos(4,3)-pos(3,3) 

     r34dotr23 = x34*x23 + y34*y23 + z34*z23
     r21dotr23 = x21*x23 + y21*y23 + z21*z23
     r23sq = x23**2 + y23**2 + z23**2

     onex = x34 - r34dotr23 * x23 / r23sq
     oney = y34 - r34dotr23 * y23 / r23sq
     onez = z34 - r34dotr23 * z23 / r23sq

     twox = x21 - r21dotr23 * x23 / r23sq 
     twoy = y21 - r21dotr23 * y23 / r23sq 
     twoz = z21 - r21dotr23 * z23 / r23sq 

     rone = sqrt(onex**2 + oney**2 + onez**2)
     rtwo = sqrt(twox**2 + twoy**2 + twoz**2)
     rprod = rone * rtwo
     cosphi = (onex*twox + oney*twoy + onez*twoz) / rprod
     if(cosphi > 1.0) cosphi = 1.0
     if(cosphi < -1.0) cosphi = -1.0
     arg = 1.0 - cosphi**2

     sinphi = sqrt(arg)
!..getting the sign of the sin(phi) ...
     sign = (oney*twoz - twoy*onez)*x23 + (twox*onez - onex*twoz)*y23 + &
            (onex*twoy - twox*oney)*z23 
     if(sign > 0.0) sinphi = -sinphi
     DihedralAngle = acos(cosphi)
     if(sinphi < 0.0) DihedralAngle = - DihedralAngle

end function DihedralAngle

!.................................................................

subroutine Hbonds
 use parameters
 implicit none
 real, parameter :: qe1=0.42, qe2=0.20
 integer :: i,j,k,str,l,m,n,stri,dot,d
 integer, dimension(npmon,npmon) :: hbmap 
 integer, dimension(npmon*2) :: hblist,hbindex
 integer, dimension(npmon,3) :: hbres
 integer, dimension(nmol) :: hb
 integer, dimension(2,4) :: si
 real, dimension(3) :: u1,u2
 real, dimension(npmon,npmon) :: hbmapE
 real :: rHN, rHO, rNO, rNC, rHC, alphaOHN, Ehb, EhbNH, EhbCO, rmin, &
         Score,ScoreNH,ScoreCO

 si(1,:)=(/-4,-3,4,5/); si(2,:)=(/-5,-4,3,4/)
 do str=1, strnum(run,rep)

  hb = 0
  hbmap = 0; hblist = 0; hbindex=0; hbmapE=0.0;  hbres=0
  do i=1,npmon
   do j=1,npmon
     if( (chain(i+nlip)/=chain(j+nlip)).or.(abs(i-j)<3) ) cycle
     if(chainterminal(i).or.chainterminal(j))cycle

     rHO = sqrt( sum( (posO(str,j,:) - posH(str,i,:))**2 ) )
     rHN = sqrt( sum( (posN(str,i,:) - posH(str,i,:))**2 ) )
     rNO = sqrt( sum( (posN(str,i,:) - posO(str,j,:))**2 ) )
     rNC = sqrt( sum( (posN(str,i,:) - posC(str,j,:))**2 ) )
     rHC = sqrt( sum( (posH(str,i,:) - posC(str,j,:))**2 ) )

     u1(:) = (posN(str,i,:)-posH(str,i,:))/rHN
     u2(:) = (posO(str,j,:)-posH(str,i,:))/rHO
     if(dot_product(u1,u2) > 1.0+eps)then
       print*,'Error: cos > 1'; stop
     endif
     alphaOHN = acos( dot_product(u1,u2) )
     Score = (alphaOHN/pi)**2 + (rDAmin/rNO)**2
     Ehb = qe1*qe2*332.*(1./rNO + 1./rHC - 1./rHO - 1./rNC)

     if(HBdefinition=='geometric')then
         if( (rNO <= rc_NO).and.(alpha <= alphaOHN) ) then 
          if( (sum(hbmap(i,:)) == 0).and.(sum(hbmap(:,j)) < 2) )then
             hbmap(i,j)=1; hbmapE(i,j)=Score
          elseif( (sum(hbmap(i,:)) == 1).and.(sum(hbmap(:,j)) < 2) )then
            if(Score > maxval(hbmapE(i,:)) ) then 
                hbmap(i,:)=0; hbmapE(i,:)=0.0
                hbmap(i,j)=1; hbmapE(i,j)=Score
            endif
          elseif( (sum(hbmap(i,:)) == 0).and.(sum(hbmap(:,j)) ==2 ) )then
            if(Score > minval(hbmapE(:,j),mask=hbmapE(:,j) > eps) ) then
                k = minval(minloc(hbmapE(:,j),mask=hbmapE(:,j) > eps))
                hbmap(k,j)=0; hbmapE(k,j)=0.0 
                hbmap(i,j)=1; hbmapE(i,j)=Score
            endif
          elseif( (sum(hbmap(i,:)) == 1).and.(sum(hbmap(:,j)) ==2 ) )then
            ScoreNH = maxval(hbmapE(i,:))
            ScoreCO = minval(hbmapE(:,j),mask=hbmapE(:,j) > eps)
            if(Score > ScoreNH + ScoreCO)then 
                hbmap(i,:)=0; hbmapE(i,:)=0.0
                k = minval(minloc(hbmapE(:,j),mask=hbmapE(:,j) > eps))
                hbmap(k,j)=0; hbmapE(k,j)=0.0
                hbmap(i,j)=1; hbmapE(i,j)=Score
            endif
          else
             print*,'Unusual HB geometry'; stop
          endif
         endif
     elseif(HBdefinition=='energetic')then
         if( Ehb <= Ehb0)  then
          if( (sum(hbmap(i,:)) == 0).and.(sum(hbmap(:,j)) < 2) )then
             hbmap(i,j)=1; hbmapE(i,j)=Ehb
          elseif( (sum(hbmap(i,:)) == 1).and.(sum(hbmap(:,j)) < 2) )then
            if(Ehb < minval(hbmapE(i,:)) ) then 
                hbmap(i,:)=0; hbmapE(i,:)=0.0
                hbmap(i,j)=1; hbmapE(i,j)=Ehb 
            endif
          elseif( (sum(hbmap(i,:)) == 0).and.(sum(hbmap(:,j)) ==2 ) )then
            if(Ehb < maxval(hbmapE(:,j),mask=hbmapE(:,j) < 0.0) ) then
                k = minval(maxloc(hbmapE(:,j),mask=hbmapE(:,j) < 0.0))
                hbmap(k,j)=0; hbmapE(k,j)=0.0 
                hbmap(i,j)=1; hbmapE(i,j)=Ehb
            endif
          elseif( (sum(hbmap(i,:)) == 1).and.(sum(hbmap(:,j)) ==2 ) )then
            EhbNH = minval(hbmapE(i,:))
            EhbCO = maxval(hbmapE(:,j),mask=hbmapE(:,j) < 0.0)
            if(Ehb < EhbNH + EhbCO)then 
                hbmap(i,:)=0; hbmapE(i,:)=0.0
                k = minval(maxloc(hbmapE(:,j),mask=hbmapE(:,j) < 0.0))
                hbmap(k,j)=0; hbmapE(k,j)=0.0
                hbmap(i,j)=1; hbmapE(i,j)=Ehb
            endif
          else
             print*,'Unusual HB geometry'; stop
          endif
         endif
     else; print*,'Invalid HBdefinition'; stop
     endif

   enddo
  enddo


!...testing if bifurcated HBs are there...........
  do i=1,npmon
   if(sum(hbmap(i,:)) > 1)then
       print*,'Error: bifurcated bond:',i; stop
   endif
  enddo
  do j=1,npmon
   if(sum(hbmap(:,j)) > 2)then
       print*,'Error: bifurcated bond:',j; stop
   endif
  enddo

!.....transferring hbmap onto hblist arrays...............

  do i=1,npmon
   do j=1,npmon
     if( (chain(i+nlip)/=chain(j+nlip)).or.(abs(i-j)<3) ) cycle
     if(hbmap(i,j)==0)cycle
     if(chain(i+nlip)==chain(j+nlip))hb(chain(i+nlip)) = hb(chain(i+nlip)) + hbmap(i,j)
     hbres(i,1) = j
     hbindex(2*i-1) = 1
     if(hbindex(2*j)==0)then 
        hbindex(2*j) = 1
        hbres(j,2) = i
     else 
        hbres(j,3) = i
     endif 
     hblist(2*i-1) = hblist(2*i-1) + hbmap(i,j)
     hblist(2*j) = hblist(2*j) + hbmap(i,j)

   enddo
  enddo

!...this test is designed to see if a hydrogen makes more than one bond...
  do i=1,2*npmon-1,2
    if( hblist(i) > 1 )then
      print*,'Error: hblist(:) > 1',str,i,hblist(i); stop
    endif
  enddo

!...this test is designed to see if an oxygen makes more than two bonds ..
  do i=2,2*npmon,2
    if( hblist(i) > 2 )then
      print*,'Error: hblist(:) > 2',str,i,hblist(i); stop
    endif
  enddo


  stri = sum(strnum(runmin:run-1,rep)) + str
  write(35,'(i6,2(1x,i2))')stri,hb
  write(36,'(i6,124(1x,i1))')stri,hblist
  write(38,'(i6,186(1x,i2))')stri,((hbres(i,j),j=1,3),i=1,npmon)

enddo

end subroutine Hbonds
