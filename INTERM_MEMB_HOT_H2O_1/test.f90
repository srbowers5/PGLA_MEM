
module parameters


integer, parameter :: nlip=4, npmon=2

end module parameters

subroutine structure
 use parameters
 integer, dimension(npmon,nlip,3) :: cmlp
 integer :: stri

cmlp=0
cmlp(1,1,1) = 1
cmlp(1,1,2) = 2
cmlp(1,1,3) = 3
cmlp(2,1,1) = 0
cmlp(2,1,2) = 0
cmlp(2,1,3) = 0

cmlp(1,2,1) = 1
cmlp(1,2,2) = 0
cmlp(1,2,3) = 3
cmlp(2,2,1) = 0
cmlp(2,2,2) = 0
cmlp(2,2,3) = 1


cmlp(1,3,1) = 1
cmlp(1,3,2) = 0
cmlp(1,3,3) = 2
cmlp(2,3,1) = 0
cmlp(2,3,2) = 0
cmlp(2,3,3) = 2

cmlp(1,4,1) = 5
cmlp(1,4,2) = 0
cmlp(1,4,3) = 1
cmlp(2,4,1) = 2
cmlp(2,4,2) = 2
cmlp(2,4,3) = 2
stri=543

write(50,'(i6,8(i1))')stri,((sum(cmlp(i,j,:)),j=1,nlip),i=1,npmon)

!6  = 1,1
!4  = 1,2
!3  = 1,3
!6  = 1,4
!0  = 2,1
!1  = 2,2
!2  = 2,3
!6  = 2,4


end subroutine structure


program interm
 use parameters
 open(50,file='test1.dat')
  call structure
  close(50)
end program interm



