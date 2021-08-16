!! GLWH Dec 2020
!! <summary>This program reads in a POSCAR file and finds the 
!! number of operations in the spacegroup.</summary>
PROGRAM symdriver
use num_types

use symmetry
implicit none

real(dp) LV(3,3)
integer nAt, nPG, i, nTypes, status, iTypes, idx
character(len=128) filename, dummy, line1, line2
real(dp), allocatable    :: d(:,:)
real(dp), allocatable    :: rots(:,:,:), shifts(:,:)
integer, allocatable     :: aTyp(:)

call getarg(1,filename)
if (filename=="") then; filename = "POSCAR"; endif
 
open(33,file=filename,status="old")
read(33,*) dummy !Skip title line
read(33,*) dummy !Skip lattice parameter
read(33,*) (LV(i,:),i=1,3)
! Count the number of atoms and types
read(33,'(A)') line1
line1 = trim(adjustl(line1))
line2 = line1
nTypes = 0; nAt = 0
do
  read(line1,*,iostat=status) i
  if (status/=0) exit
  nTypes = nTypes + 1
  nAt = nAt + i
  line1 = adjustl(line1(index(line1," "):))
  print*,line1
enddo
write(*,'(i3," types found")') nTypes
write(*,'(i3," atoms found")') nAt
allocate(d(nAt,3),aTyp(nAt))
! Read in the number of atoms *of each type*
idx = 1; aTyp = 0
do i = 1, nTypes
    read(line2,*,iostat=status) iTypes
    if (status/=0) exit
    aTyp(idx:idx+iTypes-1) = i ! should be i 
    idx = idx + iTypes
    line2 = adjustl(line2(index(line2," "):))
enddo
read(33,*) !skip cart/direct line
do i = 1, nAt 
    read(33,*) d(:,i)
    write(*,'(3(f8.4,1x))') d(:,i)
enddo
call get_spacegroup(LV, aTyp, d, rots, shifts, .true.)
close(33)

nPG=size(rots,3)
print*,"shape(rots)",shape(rots)
print*,"Number of rotations:",nPG
!print*,rots(:,:,1)
!print*,rots(:,:,2)


END PROGRAM symdriver