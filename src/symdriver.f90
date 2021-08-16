!! GLWH Dec 2020
!! <summary>This program reads in a POSCAR-like file and finds the 
!! number of operations in the spacegroup.</summary>
!! The input file should be the three lattice vectors (one in each row)
!! then the number of atoms on the next line
!! then the atomic type (integer) followed by three floats (for the position)
PROGRAM symdriver
use num_types

use symmetry
implicit none

real(dp) LV(3,3)
integer nAt, nSG, nLPG, i
character(len=128) filename
real(dp), allocatable    :: d(:,:)
real(dp), allocatable    :: rots(:,:,:), shifts(:,:)
real(dp), pointer        :: pg(:,:,:)
integer, allocatable     :: aTyp(:)

call getarg(1,filename)
if (filename=="") then; filename = "symcheck.str"; endif
 
open(33,file=filename,status="old")
read(33,*) (LV(:,i),i=1,3)
! Count the number of atoms and types
read(33,*) nAt
allocate(aTyp(nAt),d(3,nAt))
do i = 1,nAt
  read(33,*) aTyp(i), d(:,i)
enddo
call get_lattice_pointGroup(LV, pg)
call get_spacegroup(LV, aTyp, d, rots, shifts, .false.)
close(33)

nSG=size(rots,3); nLPG=size(pg,3)
write(*,'(i3,2x,"(Number of spacegroup operations)")') nSG
write(*,'(i3,2x,"(Number of lattice PG operations)")') nLPG
END PROGRAM symdriver