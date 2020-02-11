! stand alone driver to read in lattice vectors and atomic basis
! and return the spacegroup rotations and shift.
PROGRAM spacegroup
use num_types
use symmetry
implicit none

real(dp), dimension(3,3) :: LV
real(dp), allocatable    :: d(:,:)
real(dp), allocatable    :: rots(:,:,:), shifts(:,:)
real(dp)                 :: dummy
integer, allocatable     :: aTyp(:)
character(120)           :: LVfile, dFile
integer iD, nD, iLV, i, iRot, j
integer nArg, LVinID, dINid, status
character(10000) line
LVinID = 11
dINid = 12

nArg = command_argument_count()
if (nArg /= 2)  &
stop "ERROR: requires two argments. File names for lattice vectors and for atomic basis vectors"
call get_command_argument(1,LVfile)
open(LVinID,file=LVfile,status='old')
call get_command_argument(2,dFile)
open(dINid,file=dFile,status='old')

! Read the two file names for LV and d from the command line
! Skip the first line, which is a comment for fortpy
read(LVinID,*)
do iLV = 1, 3
  read(LVinID,*) LV(iLV,:)
  write(*,'("LV: ",3(f7.3,","))') LV(iLV,:)
enddo
print*

! Count the number of atoms, allocate atom positions and types
! first line is a comment, so skip that
read(dINid,*)
read(dINid,'(A)') line
line = trim(adjustl(line))
backspace(dINid)
nD = 0
do
  read(line,*,iostat=status) dummy
  if (status/=0) exit
  nD = nD + 1
  line = adjustl(line(index(line," "):))
enddo
write(*,'(i3," atoms found")') nD

allocate(d(3,nD),aTyp(nD))
aTyp = 1 ! Set all atoms to have the same type
!Read in the atomic positions, rows are x,y,z, columns are d-vectors
do iD = 1, 3
  read(dINid,*) d(iD,:)
enddo

! Call the spacegroup routine
call get_spacegroup(LV, aTyp, d, rots, shifts, .false.)
write(*,'(i3," spacegroup operations")') size(rots,3)
close(dINid)
close(LVinID)

! Write out the rotations and shifts
open(11,file="rots.out")
write(11,"('# <fortpy version=""1"" template=""real""></fortpy>')")
write(11,'("##",3(1x,i13))') shape(rots)
do i = 1, size(rots,1)
  write(11,'("##",3(1x,i13))') i,0,0
  do j = 1, size(rots,2)
    write(11,'(300(f19.16,4x))') rots(i,j,:)
  enddo
enddo
close(11)

open(11,file="shifts.out")
write(11,"('# <fortpy version=""1"" template=""real""></fortpy>')")
do i = 1, size(shifts,1)
  write(11,'(300(f19.16,4x))') shifts(i,:)
enddo
close(11)

open(11,file="nD.out")
write(11,"('# <fortpy version=""1"" template=""integer""></fortpy>')")
write(11,'(i10)') nD
close(11)

END PROGRAM spacegroup
