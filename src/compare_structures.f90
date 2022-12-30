!!<summary>The point of this module take two structures and compare them to see if the are
!! equivalent. Implicit in the approach is the assumption that the structures are derivatives
!! structures of a parent lattice (the same one if they are equivalent).
!! Gus Hart Dec. 2006
!! * The code as it stands (12/2007) occasionally finds matches that really aren't
!!(or so it seems). Rather than debugging it, I'm going add new routines that compare structures
!!using the "enum" approach---the g-representation---avoiding as much as possible, any geometric
!!comparisons.</summary>
MODULE compare_structures
use num_types
use vector_matrix_utilities, only: volume, determinant, matrix_inverse
use numerical_utilities, only: equal
use symmetry

implicit none
private
public compare_arbitrary_structures, is_derivative, is_equiv_lattice, is_lattice_point

CONTAINS

  !!<summary>COMPARE ARBITRARY STRUCTURES
  !! This subroutine takes two structures, str1 and str2, and compares them to
  !! see if they are equivalent.
  !!  - Rescales the structures to have the same volume/atom
  !!  - Checks that both structures reside on the same
  !! underlying lattice, and that the atomic sites of str2 also
  !! lie on the lattice of str1
  !!  - Uses a set of spacegroup operations to find equivalent,
  !! but mis-oriented structures</summary>
  !!<parameter name="LV1in" regular="true">Lattice vectors for the first structure.</parameter>
  !!<parameter name="LV2in" regular="true">Lattice vectors for the second stucture.</parameter>
  !!<parameter name="aTyp1in" regular="true">Atom type of the first structure.</parameter>
  !!<parameter name="aTyp2in" regular="true">Atom type of the second structure.</parameter>
  !!<parameter name="aPos1in" regular="true">Atomic positions for the first structure.
  !!</parameter>
  !!<parameter name="aPos2in" regular="true">Atomic positions for the second stucture.
  !!</parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  !!<parameter name="status" regular="true">Returns the similarity/differences in the structures:
  !! Status values/meanings:
  !! 0 -- structures are equivalent
  !! 1 -- structures are inequivalent (but the following are not true)
  !! 2 -- underlying lattices are different
  !! 3 -- atoms of str2 do not lie on the underlying lattice of str1
  !! 4 -- number of atoms in the compared structures are different
  !! 5 -- number of types of atoms are not the same in each structure</parameter>
  !!<parameter name="irot" regular="true">ith rotation that mapped structure one onto structure
  !!two.</parameter>
  !!<parameter name="mapped" regular="true">True if structure 1 and structure 2 are symmetrically
  !!equivalent.</parameter>
  !!<parameter name="SMskipout" regular="true">Number of skipped cells.</parameter>
  !!<parameter name="identical" regular="true">Logical that returns if the structures are the
  !!same under translations.</parameter>
  subroutine compare_arbitrary_structures(LV1in,aTyp1in,aPos1in,LV2in,aTyp2in,aPos2in,&
       eps,mapped,status,irot,SMskipout,identical)
    real(dp),intent(in):: LV1in(3,3), LV2in(3,3)
    real(dp),intent(in):: aPos1in(:,:), aPos2in(:,:)! Atomic positions for each structure
    integer,intent(in) :: aTyp1in(:), aTyp2in(:)    ! Atom types for each structure
    real(dp), intent(in):: eps ! Finite precision tolerance
    logical,intent(out):: mapped
    logical,optional,intent(out):: identical
    integer,optional,intent(out) :: status, irot, SMskipout

    ! Status values/meanings:
    ! 0 -- structures are equivalent
    ! 1 -- structures are inequivalent (but the following are not true)
    ! 2 -- underlying lattices are different
    ! 3 -- atoms of str2 do not lie on the underlying lattice of str1
    ! 4 -- number of atoms in the compared structures are different
    ! 5 -- number of types of atoms are not the same in each structure

    real(dp), allocatable:: rots(:,:,:), shifts(:,:)! Space group rotations/shifts for
    ! the underlying lattice
    integer :: N1, N2 ! Number of atoms in each structure
    real(dp)         :: L2C(3,3), C2L(3,3) ! latt2cart, cart2latt conversion matrices
    integer          :: j, jrot, ip  ! Loop counters
    real(dp)         :: rescale
    integer,pointer  :: pTyp(:)   ! The "p" stands for the "primitive" cell
    real(dp),allocatable :: pPos(:,:) ! of the underlying lattice
    real(dp)         :: pVecs(3,3)! of str1
    real(dp),allocatable:: rotPos2(:,:)
    integer,allocatable :: rlaTyp2(:) ! "relabeled" list of types in the basis of str2
    real(dp)         :: LV1(3,3), LV2(3,3), LV1inv(3,3)
    real(dp)         :: aPos1(3,size(aTyp1in)), aPos2(3,size(aTyp2in))
    integer          :: aTyp1(size(aTyp1in)), aTyp2(size(aTyp2in))
    real(dp)         :: vpa1, vpa2 ! Volume/atom for each structure
    real(dp), parameter:: ident(3,3)=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    logical            :: err
    integer,pointer    ::  block1(:), block2(:), perms(:,:),uqlist1(:),uqlist2(:)
    integer :: il, lab ! loop counter over label types in str1, temporary label
    integer :: SMskip  ! counts number of lattices skipped because they don't match in the
                       ! Santoro and Mighell sense.

    if(present(status)) status = 0
    N1 = size(aTyp1) ! Number of atoms in str1
    N2 = size(aTyp2) ! Number of atoms in str2
    mapped = .false. ! Flag indicating whether structures are equivalent
    SMskip = 0       ! Number of times a rotation was skipped because, for
                     ! this particular orientation, the cells of str1 and str2
                     ! aren't derivative lattices
    identical = .false.
    LV1 = LV1in; aPos1 = aPos1in; aTyp1 = aTyp1in
    LV2 = LV2in; aPos2 = aPos2in; aTyp2 = aTyp2in
    ! Check that structures have the same # atoms/cell
    if (.not.N1==N2) then
       if (present(status)) status=4
       return
    endif

    ! Get inverse of LV1 for use later
    call matrix_inverse(LV1,LV1inv,err)
    if (err) stop "Lattice vectors of first input structure are co-planar"

    ! Map all atoms in each structure into the unit cell (in case they're not
    ! input that way).
    call get_transformations(LV1,L2C,C2L)
    do j = 1,N1
       call bring_into_cell(aPos1(:,j),C2L,L2C,eps)
    enddo

    call get_transformations(LV2in,L2C,C2L)
    do j = 1,N2
       call bring_into_cell(aPos2(:,j),C2L,L2C,eps)
    enddo
    ! Find the symmetry of the underlying lattice of str1
    !  - Make copy of str1, make all atoms equivalent, then find symmetry
    allocate(pTyp(N1),pPos(3,N1))
    pTyp = 1; pPos = aPos1; pVecs = LV1
    ! Get the underlying (multi)lattice
    call make_primitive(pVecs,pTyp,pPos,.false.,eps)
    call get_spacegroup(pVecs,pTyp,pPos,rots,shifts,.false.,eps)

    ! This block could be its own subprogram
    ! Make the volume/atom the same for each structure
    vpa2 = abs(determinant(LV2in))  ! volume/atom for str2
    vpa1 = abs(determinant(LV1))    ! volume/atom for str1
    rescale = (vpa1/vpa2)**(1/3._dp)! scaling factor
    aPos2 = rescale*aPos2           ! rescale the atom positions (str2)
    LV2 = rescale*LV2in             ! rescale the lattice vectors (str2)

    ! These next two checks are redundant with main loop--early out
    ! Is lattice of str2 a derivative of str1's underlying lattice?
    if (.not. is_derivative(pVecs,LV2,eps)) then
       if (present(status)) status=2;
       deallocate(pTyp,pPos)
       return;
    endif

    ! Are the atoms of str2 on lattice sites of str1's underlying lattice?
    if (.not. are_lattice_points(pVecs,aPos2,pPos,eps)) then
       if (present(status)) status=3;
       deallocate(pTyp,pPos)
       return;
    endif

    call unique(aTyp1,uqlist1,block1)
    call unique(aTyp2,uqlist2,block2)
    if (size(uqlist1)/=size(uqlist2)) then
       if (present(status)) status = 5;
       deallocate(pTyp,pPos)
       return;
    endif

    ! Get a list of permutations here
    call get_basis_permutations(uqlist2,perms)
    ! Loop over all possible operations of the spacegroup
    ! (because the structures may be equivalent but not aligned)
    mapped = .false.; if(present(identical)) identical = .false.
    allocate(rotPos2(3,N1),rlaTyp2(N1))
    rotns_loop:do jrot=1,size(rots,3)
       ! Rotate the lattice using the jrot-th rotation of parent lattice pointgroup
       LV2=matmul(rots(:,:,jrot),rescale*LV2in)
       ! Check that str2 lattice is a derivative lattice.
       if (.not. is_equiv_lattice(LV1,LV2,eps)) then
          SMskip = SMskip + 1  ! Keep track for debugging/testing
          !write(*,*) "## Santoro-Mighell skip ##"
          cycle
       endif

       ! If we get to here then the lattice vectors are equivalent
       ! (derivative lattices). So now check atomic configuration.
       ! I think it makes sense that the rotated atomic basis vectors
       ! should be shifted by the non-symmorphic shift of the rotation,
       ! so add it to each of the atomic positions of cell 2 (the rotated one)
       rotPos2 = matmul(rots(:,:,jrot),aPos2) + spread(shifts(:,jrot),2,N2)
       do ip = 1,size(perms,1) ! Loop over permutations of atom types
          if (any(block2(perms(ip,:))/=block1)) then
             cycle
          endif
          ! We need to substitute the atom types in case a different labeling scheme
          ! was used on the two structures being compared
          do il = 1,size(uqlist1) ! Loop over labels in str1 and relabel str2
             lab = uqlist1(perms(ip,il)) ! current label under ip-th permutation
             where(aTyp1==lab); rlaTyp2 = lab; endwhere
          enddo
          ! The first checks for a mapping translation with str2's atoms relabeled
          if (mapping_translation_exists(LV2,aPos1,rotPos2,aTyp1,rlaTyp2, eps)) then
             mapped = .true.; endif
          ! This check looks for mapping translations *without* relabeling any atoms
          if (present(identical)) then
             if (mapping_translation_exists(LV2,aPos1,rotPos2,aTyp1,aTyp2, eps)) then
                identical = .true.; endif
          endif
          if(mapped .and. .not. present(identical)) exit rotns_loop
          if (present(identical)) then
             if (mapped .and. identical) exit rotns_loop
          endif
       enddo
    enddo rotns_loop
    ! Set some debugging/testing values, if present in argument list
    if (present(status)) status=0  ! Structures were found to be equivalent
    if (present(status) .and. .not. mapped) status=1 ! They are inequivalent
    if (present(irot)) irot=jrot ! current spacegroup op # when loop ended
    if (present(SMskipout)) SMskipout = SMskip

    deallocate(pTyp,pPos)
    deallocate(rotPos2,rlaTyp2)

  end subroutine compare_arbitrary_structures

  !!<summary>This function determines whether two lattices are equivalent. That is, if one
  !! is an equal-volume derivative lattice of the other. It uses the idea of
  !! Santoro and Mighell (Acta. Cryst. 1972) that, for equivalent lattices, the
  !! transformation matrix that takes the vectors of str2 to str1 has determinant
  !! of 1 and all integer elements. [LV2]=[LV1][S]=&gt;S=LV1^-1*LV2. But I had to use
  !! abs on the determinant. This is not mentioned by Santoro and Mighell (issue of
  !! opposite handedness---which doesn't matter to me).</summary>
  !!<parameter name="lat1" regular="true">The first lattice.</parameter>
  !!<parameter name="lat2" regular="true">The second lattice.</parameter>
  !!<parameter name="eps" regular="true">The finite precision tolerance.</parameter>
  function is_equiv_lattice(lat1,lat2,eps)
    real(dp), intent(in) :: lat1(3,3), lat2(3,3), eps
    logical :: is_equiv_lattice, err
    real(dp)  :: lat1inv(3,3), S(3,3)
    real(dp) :: atol  ! An absolute tolerance for the "equal" function

    atol = 5E-4_dp
    is_equiv_lattice = .false.
    call matrix_inverse(lat1,lat1inv,err,10d-13)
    if (err) stop "Problem with input vectors in function 'is_equiv_lattice'"
    S = matmul(lat1inv,lat2)
    if (equal(abs(determinant(S)),1._dp,eps,atol) .and. &
         equal(S,nint(S),eps,atol)) is_equiv_lattice = .true.
  endfunction is_equiv_lattice

  !!<summary>Same as "is_equiv_lattice" except without the volume check. That is, the
  !!derivative lattice may be a sublattice (larger volume) of the original.</summary>
  !!<parameter name="lat1" regular="true">The first lattice.</parameter>
  !!<parameter name="lat2" regular="true">The second lattice.</parameter>
  !!<parameter name="eps" regular="true">The finite precision tolerance.</parameter>
  function is_derivative(lat1,lat2,eps)
    real(dp), intent(in) :: lat1(3,3), lat2(3,3), eps

    logical :: is_derivative, err
    real(dp)  :: lat1inv(3,3), S(3,3)
    real(dp) :: atol  ! An absolute tolerance for the "equal" function
    atol = 5E-4_dp
    is_derivative = .false.
    call matrix_inverse(lat1,lat1inv,err)
    if (err) stop "Problem with input vectors in function 'is_derivative'"
    S = matmul(lat1inv,lat2)
    if (equal(S,nint(S),eps,atol)) is_derivative = .true.

  endfunction is_derivative

  !!<summary>Checks whether a point is a lattice point of the given basis.</summary>
  !!<parameter name="LV" regular="true">The lattice vectors.</parameter>
  !!<parameter name="eps" regular="true">The finite precision tolerance.</parameter>
  !!<parameter name="pt" regular="true">The location of the point.</parameter>
  FUNCTION is_lattice_point(LV,pt,eps)
    real(dp) :: LV(3,3), pt(3), LVinv(3,3),eps
    real(dp) :: lattcoords(3)
    logical :: is_lattice_point, err
    real(dp) :: atol  ! An absolute tolerance for the "equal" function

    atol = 5E-4_dp
    is_lattice_point = .false.
    call matrix_inverse(LV,LVinv,err)
    if(err) stop "Problem with matrix inverse in function 'is_lattice_point'"
    lattcoords = matmul(LVinv,pt)
    if (equal(lattcoords,nint(lattcoords),eps,atol)) is_lattice_point = .true.

  END FUNCTION is_lattice_point

  !!<summary>Checks whether a group of points are lattice points of the given
  !!multilattice.</summary>
  !!<parameter name="LV" regular="true">The lattice vectors for the unit cell.</parameter>
  !!<parameter name="pt" regular="true">The lattice points inside the unit cell.</parameter>
  !!<parameter name="pB" regular="true">Points to be checked to see if they are lattice
  !!points.</parameter>
  !!<parameter name="eps" regular="true">The finite precision tolerance.</parameter>
  FUNCTION are_lattice_points(LV,pt,pB,eps)
    real(dp) :: LV(3,3), pt(:,:), LVinv(3,3),eps,pB(:,:)
    real(dp) :: lattcoords(size(pt,1),size(pt,2))
    logical :: flag(size(pt,2)), are_lattice_points, err
    integer :: iB, nB, nPt, iPt
    real(dp) :: atol  ! An absolute tolerance for the "equal" function

    atol = 5E-4_dp
    are_lattice_points = .false.; flag = .false.
    call matrix_inverse(LV,LVinv,err)
    if(err) stop "Problem with matrix inverse in function 'are_lattice_points'"
    nB = size(pB,2); nPt = size(pt,2)
    do iB = 1, nB
       lattcoords = matmul(LVinv,pt-spread(pB(:,iB),2,nPt))
       do iPt = 1, nPt
          flag(iPt) = flag(iPt) .or. equal(lattcoords(:,iPt),nint(lattcoords(:,iPt)),eps,atol)
       enddo
    enddo
    if (all(flag)) are_lattice_points = .true.
  END FUNCTION are_lattice_points

  !!<summary>This function checks if there is a single translation that maps every atom in
  !!  cell2 to a corresponding matching atom in cell1. If this is true, and we
  !!  already determined that the two lattices are equivalent, then the
  !!  structures must be equivalent. This test won't work unless the two atomic bases
  !!  are both inside the same unit cell.</summary>
  !!<parameter name="LV2" regular="true">The lattice vectors.</parameter>
  !!<parameter name="aPos1in" regular="true">Position of first atom.</parameter>
  !!<parameter name="aPos2" regular="true">Position of second atom.</parameter>
  !!<parameter name="aTyp1" regular="true">Type of first atom.</parameter>
  !!<parameter name="aTyp2" regular="true">Type of the second atom.</parameter>
  !!<parameter name="eps" regular="true">The finite precision tolerance.</parameter>
  function mapping_translation_exists(LV2, aPos1in, aPos2, aTyp1, aTyp2, eps)
    real(dp), intent(in) :: LV2(3,3), aPos1in(:,:), aPos2(:,:), eps
    logical :: mapping_translation_exists
    integer :: aTyp1(:), aTyp2(:)

    real(dp) :: transVec(3), conVec(3) ! The translation vector; connection vector
    integer N  ! Number of atoms in the basis (assumed same for both structures)
    integer ishift, iatom  ! loop counters
    real(dp) :: L2C(3,3), C2L(3,3) ! basis conversion matrices
    logical mapped ! true if there are matching atoms at position ConVec
    real(dp) :: aPos1(size(aPos1in,1),size(aPos1in,2)) ! str1 atoms in str2's cell

    mapping_translation_exists = .false.
    N=size(aPos1,2)  ! Total number of atoms in the cell
    if (N/=size(aPos2,2)) stop "cells must contain the same # of atoms"
    call get_transformations(LV2,L2C,C2L)
    aPos1 = aPos1in
    do iatom = 1, N ! Move the atom positions of str1 into the cell of str2
       call bring_into_cell(aPos1(:,iatom),C2L,L2C,eps)
    enddo

    do ishift = 1,N ! Loop over all possible translations
       ! Since the first atom in cell 2 must be mapped to one of the
       ! atoms in cell 1 if they are equivalent, try all vectors between
       ! the first atom in cell 2 and any atom in cell 1.
       transVec = aPos1(:,ishift)-aPos2(:,1) ! Trial translation
       do iatom = 1,N ! Loop over each atom in the basis
          conVec = transVec+aPos2(:,iatom) ! Connection vector to test
          call bring_into_cell(conVec,C2L,L2C,eps)
          call does_mapping_exist(conVec,aTyp2(iatom),aPos1,aTyp1,mapped,eps)
          if(.not. mapped) exit! if any atom doesn't map, try a new trans. vector
       enddo
       if(mapped) exit ! all atoms translated to matching atoms in the other structure
    enddo
    if(mapped) mapping_translation_exists = .true.
  endfunction mapping_translation_exists

  !!<summary>Subroutine to get the permutations of the basis vectors.</summary>
  !!<parameter name="uqlist" regular="true">A list of the unique basis.</parameter>
  !!<parameter name="perm">The permutation of the basis.</parameter>
  subroutine get_basis_permutations(uqlist,perm)
    integer, intent(in)  ::  uqlist(:)
    integer, pointer :: perm(:,:)

    integer Nq, Np ! Number of unique elements, # of permutations
    integer ip     ! counter for # of permutations
    integer  :: pl(size(uqlist))! permutation list
    integer j, l, temp    ! loop counter and swap temporary

    Nq = size(uqlist)
    pl = (/(j,j=1,Nq)/)

    ! Find factorial(Nq)
    Np = 1
    do ip=2,Nq
       Np = Np*ip
    enddo
    allocate(perm(Np,Nq))

    ! Generate all permutations of pl (and store in perms)
    ip=0
    do; ip=ip+1;
       perm(ip,:)=pl
       ! index rightmost element that is in *ascending* order
       do j=Nq-1,0,-1; if(j==0)exit; if (pl(j) < pl(j+1)) exit; enddo
       if (j==0) exit ! If there aren't any then all permutations generated
       ! index rightmost element larger than j-th one
       do l=Nq,1,-1; if (pl(j) < pl(l)) exit; enddo
       ! Swap j-th and l-th elements
       temp = pl(j);
       pl(j) = pl(l);
       pl(l) = temp;
       ! reverse the order of all of the elements right of the j-th one
       pl(j+1:Nq)=pl(Nq:j+1:-1);
    enddo
    if (.not. ip==Np) stop "Bug in 'get_basis_permutations': wrong # generated"
  endsubroutine get_basis_permutations

  !!<summary>Removes all the duplicates from a list of integers to find the unique ones.
  !!</summary>
  !!<parameter name="list" regular="true">List of atom labels.</parameter>
  !!<parameter name="uqlist">List of unique atom labels.</parameter>
  !!<parameter name="block">Order that the atom labels occur in.</parameter>
  subroutine unique(list,uqlist,block)
    integer, intent(in) :: list(:)
    integer, pointer :: uqlist(:),block(:)

    integer CurMin ! Current minimum
    integer iq     ! Look counter over unique values in list
    integer Nq     ! Number of unique values in list

    ! Count the number of each unique integer in a list
    ! and make a list of the unique numbers
    ! First loop counts so we can do the allocate...
    Nq = 0; CurMin = minval(list)-1
    do
       CurMin = minval(list,MASK=list>CurMin)
       Nq = Nq + 1
       if (CurMin==maxval(list))exit
    enddo
    allocate(uqlist(Nq),block(Nq))
    ! Now store unique values (uqlist) and the # of each (block)
    CurMin = minval(list)-1
    do iq = 1,Nq
       CurMin=minval(list,MASK=list>CurMin)
       block(iq) = count(list==CurMin)
         uqlist(iq) = CurMin
    enddo
    if(sum(block)/=size(list)) stop "Bug in routine 'unique'"
  endsubroutine unique

END MODULE compare_structures
