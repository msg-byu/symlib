MODULE vector_matrix_utilities
use num_types
use utilities
use numerical_utilities

implicit none
private
public matrix_inverse, determinant, cross_product, volume, norm, reduce_C_in_ABC, &
     orthogonality_defect, minkowski_conditions_check, norm_real_vector, &
     norms_real_vector_list, determinant_real, determinant_integer,&
     minkowski_reduce_basis, find_value_in_array

INTERFACE determinant
   module procedure determinant_real, determinant_integer
END INTERFACE

INTERFACE norm
   module procedure norm_real_vector, norms_real_vector_list
END INTERFACE norm

INTERFACE find_value_in_array
  module procedure find_lvalue_in_larray1
  module procedure find_intvalue_in_intarray3
  module procedure find_sintvalue_in_sintarray3
  module procedure find_intvalue_in_intarray4
  module procedure find_sintvalue_in_sintarray4
END INTERFACE find_value_in_array

CONTAINS

  !!<summary>This function checks the minkowski conditions for a 3D
  !!lattice basis.</summary>
  !!<parameter name="basis" regular="true"></parameter>
  !!<parameter name="eps" regular="true"></parameter>
  function minkowski_conditions_check(basis,eps)
    logical  :: minkowski_conditions_check
    ! real(dp), dimension(3,3), intent(IN) :: basis
    real(dp), intent(IN) :: basis(3,3)
    real(dp), intent(IN) :: eps
    real(dp), dimension(3) :: b1, b2, b3

    b1 = basis(:,1)
    b2 = basis(:,2)
    b3 = basis(:,3)

    minkowski_conditions_check = .true.
    if (norm(b1) > norm(b2)+eps) then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 1 failed: b1 > b2"
    endif
    if (norm(b2) > norm(b3)+eps)        then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 2 failed: b2 > b3"
    endif
    if (norm(b2) > norm(b1+b2)+eps)     then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 3 failed: b2 > b1+b2"
    endif
    if (norm(b2) > norm(b1-b2)+eps)     then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 4 failed: b2 > b1-b2"
    endif
    if (norm(b3) > norm(b1+b3)+eps)     then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 5 failed: b3 > b1+b3"
    endif
    if (norm(b3) > norm(b3-b1)+eps)     then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 6 failed: b3 > b3-b1"
    endif
    if (norm(b3) > norm(b2+b3)+eps)     then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 7 failed: b3 > b2+b3"
    endif
    if (norm(b3) > norm(b3-b2)+eps)     then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 8 failed: b3 > b3-b2"
    endif
    if (norm(b3) > norm(b1+b2+b3)+eps)  then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 9 failed: b3 > b1+b2+b3"
    endif
    if (norm(b3) > norm(b1-b2+b3)+eps)  then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 10 failed: b3 > b1-b2+b3"
    endif
    if (norm(b3) > norm(b1+b2-b3)+eps)  then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 11 failed: b3 > b1+b2-b3"
    endif
    if (norm(b3) > norm(b1-b2-b3)+eps)  then
       minkowski_conditions_check = .false.
       write(*,*) "Minkowski_condition 12 failed: b3 > b1-b2-b3"
    endif
  endfunction minkowski_conditions_check

  !!<summary>This routine takes three vectors, A,B,C, defining a
  !!lattice, and reduces the last one so that it is as close as
  !!possible to the origin while remaining in an affine plane, which
  !!is parallel to the A-B plane but which passes through the end of
  !!the C vector. See Lecture notes in computer science, ISSN
  !!0302-974, ANTS - VI : algorithmic number theory, 2004, vol. 3076,
  !!pp. 338-357 ISBN 3-540-22156-5</summary>
  !!<parameter name="A" regular="true"></parameter>
  !!<parameter name="B" regular="true"></parameter>
  !!<parameter name="C" regular="true"></parameter>
  !!<parameter name="eps" regular="true"></parameter>
  subroutine reduce_C_in_ABC(A,B,C,eps)
    real(dp), intent(inout) :: A(3), B(3), C(3)
    real(dp), intent(in) :: eps
    real(dp), dimension(3) :: T  ! projection of C into the A-B plane
    real(dp), dimension(3,3) :: ABC, ABCinv, oldABC ! Matrices of ABC basis vectors and inverse
    real(dp), dimension(3)   :: cpdAB ! unit vector perpendicular to the A-B plane
    real(dp) :: dist(4) ! the distances from T to enclosing lattice
                        ! points of A,B (4 corners of the ppiped)
    integer LC(3) ! lattice coordinates of C, in the affine plane, using the A,B basis vectors
    integer idx(1) ! index of the smallest distance from T to a lattice point in A,B
    logical err
    real(dp) :: corners(4,3)  ! Array to store 4 corners of the unit cell (in lattice coordinates)
    integer i ! loop counter

    ABC = reshape((/A,B,C/),(/3,3/))
    oldABC = ABC

    ! Use Gaussian reduction to reduce the A,B 2D basis so that it is
    ! itself Minkowski reduced. If this is done, then the closest
    ! lattice point (in A,B plane) to the projection of C (into the
    ! A,B plane) is guaranteed to be one of the corners of the unit
    ! cell enclosing the projection of C
    call gaussian_reduce_two_vectors(A,B,eps)

    ! First thing to do is find the (real, not lattice) point in the affine plane A,B + C that is
    ! nearest the origin. Call this T.
    cpdAB = cross_product(A,B)/norm(cross_product(A,B))
    T = C - cpdAB*dot_product(C,cpdAB)

    if(.not. equal(dot_product(T,cross_product(A,B)),0._dp,eps)) then
       print *,dot_product(T,cross_product(A,B))
       stop "Projection of C into A,B plane failed"
    endif

    ! Now find the four points of the A,B lattice, in the affine plane, that enclose the point T
    ABC = reshape((/A,B,C/),(/3,3/))

    call matrix_inverse(ABC,ABCinv,err)
    if(err)stop "A,B,C vectors in reduce_C_in_ABC are co-planar"
    LC = floor(matmul(ABCinv,T) + eps)

!    dist(1) = norm(T-matmul(ABC,LC))
!    dist(2) = norm(T-matmul(ABC,(/LC(1)+1,LC(2),LC(3)/)))
!    dist(3) = norm(T-matmul(ABC,(/LC(1),LC(2)+1,LC(3)/)))
!    dist(4) = norm(T-matmul(ABC,(/LC(1)+1,LC(2)+1,LC(3)/)))
!
! GLWH July 10 2015. This stupid do loop business is needed to circumvent a compiler bug in ifort
    ! -O3 version  14.0.2 20140121. Old code is above.

    ! Compute the distance from T to each of the four corners of the cell and pick
    ! the one that is the closest.
    corners(1,:) =(/0,0,0/)
    corners(2,:) =(/1,0,0/)
    corners(3,:) =(/0,1,0/)
    corners(4,:) =(/1,1,0/)
    do i = 1,4
       dist(i) = norm(T-matmul(ABC,LC+corners(i,:)))
    enddo
    idx = minloc(dist)

    select case(idx(1))
    case(1)
       C = C - matmul(ABC,LC+corners(1,:))
    case(2)
       C = C - matmul(ABC,LC+corners(2,:))
    case(3)
       C = C - matmul(ABC,LC+corners(3,:))
    case(4)
       C = C - matmul(ABC,LC+corners(4,:))
    case default
       print *, "Case failed in reduce_C_in_ABC"
       write(*,'("Lattice coordinates in the A,B plane: ",2(i2,1x))') LC
       stop
    end select

    ABC = reshape((/A,B,C/),(/3,3/))
    call matrix_inverse(ABC,ABCinv,err)
    if(any(abs(matmul(ABCinv,oldABC)-nint(matmul(ABCinv,oldABC)))>eps)) stop "Lattice was not preserved &
         & in reduce_A_in_ABC"

    if(err)stop "A,B,C vectors in reduce_C_in_ABC are co-planar (after Minkowski)"

  endsubroutine reduce_C_in_ABC

  !!<summary>This routine takes two vectors (in three-space) and
  !!reduces them to form a shortest set (Minkowski reduced). The idea
  !!is to subtract multiples of U from V so that the new V is as close
  !!to the origin as any lattice point along the line that passes
  !!through U in the direction of V. The process is repeated until the
  !!new vector isn't shorter than the other. It's pretty obvious if
  !!you do an example by hand. Also see 3.1 of Lecture notes in
  !!computer science, ISSN 0302-974, ANTS - VI: algorithmic number
  !!theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5. Fixes
  !!made Apr 2012 GLWH (not sure if they made a practical difference
  !!though)</summary>
  !!<parameter name="U" regular="true"></parameter>
  !!<parameter name="V" regular="True"></parameter>
  !!<parameter name="eps" regular="true"></parameter>
  subroutine gaussian_reduce_two_vectors(U,V,eps)
    real(dp) :: U(3), V(3), R(3)
    real(dp), intent(in) :: eps

    real(dp) temp(3)
    integer it

    it = 0
    if (norm(U) > norm(V) - eps) then
       ! Make sure that the {U,V} are listed in ascending order; ||U||<||V||
       temp = U; U = V; V = temp ! Keep V as the longest vector
    end if
    do; it = it + 1
       if (it > 10) stop "gaussian_reduce_two_vectors failed to converge in 10 iterations"
       R = V-nint(dot_product(U,V)/dot_product(U,U))*U !Shorten V as much as possible
       V = U ! Swap U and V (so U remains the shortest)
       U = R
       if (norm(U) >= norm(V) - eps) exit
    end do
    ! Make sure that the {U,V} are listed in ascending order on exit; ||U||<||V||
    temp = U; U = V; V = temp

  endsubroutine gaussian_reduce_two_vectors

  !!<summary></summary>
  !!<parameter name="IN" regular="true"></parameter>
  !!<parameter name="OUT" regular="true"></parameter>
  !!<parameter name="eps" regular="true"></parameter>
  !!<parameter name="aeps_" regular="true"></parameter>
  SUBROUTINE minkowski_reduce_basis(IN, OUT, eps, aeps_)
    real(dp), intent(in) :: IN(3,3)
    real(dp), intent(out) :: OUT(3,3)
    real(dp), intent(in) :: eps
    real(dp), optional, intent(in) :: aeps_
    real(dp)             :: norms(3), temp(3,3), aeps
    integer              :: i, idx(1), it, limit

    if (present(aeps_)) then; aeps = aeps_; else; aeps = 5e-5_dp; endif
    limit = 10
    if (equal(determinant(IN),0._dp,eps, atolerance_=aeps)) stop "Input basis for 'minkowski_reduce_basis' was not linearly independent"
    OUT = IN

    ! print*,"in mink routine"
    ! Keep applying the greedy algorithm until the vectors come out already sorted
    do it = 1, limit
       ! Sort the three vectors into ascending order
       temp = OUT
       norms = norm(temp)
       do i = 3,1,-1
          idx = maxloc(norms)
          temp(:,i) = OUT(:,idx(1))
          norms(idx(1)) = 0
       enddo
       OUT = temp ! Copy the sorted vectors back to OUT
       !   write(*,'(f7.3)') orthogonality_defect(OUT)
       call reduce_C_in_ABC(OUT(:,1),OUT(:,2),OUT(:,3),eps)
       !   write(*,'(f7.3)') orthogonality_defect(OUT)
       if (norm(OUT(:,3))>=norm(OUT(:,2))-eps) exit

    end do
    ! print *,"mink past first loop"
    !if (it>limit+1) stop "Too many iterations in 'minkowski_reduce_basis'"
    if (.not. minkowski_conditions_check(OUT,eps)) then
       write(*,'("ERROR in minkowski_reduce_basis: Minkowski conditions not met.")')
       write(*,'("Number of iterations:",i3)') limit
       stop
    end if

    ! we want to make sure that the det is positive.
    ! NOTE: This *destroys* the mathematical picture of a "greedy reduced basis" (Minkowski), but
    !       from a physical point of view we don't care ;-)
    !       Either way, the basis is as orthogonal as possible.
    if (determinant(OUT)<0) then
       temp(:,1) = OUT(:,2)
       OUT(:,2) = OUT(:,3)
       OUT(:,3) = temp(:,1)
    endif

  END SUBROUTINE minkowski_reduce_basis

  !!<summary>This function calculates the "orthogonality defect" of
  !!the given basis of a 3D lattice.</summary>
  !!<parameter name="basis" regular="true"></parameter>
  function orthogonality_defect(basis)
    real(dp), intent(in) :: basis(3,3)
    real(dp) :: orthogonality_defect
    real(dp) :: od
    integer j
    od = 1._dp
    do j=1,3
       od = od*norm(basis(:,j))
    enddo
    od = od/abs(determinant(basis))
    orthogonality_defect = od
  endfunction orthogonality_defect

  !!<summary>Given the matrix a, finds its inverse b</summary>
  !!<parameter name="a" regular="true"></parameter>
  !!<parameter name="b" regular="true"></parameter>
  !!<parameter name="err_" regular="true"></parameter>
  !!<parameter name="eps_" regular="true"></parameter>
  subroutine matrix_inverse(a, b, err_, eps_)
    real(dp), intent(in):: a(3,3)
    real(dp),intent(out):: b(3,3)
    real(dp) :: c,avec(9)
    real(dp)           :: eps

    logical,  optional :: err_
    real(dp), optional :: eps_

    if(present(err_)) err_ = .false.
    if(present(eps_)) then; eps=eps_; else; eps=10d-14; endif

    c=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
         a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
         a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
    avec=(/-a(2,3)*a(3,2)+a(2,2)*a(3,3),a(1,3)*a(3,2)-a(1,2)*a(3,3), &
         -a(1,3)*a(2,2)+a(1,2)*a(2,3),a(2,3)*a(3,1)-a(2,1)*a(3,3), &
         -a(1,3)*a(3,1)+a(1,1)*a(3,3),a(1,3)*a(2,1)-a(1,1)*a(2,3), &
         -a(2,2)*a(3,1)+a(2,1)*a(3,2),a(1,2)*a(3,1)-a(1,1)*a(3,2), &
         -a(1,2)*a(2,1)+a(1,1)*a(2,2)/)
    if(abs(c) < eps) then; if(present(err_)) err_ = .true.
    else; b=1/c*transpose(reshape(avec,(/3,3/))); endif
    ! transpose due to column major storage
  end subroutine matrix_inverse

  !!<summary>Given vectors a and b, c = a x b</summary>
  !!<parameter name="a" regular="true"></parameter>
  !!<parameter name="b" regular="true"></parameter>
  function cross_product(a,b)
    real(dp) :: a(3), b(3), cross_product(3)

    cross_product = (/a(2)*b(3) - a(3)*b(2), &
         a(3)*b(1) - a(1)*b(3), &
         a(1)*b(2) - a(2)*b(1)/)
  end function cross_product

  !!<summary>This routine takes the norm of a vector</summary>
  !!<parameter name="vector" regular="true"></parameter>
  pure function norm_real_vector(vector)
    real(dp) :: norm_real_vector
    real(dp), intent(in):: vector(:)

    norm_real_vector = sqrt(dot_product(vector,vector))
  end function norm_real_vector

  !!<summary>This routine takes the norms of vectors in a list</summary>
  !!<parameter name="vector_list" regular="true"></parameter>
  pure function norms_real_vector_list(vector_list)
    real(dp), intent(in)               :: vector_list(:,:)
    real(dp) :: norms_real_vector_list(size(vector_list,2))
    integer :: i

    do i = 1, size(vector_list,2)
       norms_real_vector_list(i) = sqrt(dot_product(vector_list(:,i),vector_list(:,i)))
    end do
  end function norms_real_vector_list

  !!<summary>This routine finds the determinant of 2x2 or 3x3 matrices</summary>
  !!<parameter name="a" regular="true"></parameter>
  function determinant_real(a)
    real(dp) :: determinant_real,a(:,:)

    integer :: n
    n=size(a,dim=1)
    if (n==2) then
       determinant_real=-a(1,2)*a(2,1)+a(1,1)*a(2,2)
    else if (n==3) then
       determinant_real=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
            a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
            a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
    else
       stop "DETERMINANT ERROR: Matrix dimension exceeds 3."
    end if
  end function determinant_real

  !!<summary>This routine finds the determinant of 2x2 or 3x3 matrices</summary>
  !!<parameter name="a" regular="true"></parameter>
  function determinant_integer(a)
    integer :: determinant_integer,a(:,:)
    integer :: n
    n=size(a,dim=1)
    if (n==2) then
       determinant_integer=-a(1,2)*a(2,1)+a(1,1)*a(2,2)
    else if (n==3) then
       determinant_integer=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
            a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
            a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
    else
       stop "DETERMINANT ERROR: Matrix dimension exceeds 3."
    end if
  end function determinant_integer


  !!<summary>This function takes three vectors and returns a "signed"
  !!volume of the parallelpiped that they form.</summary>
  !!<parameter name="a1" regular="true"></parameter>
  !!<parameter name="a2" regular="true"></parameter>
  !!<parameter name="a3" regular="true"></parameter>
  function volume(a1, a2, a3)
    real(dp) :: volume
    real(dp) :: a1(3), a2(3), a3(3)

    volume = dot_product(a1, cross_product(a2,a3))
  end function volume

  !!<summary>Find a int value inside a 3-dim array</summary>
  !!<parameter name="array" regular="true"></parameter>
  !!<parameter name="vals" regular="true"></parameter>
  !!<parameter name="loc" regular="true">intent(out), { #found_value }</parameter>
  !!<parameter name="nfound" regular="true">number of val occurrences in array</parameter>
  subroutine find_lvalue_in_larray1(array,vals,loc,nfound)
    logical, intent(in)     :: array(:)
    logical, intent(in)     :: vals(:)
    integer, intent(out)    :: loc(:)
    integer, intent(out)    :: nfound

    integer, parameter :: dim=1
    integer :: count, l

    count = 0
    do l=1,size(array,1)
       if (any(array(l).eqv.vals(:))) then
          count = count+1
          loc(count) = l
       endif
    enddo

    nfound=count

  end subroutine find_lvalue_in_larray1



  !!<summary>Find a int value inside a 3-dim array</summary>
  !!<parameter name="array" regular="true"></parameter>
  !!<parameter name="vals" regular="true"></parameter>
  !!<parameter name="loc">intent(out), { #coordinate, #found_value }</parameter>
  !!<parameter name="nfound" regular="true">number of val occurrences in array</parameter>
  subroutine find_intvalue_in_intarray3(array,vals,loc,nfound)
    integer, intent(in)     :: array(:,:,:)
    integer, intent(in)     :: vals(:)
    integer, pointer        :: loc(:,:)
    integer, intent(out)    :: nfound

    integer, parameter :: dim=3

    integer :: j,k,l
    integer :: count

    count = 0
    do j=1,size(array,3)
       do k=1,size(array,2)
          do l=1,size(array,1)
             if (any(array(l,k,j)==vals(:))) then
                count = count+1
                loc(:dim,count) = (/l,k,j/)
             endif
          enddo
       enddo
    enddo

    nfound=count

  end subroutine find_intvalue_in_intarray3

  !!<summary>Find a int value inside a 3-dim array, si version</summary>
  !!<parameter name="array" regular="true"></parameter>
  !!<parameter name="vals" regular="true"></parameter>
  !!<parameter name="loc">intent(out), { #coordinate, #found_value }</parameter>
  !!<parameter name="nfound" regular="true">number of val occurrences in array</parameter>
  subroutine find_sintvalue_in_sintarray3(array,vals,loc,nfound)
    integer(si), intent(in)     :: array(:,:,:)
    integer(si), intent(in)     :: vals(:)
    integer, pointer        :: loc(:,:)
    integer, intent(out)    :: nfound

    integer, parameter :: dim=3

    integer :: j,k,l
    integer :: count

    count = 0
    do j=1,size(array,3)
       do k=1,size(array,2)
          do l=1,size(array,1)
             if (any(array(l,k,j)==vals(:))) then
                count = count+1
                loc(:dim,count) = (/l,k,j/)
             endif
          enddo
       enddo
    enddo

    nfound=count

  end subroutine find_sintvalue_in_sintarray3


  !!<summary>Find a int value inside a 4-dim array</summary>
  !!<parameter name="array" regular="true"></parameter>
  !!<parameter name="vals" regular="true"></parameter>
  !!<parameter name="loc">intent(out), { #coordinate, #found_value }</parameter>
  !!<parameter name="nfound" regular="true">number of val occurrences in array</parameter>
  subroutine find_intvalue_in_intarray4(array,vals,loc,nfound)
    integer, intent(in)     :: array(:,:,:,:)
    integer, intent(in)     :: vals(:)
    integer, pointer        :: loc(:,:) ! intent(out), { #coordinate, #found_value }
    integer, intent(out)    :: nfound   ! number of val occurrences in array

    integer, parameter :: dim=4

    integer :: i,j,k,l
    integer :: count

    count = 0
    do i=1,size(array,4)
       do j=1,size(array,3)
          do k=1,size(array,2)
             do l=1,size(array,1)
                if (any(array(l,k,j,i)==vals(:))) then
                   count = count+1
                   loc(:dim,count) = (/l,k,j,i/)
                endif
             enddo
          enddo
       enddo
    enddo

    nfound=count

  end subroutine find_intvalue_in_intarray4

  !!<summary>Find a int value inside a 4-dim array, si version</summary>
  !!<parameter name="array" regular="true"></parameter>
  !!<parameter name="vals" regular="true"></parameter>
  !!<parameter name="loc">intent(out), { #coordinate, #found_value }</parameter>
  !!<parameter name="nfound" regular="true">number of val occurrences in array</parameter>
  subroutine find_sintvalue_in_sintarray4(array,vals,loc,nfound)
    integer(si), intent(in)     :: array(:,:,:,:)
    integer(si), intent(in)     :: vals(:)
    integer, pointer        :: loc(:,:)
    integer, intent(out)    :: nfound

    integer, parameter :: dim=4

    integer :: i,j,k,l
    integer :: count

    count = 0
    do i=1,size(array,4)
       do j=1,size(array,3)
          do k=1,size(array,2)
             do l=1,size(array,1)
                if (any(array(l,k,j,i)==vals(:))) then
                   count = count+1
                   loc(:dim,count) = (/l,k,j,i/)
                endif
             enddo
          enddo
       enddo
    enddo

    nfound=count

  end subroutine find_sintvalue_in_sintarray4

END MODULE
