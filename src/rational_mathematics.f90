MODULE rational_mathematics
  use num_types
  use vector_matrix_utilities
  use utilities, only: ralloc

implicit none
private
public gcd, SmithNormalForm, HermiteNormalForm, is_a_rational_in_range,&
       SmithNormalForm_li, get_rationals_in_range

! Overloaded procedure for computing the greatest common denominator
INTERFACE gcd
   MODULE PROCEDURE gcd_2ints, gcd_rank1, gcd_3ints, gcd_4ints
END INTERFACE

CONTAINS

  !!<summary>This routine takes an integer 3x3 matrix and computes its
  !!Smith Normal Form.</summary>
  !!<parameter name="H" regular="True">Input matrix.</parameter>
  !!<parameter name="A" regular="True">Left Transform.</parameter>
  !!<parameter name="M" regular="true">Smith Normal Form matrix.</parameter>
  !!<parameter name="B" regular="True">Right Transform.</parameter>
  !!<parameter name="err_" regular="True">Returns 1 if overflow occures.</parameter>
  subroutine SmithNormalForm_li(H,A,M,B,err_)
    integer, intent(in) :: H(3,3)
    integer, intent(out) :: M(3,3)
    integer(li), intent(out), dimension(3,3) :: A, B
    integer, optional, intent(out) :: err_

    integer :: i, row, col, min_val, j, check(9)
    integer :: multiple, tmpVec(3)
    integer(li) :: tmpVec_li(3)
    logical :: is_snf, new_pivot, OverFlowCheck
    integer :: itCnt, min_row

    OverFlowCheck = .False.
    if (present(err_)) OverFlowCheck = .True.

    if(determinant(H)<1) stop "SmithNormalForm routine failed because the input matrix had a negative determinant"
    A = 0; B = 0; M = H ! M starts out as H, the input matrix
    forall(i=1:3); A(i,i) = 1; B(i,i) = 1; end forall ! A & B = identity

    j=1
    itCnt = 0
    is_snf = .False.
    new_pivot = .True.
    do while (is_snf .eqv. .False. .and. j<4)
       itCnt = itCnt + 1
       if (itCnt>=100) stop "ERROR bad programming in SmithNormalForm"

       if (new_pivot) then
          call get_min_val(M, j, min_val, row, col)
       end if

       do i=1,3
          if (i==col) cycle
          multiple = nint(real(M(row,i),dp)/real(min_val,dp))
          if (multiple==0) cycle
          M(:,i) = M(:,i)-multiple*M(:,col)
          B(:,i) = B(:,i)-multiple*B(:,col)
       end do

       do i=1,3
          if (i==row) cycle
          multiple = nint(real(M(i,col),dp)/real(min_val,dp))
          if (multiple==0) cycle
          M(i,:) = M(i,:)-multiple*M(row,:)
          A(i,:) = A(i,:)-multiple*A(row,:)
       end do

       new_pivot = .True.
       if ((count(M(:,col)==0)==2) .and. (count(M(row,:)==0)==2)) then
          if (all(mod(M(j:,j:),min_val)==0)) then
             if (j < col) then
                tmpVec_li = B(:,j); B(:,j) = B(:,col); B(:,col) = tmpVec_li
                tmpVec = M(:,j); M(:,j) = M(:,col); M(:,col) = tmpVec
             end if
             if (j < row) then
                tmpVec_li = A(j,:); A(j,:) = A(row,:); A(row,:) = tmpVec_li
                tmpVec = M(j,:); M(j,:) = M(row,:); M(row,:) = tmpVec
             end if
             j = j + 1
          else
             new_pivot = .False.
             call get_min_loc(M, min_val, j, min_row)
             M(row,:) = M(row,:) + M(min_row,:)
             A(row,:) = A(row,:) + A(min_row,:)
          end if
       end if

       check = reshape(M,(/9/))
       if (all(check((/2,3,4,6,7,8/))==0) .and. all(check((/1,5,9/))/=0)) then
          if (mod(M(2,2),M(1,1))==0 .and. mod(M(3,3),M(2,2))==0) then
             is_snf = .True.
          end if
       end if
    end do

    do i=1,3
       if (M(i,i) < 0) then
          M(i,:) = -M(i,:)
          A(i,:) = -A(i,:)
       end if
    end do

    if (any(matmul(matmul(A,H),B)/=M)) stop "END: Transformation matrices didn't work"
    check = reshape(M,(/9/))
    if (any(check((/2,3,4,6,7,8/))/=0)) stop "Not diagonal"
    if (mod(M(2,2),M(1,1))/=0 .or. mod(M(3,3),M(2,2))/=0) stop "SNF conditions not met"
    if (OverFlowCheck) then
       if ((any(abs(real(A,dp)) > 1E17)) .or. (any(abs(real(B,dp)) > 1E17))) then
          write(*,*) "Warning Values in SmithNormalForm overflowing standard ints."
          err_ = 1
       else
          err_ = 0
       end if
    else
       if ((any(abs(real(A,dp)) > 1E17)) .or. (any(abs(real(B,dp)) > 1E17))) stop "Warning Values in SmithNormalForm overflowing standard ints."
    end if
  ENDSUBROUTINE SmithNormalForm_Li

  !!<summary>Finds the minimal value in the sub matrix of A where the
  !!sub matrix consits of every row and colum greater than or equal to
  !!diag.</summary>
  !!<parameter name="A" regular="true">The input matrix.</parameter>
  !!<parameter name="diag" regular="true">The row/column number we
  !!want to be below.</parameter>
  !!<parameter name="min_val" regular="true">The smallest value in A
  !!that's not in an empty row or column.</parameter>
  !!<parameter name="row" regular="true">The row with the lowest value
  !!in it.</parameter>
  !!<parameter name="col" regular="true">The column with the lowest
  !!value in it.</parameter>
  subroutine get_min_val(A, diag, min_val, row, col)
    integer, intent(in) :: A(3,3), diag
    integer, intent(out) :: min_val, row, col

    integer :: temp_A(3,3), i, tmploc(2)
    logical :: found

    found = .False.
    temp_A = abs(A)
    do i=1,9
       tmploc = minloc(temp_A, temp_A>0)
       row = tmploc(1)
       col = tmploc(2)
       min_val = A(row,col)
       if (row >= diag .and. col >= diag) then
          found = .True.
          exit
       else
          temp_A(row, col) = -1
       end if
    end do

    if (found .eqv. .False.) stop "Failed to find minimal value in get_min_val."
  end subroutine get_min_val

  !!<summary>Finds the row of the smallest value not divisible by a
  !!pivot within the sub matrix of A consisting of the row and column
  !!numbers greater than or equal to diag.</summary>
  !!<parameter name="A" regular="true">The input matrix.</parameter>
  !!<parameter name="pivot" regular="true">The pivot.</parameter>
  !!<parameter name="diag" regular="true">The row/column number we
  !!want to be below.</parameter>
  !!<parameter name="row" regular="true">The row with the lowest value
  !!in it.</parameter>
  subroutine get_min_loc(A, pivot, diag, row)
    integer, intent(in) :: A(3,3), pivot, diag
    integer, intent(out) :: row

    integer :: temp_A(3,3), i, tmploc(2), mod_p(3,3), col
    logical :: found

    found = .False.
    temp_A = abs(A)
    mod_p = mod(A,pivot)
    do i=1,9
       tmploc = minloc(temp_A, temp_A>0)
       row = tmploc(1)
       col = tmploc(2)
       if (row >= diag .and. col >= diag .and. mod_p(row,col) /= 0) then
          found = .True.
          exit
       else
          temp_A(row, col) = -1
       end if
    end do

    if (found .eqv. .False.) stop "Failed to find minimal value in get_min_val."
  end subroutine get_min_loc

  !!<summary>This routine takes an integer 3x3 matrix and computes its
  !!Smith Normal Form.</summary>
  !!<parameter name="H" regular="True">Input matrix.</parameter>
  !!<parameter name="A" regular="True">Left Transform.</parameter>
  !!<parameter name="M" regular="true">Smith Normal Form matrix.</parameter>
  !!<parameter name="B" regular="True">Right Transform.</parameter>
  !!<parameter name="err_" regular="True">Returns 1 if overflow occures.</parameter>
  subroutine SmithNormalForm(H,A,M,B,err_)
    integer, intent(in) :: H(3,3)
    integer, intent(out), dimension(3,3) :: M, A, B
    integer, optional, intent(out) :: err_

    integer :: i, row, col, min_val, j, check(9)
    integer :: multiple, tmpVec(3)
    logical :: is_snf, new_pivot, OverFlowCheck
    integer :: itCnt, min_row

    OverFlowCheck = .False.
    if (present(err_)) OverFlowCheck = .True.

    if(determinant(H)<1) stop "SmithNormalForm routine failed because the input matrix had a negative determinant"
    A = 0; B = 0; M = H ! M starts out as H, the input matrix
    forall(i=1:3); A(i,i) = 1; B(i,i) = 1; end forall ! A & B = identity

    j=1
    itCnt = 0
    is_snf = .False.
    new_pivot = .True.
    do while (is_snf .eqv. .False. .and. j<4)
       itCnt = itCnt + 1
       if (itCnt>=100) stop "ERROR bad programming in SmithNormalForm"

       if (new_pivot) then
          call get_min_val(M, j, min_val, row, col)
       end if

       do i=1,3
          if (i==col) cycle
          multiple = nint(real(M(row,i),dp)/real(min_val,dp))
          if (multiple==0) cycle
          M(:,i) = M(:,i)-multiple*M(:,col)
          B(:,i) = B(:,i)-multiple*B(:,col)
       end do

       do i=1,3
          if (i==row) cycle
          multiple = nint(real(M(i,col),dp)/real(min_val,dp))
          if (multiple==0) cycle
          M(i,:) = M(i,:)-multiple*M(row,:)
          A(i,:) = A(i,:)-multiple*A(row,:)
       end do

       new_pivot = .True.
       if ((count(M(:,col)==0)==2) .and. (count(M(row,:)==0)==2)) then
          if (all(mod(M(j:,j:),min_val)==0)) then
             if (j < col) then
                tmpVec = B(:,j); B(:,j) = B(:,col); B(:,col) = tmpVec
                tmpVec = M(:,j); M(:,j) = M(:,col); M(:,col) = tmpVec
             end if
             if (j < row) then
                tmpVec = A(j,:); A(j,:) = A(row,:); A(row,:) = tmpVec
                tmpVec = M(j,:); M(j,:) = M(row,:); M(row,:) = tmpVec
             end if
             j = j + 1
          else
             new_pivot = .False.
             call get_min_loc(M, min_val, j, min_row)
             M(row,:) = M(row,:) + M(min_row,:)
             A(row,:) = A(row,:) + A(min_row,:)
          end if
       end if

       check = reshape(M,(/9/))
       if (all(check((/2,3,4,6,7,8/))==0) .and. all(check((/1,5,9/))/=0)) then
          if (mod(M(2,2),M(1,1))==0 .and. mod(M(3,3),M(2,2))==0) then
             is_snf = .True.
          end if
       end if
    end do

    do i=1,3
       if (M(i,i) < 0) then
          M(i,:) = -M(i,:)
          A(i,:) = -A(i,:)
       end if
    end do
    if (any(matmul(matmul(A,H),B)/=M)) stop "END: Transformation matrices didn't work"
    check = reshape(M,(/9/))
    if (any(check((/2,3,4,6,7,8/))/=0)) stop "Not diagonal"
    if (mod(M(2,2),M(1,1))/=0 .or. mod(M(3,3),M(2,2))/=0) stop "SNF conditions not met"
    if (OverFlowCheck) then
       if ((any(abs(A) > 1E9)) .or. (any(abs(B) > 1E9))) then
          write(*,*) "Warning Values in SmithNormalForm overflowing standard ints."
          err_ = 1
       else
          err_ = 0
       end if
    else
       if ((any(abs(A) > 1E9)) .or. (any(abs(B) > 1E9))) stop "Warning Values in SmithNormalForm overflowing standard ints."
    end if
  ENDSUBROUTINE SmithNormalForm

  !!<summary>Find the Hermite normal form of a a given integer
  !!matrix. Similar to the SNF finder above but a little simpler. This
  !!routine is not very elegant, just brute force. Don't be
  !!disappointed.</summary>
  !!<parameter name="S" regular="true">The 3x3 integer matrix
  !!describing the relationship between two commensurate
  !!lattices.</parameter>
  !!<parameter name="H" regular="true">The resulting HNF
  !!matrix.</parameter>
  !!<parameter name="B" regular="true">The transformation matrix such
  !!that $H = SB$.</parameter>
  subroutine HermiteNormalForm(S,H,B)
    integer, intent(in) :: S(3,3)
    integer, intent(out), dimension(3,3) :: H, B

    !!<local name="tempcol">When two columns need to be swapped,
    !!tempcol mediates the swap.</local>
    !!<local name="check">Stores the HNF as a vector; only used as
    !!failsafe for meeting HNF criteria.</local>
    !!<local name="maxidx, minidx">The indices of max/min values in a
    !!column.</local>
    !!<local name="multiple">When reducing the columns to HNF form,
    !!the multiplicative relationship between two elements of two
    !!columns; reused.</local>
    !!<local name="minm">The smallest element in row 1 of the HNF
    !!being constructed.</local>

    integer :: i, minm, maxidx, minidx, multiple, j,  check(9), tempcol(3)

    if (determinant(S) == 0) stop "Singular matrix passed to HNF routine"
    B = 0; H = S ! H starts out as S, the input matrix
    forall(i=1:3); B(i,i) = 1; end forall ! B = identity

    do ! Keep doing column operations until all elements in row 1 are
       ! zero except the one on the diagonal.
       ! Divide the column with the smallest value into the largest
       do while (count(H(1,:)/=0) > 1) ! Keep going until only zeros beyond first element
          call get_minmax_indices(H(1,:),minidx,maxidx)
          minm = H(1,minidx)
          ! Subtract a multiple of the column containing the smallest element from
          ! the row containing the largest element
          multiple = H(1,maxidx)/minm ! Factor to multiply by
          H(:,maxidx) = H(:,maxidx) - multiple*H(:,minidx)
          B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
          if (any(matmul(S,B)/=H)) stop "COLS: Transformation matrices didn't work"
       enddo ! End of step row 1
       if (H(1,1) == 0) call swap_column(H,B,1) ! swap columns if (1,1) is zero
       if (H(1,1)<0)then; H(:,1) = -H(:,1); B(:,1) = -B(:,1) ! Change sign if (1,1) is negative
       endif
       if (count(H(1,:)/=0) > 1) stop "Didn't zero out the rest of the row"
       if (any(matmul(S,B)/=H)) stop "COLSWAP: Transformation matrices didn't work"
       exit
    enddo
    ! Now work on element H(2,3)
    do
       do while (H(2,3)/=0)
          if (H(2,2) == 0) then
             tempcol = H(:,2); H(:,2) = H(:,3); H(:,3) = tempcol
             tempcol = B(:,2); B(:,2) = B(:,3); B(:,3) = tempcol
             if (H(2,3) == 0) exit
          endif
          if (abs(H(2,3))<abs(H(2,2))) then; maxidx = 2; minidx = 3
          else; maxidx = 3; minidx = 2; endif
          multiple = H(2,maxidx)/H(2,minidx)
          H(:,maxidx) = H(:,maxidx) - multiple*H(:,minidx)
          B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
          if (any(matmul(S,B)/=H)) stop "COLS: Transformation matrices didn't work"
       enddo
       if (H(2,2) == 0) then
          tempcol = H(:,2); H(:,2) = H(:,3); H(:,3) = tempcol
       endif
       if (H(2,2)<0) then ! Change signs
       H(:,2) = -H(:,2); B(:,2) = -B(:,2); endif
       if (H(2,3)/=0) stop "Didn't zero out last element"

       if (any(matmul(S,B)/=H)) stop "COLSWAP: Transformation matrices didn't work"
       exit
    enddo
    if (H(3,3)<0) then ! Change signs
       H(:,3) = -H(:,3); B(:,3) = -B(:,3);
    endif
    check = reshape(H,(/9/))
    if (any(check((/4,7,8/))/=0)) stop "Not lower triangular"
    if (any(matmul(S,B)/=H)) stop "END PART1: Transformation matrices didn't work"

    ! Now that the matrix is in lower triangular form, make sure the lower off-diagonal
    ! elements are non-negative but less than the diagonal elements
    do while (H(2,2) <= H(2,1) .or. H(2,1)<0)
       if (H(2,2) <= H(2,1)) then; multiple = 1
       else; multiple = -1; endif
       H(:,1) = H(:,1) - multiple*H(:,2)
       B(:,1) = B(:,1) - multiple*B(:,2)
    enddo
    do j = 1,2
       do while (H(3,3) <= H(3,j) .or. H(3,j)<0)
          if (H(3,3) <= H(3,j)) then; multiple = 1
          else; multiple = -1; endif
          H(:,j) = H(:,j) - multiple*H(:,3)
          B(:,j) = B(:,j) - multiple*B(:,3)
       enddo
    enddo
    if (any(matmul(S,B)/=H)) stop "END: Transformation matrices didn't work"
    check = reshape(H,(/9/))
    if (any(check((/4,7,8/))/=0)) stop "Not lower triangular"
    if (any(check((/2,3,6/))<0)) stop "Negative elements in lower triangle"
    if (check(2) > check(5) .or. check(3) > check(9) .or. check(6) > check(9)) stop "Lower triangular elements bigger than diagonal"
  ENDSUBROUTINE HermiteNormalForm

  ! !!<summary>Support routines for the Smith and Hermite normal form finders.</summary>
  ! !!<parameter name="A" regular="true"></parameter>
  ! !!<parameter name="M" regular="true"></parameter>
  ! !!<parameter name="k" regular="true"></parameter>
  ! subroutine swap_row(A,M,k) ! Swap rows of M (and A)
  !   integer, intent(inout) :: M(3,3), A(3,3)
  !   integer, intent(in) :: k
  !   integer :: tmpRow(3), maxidx(1)

  !   maxidx = maxloc(abs(M(k:,k)))+k-1  ! find index of the non-zero element in col k
  !   tmpRow = A(k,:); A(k,:) = A(maxidx(1),:); A(maxidx(1),:) = tmpRow
  !   tmpRow = M(k,:); M(k,:) = M(maxidx(1),:); M(maxidx(1),:) = tmpRow
  ! endsubroutine swap_row

  !!<summary>Swaps the column 'k' with whichever column has the
  !!highest value (out of the columns to the right of 'k' in row
  !!'k'). The swap is performed in both matrices 'M' and
  !!'B'.</summary>
  !!<parameter name="M" regular="true">The matrix being
  !!transformed.</parameter>
  !!<parameter name="B" regular="true">Usually a matrix to keep track
  !!of the transformation steps on 'M'.</parameter>
  !!<parameter name="k" regular="true">The column to swap, as
  !!described in summary.</parameter>
  subroutine swap_column(M,B,k)
    integer, intent(inout) :: M(3,3), B(3,3)
    integer, intent(in) :: k
    integer :: tmpCol(3), maxidx(1)

    maxidx = maxloc(abs(M(k,k:)))+k-1 ! find index of the non-zero element in row k
    tmpCol = B(:,k); B(:,k) = B(:,maxidx(1)); B(:,maxidx(1)) = tmpCol
    tmpCol = M(:,k); M(:,k) = M(:,maxidx(1)); M(:,maxidx(1)) = tmpCol
  endsubroutine swap_column

  ! !!<summary>Finds the indices corresponding the maximum and second
  ! !!maximum Values in an integer vector.</summary>
  ! !!<parameter name="invec" regular="true">Input vector</parameter>
  ! !!<parameter name="max" regular="true">Location of maximum
  ! !!value</parameter>
  ! !!<parameter name="smax" regular="true">Location of second maximum
  ! !!value</parameter>
  ! subroutine get_max2max_indices(invec,smax,max)
  !   integer, intent(in) :: invec(3)
  !   integer, intent(out) :: max, smax

  !   integer :: tmpmax(1), tmpmax2(1), vec(3)
  !   vec = abs(invec)
  !   ! Search from the right for the max, this prevents inifinite
  !   ! looping in the SNF routine.
  !   tmpmax = 4 - maxloc(vec(3:1:-1),vec(3:1:-1)>0)
  !   vec(tmpmax(1)) = 0
  !   tmpmax2 = maxloc(vec, vec>0)!4 - maxloc(vec(3:1:-1),vec(3:1:-1)>0)
  !   smax = tmpmax2(1)
  !   max = tmpmax(1)
  ! endsubroutine get_max2max_indices

  !!<summary>Finds the indices corresponding the minimum and maximum
  !!values in an integer vector.</summary>
  !!<parameter name="invec" regular="true"></parameter>
  !!<parameter name="min" regular="true"></parameter>
  !!<parameter name="max" regular="true"></parameter>
  subroutine get_minmax_indices(invec,min,max)
    integer, intent(in) :: invec(3)
    integer, intent(out) :: min, max

    integer :: tmpmin(1), tmpmax(1), vec(3)
    vec = abs(invec)
    tmpmin = minloc(vec,vec>0)
    ! Search from the right for the max so it will be different from the minimum
    ! even if the min and max are the same value
    tmpmax = 4 - maxloc(vec(3:1:-1),vec(3:1:-1)>0)
    min = tmpmin(1)
    max = tmpmax(1)
  endsubroutine get_minmax_indices

  !!<summary>This function finds the greatest common denominator of
  !!several integers. This case works for two integers, given as
  !!separate arguments</summary>
  !!<parameter name="x1" ></parameter>
  !!<parameter name="x2" ></parameter>
  function gcd_2ints(x1, x2) result(divisor)
    integer, intent(in) :: x1, x2
    integer :: divisor

    integer :: a, b
    a = abs(x1); b = abs(x2) ! Make sure inputs are positive
    if (b>a) call swap(a,b)

    do ! Keep dividing a by b, until one of them is zero
       if (b>a) call swap(a,b) ! Keep the bigger number in a's place
       if (b == 0) exit ! we're done when b == 0
       a = mod(a,b) ! Divide a by b and keep only the remainder
    enddo
    divisor = a

  contains
    !!<summary>Swaps x and y.</summary>
    !!<parameter name="x" regular="true"></parameter>
    !!<parameter name="y" regular="true"></parameter>
    subroutine swap(x,y) ! Swap two values
      integer :: x,y,tmp
      tmp = x; x = y; y = tmp
    endsubroutine swap
  end function gcd_2ints

  !!<summary>This function finds the greatest common denominator of
  !!several integers.This case works on a list of integers (a rank-1
  !!array). </summary>
  !!<parameter name="x" regular="true"></parameter>
  function gcd_rank1(x) result(divisor)
    integer, intent(in) :: x(:)
    integer :: divisor

    integer :: a(size(x)), N, indx(1), big2

    N = size(x); a = abs(x)
    if (any(a<0)) stop "GCD requires non-negative integers"
    do ! Divide the biggest number by the second biggest until
       ! the second biggest is zero
       indx = maxloc(a)  ! Find the location of the biggest number
       if (all(a == a(indx(1)))) then ! check if all numbers are the same
          big2 = a(indx(1))
       else   ! The "real" around 'a' is a workaround for a problem in the Absoft compiler
          big2 = int(maxval(real(a),mask=a < a(indx(1)))) ! Find the size of the 2nd biggest number
       endif
       if (big2 == 0) exit
       a(indx(1)) = mod(a(indx(1)),big2)
    enddo
    divisor = a(indx(1)) ! The divisor is the number left when every other member==0
  endfunction gcd_rank1

  !!<summary>This function finds the greatest common denominator of
  !!several integers. This case works on 3 integers, not in an array.</summary>
  !!<parameter name="x1" ></parameter>
  !!<parameter name="x2" ></parameter>
  !!<parameter name="x3" ></parameter>
  function gcd_3ints(x1,x2,x3)
    integer, intent(in) :: x1,x2,x3
    integer :: gcd_3ints
    gcd_3ints = gcd_rank1((/x1,x2,x3/))
  end function gcd_3ints

  !!<summary>This function finds the greatest common denominator of
  !!several integers. This case works on 4 integers, not in an array.</summary>
  !!<parameter name="x1" ></parameter>
  !!<parameter name="x2" ></parameter>
  !!<parameter name="x3" ></parameter>
  !!<parameter name="x4" ></parameter>
  function gcd_4ints(x1,x2,x3,x4)
    integer, intent(in) :: x1,x2,x3,x4
    integer :: gcd_4ints
    gcd_4ints = gcd_rank1((/x1,x2,x3,x4/))
  end function gcd_4ints

  !!<summary>This function checks to see if there is a rational number
  !!in a range (assume a 1-D range. This isn't general for ternaries,
  !!etc.)</summary>
  !!<parameter name="range" regular="true"></parameter>
  !!<parameter name="n" regular="true"></parameter>
  !!<parameter name="eps_" regular="true"></parameter>
  function is_a_rational_in_range(range,n,eps_)
    integer :: range(3)
    integer :: n
    real(dp), optional :: eps_
    logical :: is_a_rational_in_range
    integer :: j2
    real(dp) :: eps,left,right

    is_a_rational_in_range = .false.

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1e-12_dp
    endif
    left  = real(range(1),dp)/range(3) - eps
    right = real(range(2),dp)/range(3) + eps
    do j2 = 0, n ! numerators
       if      (real(j2,dp)/n > left &
            .and. real(j2,dp)/n < right) then ! there is a rational in the range
          is_a_rational_in_range = .true.
          return
       endif
    enddo
    ! If the loops end, then the if condition was never met and
    ! there isn't a rational number in the range
  endfunction is_a_rational_in_range

  !!<summary>This routine generates all of the rational numbers of a
  !!given numerator that lie within a specified range. (assume a 1-D
  !!range. This isn't general for ternaries, etc.)</summary>
  !!<parameter name="range" regular="true"></parameter>
  !!<parameter name="n" regular="true"></parameter>
  !!<parameter name="numerators"></parameter>
  !!<parameter name="eps_" regular="true"></parameter>
  subroutine get_rationals_in_range(range,n,numerators,eps_)
    integer, intent(in) :: range(3)
    integer, intent(in) :: n
    integer, pointer :: numerators(:)
    real(dp), optional, intent(in) :: eps_
    integer :: j, cR
    real(dp) :: eps,left,right

    nullify(numerators)

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1e-12_dp
    endif
    left  = real(range(1),dp)/range(3) - eps
    right = real(range(2),dp)/range(3) + eps
    allocate(numerators(n+1))
    numerators = 0

    cR = 0
    do j = 0, n ! numerators
       if      (real(j,dp)/n > left &
            .and. real(j,dp)/n < right) then ! there is a rational in the range
          cR = cR + 1
          numerators(cR) = j
       endif
    enddo
    if (cR==0) stop "get_rationals_in_range didn't find any in the range"
    numerators => ralloc(numerators,cR)
  end subroutine get_rationals_in_range

END MODULE rational_mathematics
