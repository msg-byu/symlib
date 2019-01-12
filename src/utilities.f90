!!<summary> file : utilities.f90 function : Collection of utility
!!subroutines.</summary>

module utilities
  use num_types
  private
  public ucase, ralloc

  INTERFACE ralloc
     MODULE PROCEDURE ralloc_real, ralloc_integer, ralloc_integer_matrix_list, &
          ralloc_integer_table, ralloc_real_table
  END INTERFACE ralloc

contains
       
  !!<summary>convert string to UPPERCASE</summary>
  !!<parameter name="string" regular="true"></parameter>
  subroutine ucase(string)
    implicit none
    character(len=*), intent(INOUT) :: string
    integer :: i
    integer :: length

    length=len(string)
    do i=1,length
       if (lge(string(i:i),'a').and.(lle(string(i:i),'z'))) then
          string(i:i)=achar(iachar(string(i:i))-32)
       endif
    enddo
  end subroutine ucase

  !!<summary>Reallocates an array.</summary>
  !!<parameter name="p">Array to be extended</parameter>
  !!<parameter name="n" regular="true">Length of new array</parameter>
  FUNCTION ralloc_real(p, n)              

    REAL(dp), POINTER, DIMENSION(:) :: p, ralloc_real
    INTEGER, intent(in) :: n
    INTEGER :: nold, ierr
    
    ALLOCATE(ralloc_real(1:n), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    IF(.NOT. ASSOCIATED(p)) RETURN
    nold = MIN(SIZE(p), n)
    ralloc_real = 0._dp
    ralloc_real(1:nold) = p(1:nold)
    DEALLOCATE(p) 

  END FUNCTION RALLOC_REAL

  !!<summary>Reallocates an array.</summary>
  !!<parameter name="p">Array to be extended</parameter>
  !!<parameter name="n" regular="true">Length of new array</parameter>
  FUNCTION ralloc_integer(p, n)               

    integer, POINTER, DIMENSION(:) :: p, ralloc_integer
    INTEGER, intent(in) :: n
    INTEGER :: nold, ierr

    ALLOCATE(ralloc_integer(1:n), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    IF(.NOT. ASSOCIATED(p)) RETURN
    nold = MIN(SIZE(p), n)
    ralloc_integer = 0
    ralloc_integer(1:nold) = p(1:nold)
    DEALLOCATE(p) 

  END FUNCTION RALLOC_INTEGER

  !!<summary>Reallocates an array. GLWH Added 10/2010</summary>
  !!<parameter name="p">Array to be extended</parameter>
  !!<parameter name="n" regular="true">Length of new array</parameter>
  FUNCTION ralloc_integer_matrix_list(p, n)          
    
    integer, POINTER, DIMENSION(:,:,:) :: p, ralloc_integer_matrix_list
    INTEGER, intent(in) :: n
    INTEGER :: nold, ierr, r, c

    r = size(p,1); c = size(p,2)
    ALLOCATE(ralloc_integer_matrix_list(1:r,1:c,1:n), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error in reallocate.f90 (integer_matrix)"
    IF(.NOT. ASSOCIATED(p)) RETURN
    nold = MIN(SIZE(p,3), n) ! Number of rows
    ralloc_integer_matrix_list = 0
    ralloc_integer_matrix_list(1:r,1:c,1:nold) = p(1:r,1:c,1:nold)
    DEALLOCATE(p) 

  END FUNCTION RALLOC_INTEGER_MATRIX_LIST

  !!<summary>Reallocates an array.</summary>
  !!<parameter name="p">Table to be extended></parameter>
  !!<parameter name="n" regular="true">Rom size of new table.</parameter>
  !!<parameter name="m" regular="true">Column size of new table</parameter>
  FUNCTION ralloc_integer_table(p, n, m)
    integer, pointer, dimension(:,:) :: p, ralloc_integer_table
    integer, intent(in) :: n, m ! row and column size of new table
    integer :: nold, mold, ierr
    allocate(ralloc_integer_table(1:n,1:m),STAT=ierr)
    if (ierr/=0) stop "Allocate error in ralloc_integer_table in utilities module"
    if (.not.associated(p)) return
    nold = min(size(p,1),n); mold = min(size(p,2),m)
    ralloc_integer_table = 0
    ralloc_integer_table(1:nold,1:mold) = p(1:nold,1:mold)
    deallocate(p)
  END FUNCTION ralloc_integer_table

  !!<summary>Reallocates an array.</summary>
  !!<parameter name="p">Table to be extended></parameter>
  !!<parameter name="n" regular="true">Rom size of new table.</parameter>
  !!<parameter name="m" regular="true">Column size of new table</parameter>
  FUNCTION ralloc_real_table(p, n, m)
    real(dp), pointer, dimension(:,:) :: p, ralloc_real_table
    integer, intent(in) :: n, m ! row and column size of new table
    integer :: nold, mold, ierr
    allocate(ralloc_real_table(1:n,1:m),STAT=ierr)
    if (ierr/=0) stop "Allocate error in ralloc_integer_table in utilities module"
    if (.not.associated(p)) return
    nold = min(size(p,1),n); mold = min(size(p,2),m)
    ralloc_real_table = 0
    ralloc_real_table(1:nold,1:mold) = p(1:nold,1:mold)
    deallocate(p)
  END FUNCTION ralloc_real_table

end module utilities
