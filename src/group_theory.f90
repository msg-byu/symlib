MODULE group_theory
use utilities, only : ralloc
implicit none
private
public grouper

CONTAINS

  !!<summary>Generates a group from a list of generators</summary>
  !!<usage>Pass in 2D array where the rows are the permutations. The
  !!permutations will be taken in all possible pairs to generate new
  !!group elements. The array of elements will be added to until the
  !!list of permutations has closure and constitutes a group.</usage>
  !!<comments>Permutations are indexed 1..n, not 0..n-1</comments>
  !!<parameter name="g" regular="true">On entry contains the
  !!generators, on exit contains the whole group </parameter>
  subroutine grouper(g)
    integer, pointer:: g(:,:) 

    integer :: nG ! Number of elements in the group (on exit)
    logical :: growing, new ! Flags: group is still being added to,
    ! new group element was found
    integer :: mG ! Number of elements already checked (don't need to be re-checked)
    integer :: ig, jg, kg, cg ! Generic loop variables
    integer :: n ! Number of permutation elements
    integer,allocatable :: newg(:), vs(:) ! permuted group element, vector subscript

    n = size(g,2) ! Number of permutation elements
    if (n<1) stop "ERROR: Empty list passed to 'grouper' subroutine in 'group_theory' module (celib)"
    allocate(newg(n),vs(n))
    if(any(g>n) .or. any(g<1)) then
       write(*,'("n: ",i4)') n
       write(*,'("g: ",5000(i3,1x))') 
       stop "ERROR: generators passed to 'grouper' (group_theory in celib) contain indices outside of 1..n"
    endif
    growing = .true.
    mG = 0
    nG = size(g,1)
    do while (growing)
       growing = .false.
       g => ralloc(g,nG+nG**2,n) ! Make room for possible new elements
       ! when the existing elements are combined
       ! Loop over each pair of group elements...
       cg = 0 ! Keep track of the number of new elements of the group
       ! that have been found in this iteration
       do ig = 1, nG
          do jg = 1, nG
             if (ig <= mG .and. jg <= mG) cycle !...but skip those
             !already tested in a previous iteration
             ! Construct a permuted element, that is, multiply two group elements
             vs = g(jg,:) ! Make a vector subscript for the permuted element
             newg = g(ig,vs) ! Permute the element
             ! Is this new g unique in the list of existing group elements?
             new = .true.
             newcheck: do kg = 1, nG+cg
                if (all(newg == g(kg,:))) then
                   new = .false.
                   exit newcheck
                endif
             enddo newcheck
             if (new) then
                growing = .true.
                cg = cg + 1
                g(nG+cg,:) = newg
             end if
          end do
       end do ! Loops over pairs
       mG = nG
       nG = nG + cg
    end do ! Loop for group still growing
    g => ralloc(g,nG,n)
  end subroutine grouper
END MODULE group_theory
