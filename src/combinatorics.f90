MODULE combinatorics
  use num_types
  use rational_mathematics, only: gcd
  implicit none
  private
  public get_permutations, m_partitions_of_n, factorial, nchoosek, k_ary_counter, &
       multinomial, binomial, permutation_parity

  type subset_mask
     private
     logical          :: init, done
     logical, pointer :: mask(:) => null()
  end type subset_mask

  INTERFACE factorial
     MODULE PROCEDURE factorial_int_scalar, factorial_int_rank1, &
          factorial_dp_rank1, factorial_dp_scalar, factorial_long_int
  END INTERFACE factorial
CONTAINS

  !!<summary>Binomial coefficient. Computes the number of permutations
  !!of n things, taken k at a time</summary>
  !!<comments>! This implementation was taken from "Binomial
  !!Coefﬁcient Computation: Recursion or Iteration?" by Yannis
  !!Manolopoulos, ACM SIGCSE Bulletin InRoads, Vol.34, No.4, December
  !!2002. http://delab.csd.auth.gr/papers/SBI02m.pdf It is supposed to
  !!be robust against large, intermediate values and to have optimal
  !!complexity. This replaces the original function that was used
  !!until revision 125 GLWH June 2014 </comments>
  FUNCTION nchoosek(n,k)
    integer(li) :: nchoosek
    integer, intent(in) :: n,k
    integer(li) :: t ! Variable to accumulate answer in
    integer :: i ! Counter
    if (k < 0 .or. k > n) then
       nchoosek = 0
       return
    end if
    if (k==0 .and. n == 0) then
       nchoosek = 1
       return
    end if
    t = 1
    if (k<n-k) then
       do i = n, n-k+1, -1
          t = t*i/(n-i+1)
       end do
    else
       do i = n, k+1, -1
          t = t*i/(n-i+1)
       end do
    end if
    nchoosek = t
  END FUNCTION nchoosek

  !!<summary>Finds the binomial coefficient of two different numbers.</summary>
  FUNCTION binomial(n,k)
    integer(li) :: binomial
    integer, intent(in) :: n,k
    if (k > n) stop "Error in arguments to binomial"
    binomial = nchoosek(n,k)

  END FUNCTION binomial

  !!<summary>Calculates the multinomial for an input list of
  !!numbers.</summary>
  !!<parameter name="n" regular="true">An array of the exponents of
  !!the term we are taking the multinomial of.</parameter>
  FUNCTION multinomial(n)
    integer(li) :: multinomial
    integer, intent(in) :: n(:)
    integer(li) :: binomials(size(n,1),2), bins(size(n,1))
    integer :: i

    binomials(1,1) = sum(n)
    binomials(1,2) = n(1)
    do i = 2, size(binomials,1)
       binomials(i,1) = binomials(i-1,1) - binomials(i-1,2)
       binomials(i,2) = n(i)
    enddo

    do i = 1, size(binomials,1)
       bins(i) = binomial( int(binomials(i,1)) , int(binomials(i,2)) )
    enddo

    multinomial = product(bins)

    !nt = n; ft = 1
    !ft=factorial(sum(nt))
    !do i = 1, size(n)
    !   if (mod(ft,factorial(nt(i)))/=0) stop "Bug in multinomial"
    !   ft = ft/factorial(nt(i))
    !enddo
    !multinomial = ft
  ENDFUNCTION multinomial

  !!<summary>This subroutine generates all permutations of the
  !! list. Swapping identical items does *not* give another
  !! permutation. List is generated in lexicographical order. Input
  !! list must be in order (or the routine quits with an
  !! error)</summary>
  !!<parameter name="list" regular="true">The input list for which
  !! permutations are neede.</parameter>
  !!<parameter name="perms">The output 2d array of possible
  !! permutations.</parameter>
  SUBROUTINE get_permutations(list,perms)
    integer, intent(in) :: list(:)
    integer, pointer :: perms(:,:)

    integer, allocatable:: tempPerms(:,:)
    integer :: Nitems, Np, a(size(list)),i,nuq(size(list))
    integer :: j, l, temp, ip, status
    integer  :: ptr ! "pointer" to current position in the list of unique items in list
    real(dp) :: tot

    Nitems = size(list)
    ! Determine how much temporary storage we will need. Tricky to do this
    ! without overflowing an integer (or even a real)
    call find_unique(list,nuq)
    ! Make a list of terms in the numerator and denominator of the multinomial
    ! and keep a running product
    ptr = 1; tot = 1
    do i = Nitems,1,-1 ! Loop over each term, start at N in the numerator and work down
       if (nuq(ptr)==0) ptr = ptr + 1 ! Keep track of the number in the denominator
       tot = tot*real(i)/nuq(ptr) ! compute multinomial one pair (numer/denom) at a time
       nuq(ptr) = nuq(ptr) - 1
    enddo
    Np = ceiling(tot)
    if (Np < 0) stop "combinatorics.f90: temporary storage overflow. Cutoffs too big."

    allocate(tempPerms(Nitems,Np))
    a = list ! Copy of the original input list
    ! Make sure the items are initially ordered
    do i = 1, Nitems-1
       if (a(i+1)<a(i)) stop "ERROR in 'get_permutations'---incorrectly ordered input"
    enddo
    ip = 0
    outer:do
       ip = ip + 1;
       tempPerms(:,ip)=a
       do j = Nitems-1,0,-1
          if (j==0) exit outer
          if (a(j) < a(j+1)) exit
       enddo
       do l = Nitems,1,-1
          if (a(j) < a(l)) exit
       enddo
       temp = a(j)
       a(j) = a(l)
       a(l) = temp
       a(j+1:Nitems)=a(Nitems:j+1:-1)
    enddo outer
    !if(associated(perms)) deallocate(perms)
    allocate(perms(Nitems,ip),STAT=status)
    if(status/=0) stop "problem allocating 'perms' array in get_permutations in combinatorics.f90 (celib)"
    perms = tempPerms(:,:ip) ! There will be fewer permutations than N! if there
    deallocate(tempPerms)    ! are identical items in the permuted list

  contains
    !!<summary>Finds the number of unique entries of each type in a
    !!list.</summary>
    !!<parameter name="list">A 1D array of for which the unique
    !!entries are to be found.</parameter>
    !!<parameter name="uqlist">A 1D array of the unique elemements of
    !!the input list.</parameter>
    subroutine find_unique(list,uqlist)
      integer, intent(in) :: list(:)
      integer, intent(out):: uqlist(:)
      integer :: sm
      uqlist = 0
      sm = minval(list) ! smallest value in the list
      do i = 1, size(list) ! loop over the possible number of unique values

         uqlist(i) = count(sm==list) ! How many of this size?
         sm = minval(list,mask=(list>sm))    ! Next unique value bigger than the last
         if (count(sm==list)==0) exit ! If there aren't any of this
         ! size then we took the largest already
      enddo
    endsubroutine find_unique
  END SUBROUTINE get_permutations

  !!<summary>This routine generates a list of all numbers between 0
  !!and n^k-1</summary>
  !!<parameter name="base" regular="true">Base is the intereger that
  !!will be used as the base of the exponential.</parameter>
  !!<parameter name="n" regular="true">n is an integer that is being
  !!exponentiated.</parameter>
  !!<parameter name="list">Resultant 2d list of numbers</parameter>
  SUBROUTINE k_ary_counter(list,base,n)
    integer, pointer :: list(:,:)
    integer, intent(in) :: base
    integer, intent(in) :: n ! number of digits
    integer :: k, j, il
    integer :: a(n)

    k = base; allocate(list(n,k**n))
    il = 1; a = 0
    list(:,il) = a
    do
       j = n ! Number of digits
       do ! Check to see if we need to roll over any digits, start at the right
          if (a(j) /= k - 1) exit
          a(j) = 0  ! Rolling over so set to zero
          j = j - 1;       ! Look at the next (going leftways) digit
          if (j < 1) exit  ! If we updated the leftmost digit then we're done
       enddo
       if (j < 1) exit ! We're done counting, exit
       a(j) = a(j) + 1 ! Update the next digit
       il = il + 1
       list(:,il) = a
    enddo

  ENDSUBROUTINE k_ary_counter

  !!<summary>This subroutine generates all partions of n into m
  !! blocks. This is needed to generate all possible concentration
  !! vectors for each structures (then within each structure we need
  !! all permutations of atomic configurations...)  See page 38 in
  !! Knuth V. 4 Fascicle 3.</summary>
  !!<parameter name="n" regular="true">Number to be partitioned.</parameter>
  !!<parameter name="m" regular="true">Number of the block.</parameter>
  !!<parameter name="part" regular="true">Resultant array of partition.</parameter>
  SUBROUTINE m_partitions_of_n(n,m,part)
    integer, intent(in) :: n, m ! number to be partitioned, number of block
    integer, pointer   :: part(:,:)

    integer :: a(m)
    integer          :: ip  ! counter for number of partitions
    integer          :: x, s ! Temporary variables for swapping/reassigning values
    integer          :: j    ! index pointer for next position to be incremented

    ! First we need to count the number of partitions to know how big
    ! to allocate the "part" array (is there a closed form expression
    ! for this? I don't think so but...)
    if (m>n) then; stop "Bad input for 'm_partitions_of_n' routine";endif
    if (m<2) then; stop "Trivial case'm_partitions_of_n' routine";endif
    ! Initialize the first partition
    a(1) = n - (m-1); a(2:m) = 1
    ip = 0
    do; ip = ip + 1;
       if (a(2) < a(1)-1) then ! "crumble" at left if possible
          a(1) = a(1)-1
          a(2) = a(2)+1
          cycle; endif
       if (m==2) exit ! For binary case, we're done here
       ! Find the leftmost position that can be increased
       j = 3 ! This is the leftmost possibility but...
       s = a(1)+a(2)-1
       do while (a(j)>=a(1)-1)
          s = s+a(j)
          j = j+1
          if (j>m) exit
       enddo ! Now s = part(1)+part(2)+...+part(j-1)-1)
       if (j > m) exit ! We're done counting
       x = a(j)+1 ! Store the value of the j-th position incremented by one
       a(j) = x   ! Make the j-th position have this value
       j = j-1    ! Now look one to the left
       do while (j>1)
          a(j) = x  ! Make this next-left position match the slot that was updated
          s = s - x ! Keep track of how many of the original 'n' are left
          j = j-1   ! Now look to the left again...
       enddo
       a(1) = s  ! Now make the leftmost slot take all the leftover...
    enddo
    !if (allocated(part)) deallocate(part)
    allocate(part(ip,m))
    ip = 0; a(1) = n - (m-1); a(2:m) = 1
    do; ip = ip + 1
       ! Visit the partition
       part(ip,:) = a ! store the ip-th partition
       if (a(2) < a(1)-1) then ! "crumble" at left if possible
          a(1) = a(1)-1
          a(2) = a(2)+1
          cycle; endif
       if (m==2) exit ! For binary case, we're done here
       ! Find the leftmost position that can be increased
       j = 3 ! This is the leftmost possibility but...
       s = a(1)+a(2)-1
       do while (a(j)>=a(1)-1)
          s = s+a(j)
          j = j+1
          if (j>m) exit
       enddo ! Now s = part(1)+part(2)+...+part(j-1)-1)
       if (j > m) exit ! We're done counting
       x = a(j)+1
       a(j) = x
       j = j-1
       do while (j>1)
          a(j) = x
          s = s - x
          j = j-1
       enddo
       a(1) = s
    enddo
  END SUBROUTINE m_partitions_of_n

  !!<summary>Factorial of the long int type.</summary>
  FUNCTION factorial_long_int(N)

    integer(li) :: factorial_long_int, N, i
    factorial_long_int = 1
    do i=2,N;factorial_long_int=factorial_long_int*i;enddo
  END FUNCTION factorial_long_int

  !!<summary>Factorial of the int type.</summary>
  FUNCTION factorial_int_scalar(N)
    integer :: N,i,factorial_int_scalar
    factorial_int_scalar = 1
    do i=2,N; factorial_int_scalar=factorial_int_scalar*i;enddo
  END FUNCTION factorial_int_scalar

  !!<summary>Factorial of each element of 1D array of type int.</summary>
  !!<parameter name="N" regular="true"></parameter>
  FUNCTION factorial_int_rank1(N)
    integer :: N(:),i,factorial_int_rank1(size(N)),j
    factorial_int_rank1 = 1
    do j=1,size(N)
       do i=2,N(j); factorial_int_rank1(j)=factorial_int_rank1(j)*i;enddo;enddo
  END FUNCTION factorial_int_rank1

  !!<summary>Factorial of each element 1D array of tpye real.</summary>
  !!<parameter name="N" regular="true"></parameter>
  FUNCTION factorial_dp_rank1(N)
    real(dp) :: N(:),factorial_dp_rank1(size(N))
    integer :: i,j
    factorial_dp_rank1 = 1._dp
    do j=1,size(N)
       do i=2,nint(N(j)); factorial_dp_rank1(j)=factorial_dp_rank1(j)*i;enddo;enddo
  END FUNCTION factorial_dp_rank1


  !!<summary>Factorial of the tpye real.</summary>
  FUNCTION factorial_dp_scalar(N)
    integer :: i
    real(dp):: N,factorial_dp_scalar
    factorial_dp_scalar = 1._dp
    do i=2,nint(N); factorial_dp_scalar=factorial_dp_scalar*i;enddo
  END FUNCTION factorial_dp_scalar


  !!<summary>Finds the parity of a given permuted list aP (must start at 1)
  !! algorithm from: www.cap-lore.com/code/ocaml/parity.html</summary>
  !!<parameter name="aP" regular="true"></parameter>
  function permutation_parity(aP)
    integer :: permutation_parity
    integer, intent(in) :: aP(:)        ! the permuted list
    integer             :: v(size(aP))  ! some helper array
    integer :: j,x,p,n

    n=size(aP)

    ! check if every number i
    do j=1,n
       if (.not. any(aP==j)) stop "ERROR: invalid list in permuation_parity."
    enddo

    v = 0
    p = 0

    do j=n,1,-1
       if (v(j)>0) then
          p=p+1
       else
          x=j
          do
             x = aP(x)
             v(x) = 1
             if (x==j) exit
          enddo
       endif
    end do

    permutation_parity = -(mod(p,2)*2-1)
  end function permutation_parity
END MODULE combinatorics
