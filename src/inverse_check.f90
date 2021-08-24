PROGRAM hello
use num_types
use vector_matrix_utilities
implicit none

real(dp) mat(1000000,3,3)
real(dp) ans(3,3)
integer i, n, nM
integer, allocatable :: state(:)
logical errflag
real(dp) start, stop

nM = 1000000
!nM = 10
call random_seed(size=n)
allocate(state(n))
state = 20210518
call random_seed(put=state)
call random_number(mat)
call cpu_time(start)
do i = 1,nM
    call matrix_inverse(mat(i,:,:),ans,err_=errflag,eps_=1e-8_dp)
    if (errflag) then
        print*,"epsilon triggered stop"
        print*, i
        print*,matmul(mat(i,:,:),ans)
        stop
    endif
enddo
call cpu_time(stop)
print*,"Time/inverse: ",(stop - start)/nM
end
