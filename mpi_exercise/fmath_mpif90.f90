subroutine scal(x, alpha, N, comm) bind(C)
implicit none
integer(8), intent(in) :: N
real(8), dimension(N), intent(inout) :: x 
real(8), intent(in) :: alpha
real(8), dimension(N) :: y
include "mpif.h"
integer(4) :: comm, size, rank, ierr, len 
integer(4) :: start, end
call MPI_Comm_size(comm, size, ierr) 
call MPI_Comm_rank(comm, rank, ierr) 
len = N/size
start = rank * len + 1
end = start + len - 1
x(start:end) = x(start:end) * alpha
call MPI_Allgather(x(start), len, MPI_REAL8, y, len, MPI_REAL8, comm, ierr)
x(:) = y(:)
end subroutine scal
