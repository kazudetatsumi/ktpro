program mpi
implicit none
include "mpif.h"

integer :: ierr,myid,nprocs,status(mpi_status_size),i,j,k,&
&          starts(3),newsize(3),oldsize(3)
real    :: x(1:5,1:5,1:5),y(1:5,1:5,1:5),z(2:3,2:5,4:5)
integer :: arr

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)

if(myid == 0) then
  x = 0.0
  call random_number(x)
  starts = (/1,1,3/)
  newsize = (/2,4,2/)
  oldsize = (/5,5,5/)
  call mpi_type_create_subarray(3, oldsize, newsize, starts, mpi_order_fortran, mpi_real, arr, ierr)
  call mpi_type_commit(arr, ierr)
  call mpi_send(x, 1, arr, 1, 1, mpi_comm_world, ierr)
  do i = 2, 3
    do j = 2, 5
      do k = 4, 5
        print*, '#1', x(i, j, k)
      enddo
    enddo
  enddo
  print*,' '
else
  y = 0.0
  call mpi_recv(z, 16, mpi_real, 0, 1, mpi_comm_world, status, ierr)
  print *, z
  !do i = 2, 3
  !  do j = 2, 5
  !    do k = 4, 5
  !      print*, '#2', z(i, j, k)
  !    enddo
  !  enddo
  !enddo
endif

call mpi_finalize(ierr)

stop
end
