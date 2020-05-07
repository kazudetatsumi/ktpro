!subroutine hello(comm) bind(C, name="hello")
!  use mpi
!  implicit none
!  integer(4) :: comm, PETOT, my_rank, ierr
!
!  call MPI_COMM_SIZE(comm, PETOT, ierr)
!  call MPI_COMM_RANK(comm, my_rank, ierr)
!  print *, 'Hello', my_rank, "of", PETOT
!end subroutine hello
function hello(comm) bind(C)
  use mpi
  implicit none
  integer :: comm, PETOT, my_rank, ierr
  integer :: hello

  call MPI_COMM_SIZE(comm, PETOT, ierr)
  call MPI_COMM_RANK(comm, my_rank, ierr)
  print *, 'Hello', my_rank, "of", PETOT, "from fortran hello.so"
  if (my_rank == 0) hello = PETOT
  !call mpi_finalize(ierr)
end function hello
