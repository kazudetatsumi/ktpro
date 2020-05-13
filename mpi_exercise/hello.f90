program hello
  use mpi 
  implicit none
  integer :: PETOT, my_rank, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  print *, 'Hello', my_rank, "of", PETOT
  call MPI_FINALIZE(ierr)
  stop
  end
