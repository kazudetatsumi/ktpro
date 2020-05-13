program sum_by_two_parallel_procs
  use mpi 
  implicit none
  integer :: nprocs, myrank, ierr
  integer :: istat(MPI_STATUS_SIZE)
  integer :: i, iini, ifin
  integer,parameter :: nmax = 100
  integer suml, sumt

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  iini = 1 + myrank*nmax/nprocs
  ifin = iini + nmax/nprocs - 1
  print *, "iini, ifin:", iini, ifin
  suml = 0
  do i = iini, ifin
    suml = suml + i 
  enddo
  print *, "suml:", suml
  call MPI_REDUCE(suml, sumt, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (myrank == 0) print *, "sum:", sumt

  call MPI_FINALIZE(ierr)
  stop
  end  program  sum_by_two_parallel_procs
