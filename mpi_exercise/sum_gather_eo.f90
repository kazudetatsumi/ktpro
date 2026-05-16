program sum_by_two_parallel_procs
  use mpi 
  implicit none
  integer :: nprocs, myrank, ierr
  integer :: istat(MPI_STATUS_SIZE)
  integer :: i, iini, ifin
  integer,parameter :: nmax = 100
  integer suml, sumt
  integer, allocatable :: sumv(:)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  iini = 1 + myrank
  ifin = 99 + myrank
  print *, "iini, ifin:", iini, ifin
  suml = 0
  if (myrank == 0) allocate(sumv(nprocs))
  do i = iini, ifin, 2
    suml = suml + i 
  enddo
  print *, "suml:", suml
  call MPI_GATHER(suml, 1, MPI_INTEGER, sumv(1), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (myrank == 0) print *, "sum:", sum(sumv), "from sumv:", sumv

  call MPI_FINALIZE(ierr)
  stop
  end  program  sum_by_two_parallel_procs
