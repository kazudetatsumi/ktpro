!fortran 90 program for a kernel density estimation with localy optimized bandwidths.
!openmpi version: M and tinsize must be multiples of psize. 
!Kazuyoshi TATSUMI 2022/04/08

program testmpi
  include 'mpif.h'
  double precision, dimension(16) :: A
  double precision :: b(8)
  integer :: rank, psize, ierr
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  A = 0.
  do i = 1, 8
    b(i) = i + rank*32  
  enddo
  call mpi_allgather(b, 8, mpi_double_precision, A, 8, mpi_double_precision, mpi_comm_world, ierr)
  if (rank == 0) then
      print *, A
  endif
  stop
  end



