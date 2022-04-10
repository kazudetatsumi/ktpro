!fortran 90 program for a kernel density estimation with localy optimized bandwidths.
!openmpi version: M and tinsize must be multiples of psize. 
!Kazuyoshi TATSUMI 2022/04/08

program testmpi
  include 'mpif.h'
  double precision, dimension(4, 4) :: A
  double precision :: b(4,2)
  integer :: rank, psize, ierr
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  A = 0.
  do i = 1, 4
  do j = 1, 2
    b(i,j) = i*j + rank*32  
  enddo
  enddo
  call mpi_allgather(b, 8, mpi_double_precision, A, 8, mpi_double_precision, mpi_comm_world, ierr)
  if (rank == 0) then
      print *, A(:,1)
      print *, A(:,2)
      print *, A(:,3)
      print *, A(:,4)
  endif
  call MPI_Finalize(ierr)
  stop
  end



