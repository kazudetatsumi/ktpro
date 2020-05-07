subroutine pi_reduce_mpif90_by_py(comm, nmax, pi) bind(C)
  use mpi 
  use ISO_C_binding
  implicit none
  integer, intent(in) :: nmax
  real(c_double), intent(inout) :: pi
  integer(c_int), intent(in) :: comm
  integer :: nprocs, myrank, ierr
  integer :: i, iini, ifin
  double precision  :: x, dx, pil

  !call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  iini = 1 + myrank*nmax/nprocs
  ifin = iini + nmax/nprocs - 1

  dx = 1/nmax
  pil = 0
  do i = iini, ifin
    x = i * dx
    pil = pil + 4/(1 + x ** 2)*dx
  end do

  print *, "pil:", pil
  call MPI_REDUCE(pil, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  !if (myrank == 0) print *, "pi:", pi

  !call MPI_FINALIZE(ierr)
  !stop
end subroutine pi_reduce_mpif90_by_py

