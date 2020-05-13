subroutine pifort(nmax, comm, pi) bind(C, name="pif")
  use mpi 
  !use ISO_C_binding
  implicit none
  integer, intent(in) :: nmax
  !real(8), intent(inout) :: pi
  double precision, intent(inout) :: pi
  integer(4)  :: comm
  integer :: nprocs, myrank, ierr
  integer :: i, iini, ifin
  double precision :: x, dx, pil, pi_tmp

  !call MPI_INIT(ierr)
  print *, nmax
  call MPI_COMM_SIZE(comm, nprocs, ierr)
  call MPI_COMM_RANK(comm, myrank, ierr)

  iini = 1 + myrank*nmax/nprocs
  ifin = iini + nmax/nprocs - 1

  dx = 1/real(nmax)
  print *, dx
  pil = 0
  do i = iini, ifin
    x = i * dx
    pil = pil + 4/(1 + x ** 2)*dx
  end do

  print *, "pil:", pil
  !call MPI_REDUCE(pil, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(pil, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  print *, pi

  !call MPI_FINALIZE(ierr)
  !stop
end subroutine pifort

