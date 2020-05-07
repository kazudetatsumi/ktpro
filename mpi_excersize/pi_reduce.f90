program pi_parallel
  use mpi 
  implicit none
  integer :: nprocs, myrank, ierr
  integer :: i, iini, ifin
  integer, parameter :: SP = kind(2.0)
  integer, parameter :: DP = selected_real_kind(2*precision(1.0_SP))
  integer, parameter :: nmax = 1000000000
  real(DP) :: x, dx, pil, pi

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  iini = 1 + myrank*nmax/nprocs
  ifin = iini + nmax/nprocs - 1

  dx = 1.0_DP/real(nmax,DP)
  pil = 0.0_DP
  do i = iini, ifin
    x = real(i, DP) * dx
    pil = pil + 4.0_DP/(1.0_DP + x ** 2)*dx
  end do

  print *, "pil:", pil
  call MPI_REDUCE(pil, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (myrank == 0) print *, "pi:", pi

  call MPI_FINALIZE(ierr)
  stop
end program pi_parallel

