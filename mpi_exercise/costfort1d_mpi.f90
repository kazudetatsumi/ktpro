module costfort
  use ISO_C_binding
  use mpi
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    type(c_ptr) :: arr
    type(c_ptr) :: kavearr
    type(c_ptr) :: darr
  end type result

contains

  function cost1d(maxw, Al0, A) bind(C, name="cost1d")
    integer :: nprocs, myrank, ierr
    integer :: istat(MPI_STATUS_SIZE)
    integer :: iini, ifin
    integer(c_int), intent(in) :: maxw
    integer(c_int), intent(in) :: Al0
    real(c_double), intent(in) :: A(Al0)                         
    type(result) :: cost1d                                  
    real(c_double), pointer :: k(:)                    
    real(c_double), pointer :: Cn(:)                    
    real(c_double), pointer :: kaves(:)                    
    real(c_double), pointer :: deltas(:)                    
    integer nw0
    integer N0
    integer i, ihead
    real kave,v
    allocate(k(Al0))
    allocate(Cn(maxw))
    allocate(kaves(maxw))
    allocate(deltas(maxw))
    k(:) = 0.0
    Cn(:) = 0.0
    kaves(:) = 0.0
    deltas(:) = 0.0
    
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if ( mod(maxw, nprocs) == 0 ) then
       iini = 1 + myrank*maxw/nprocs
       ifin = iini + maxw/nprocs - 1
    else
       iini = 1 + myrank*(maxw-mod(maxw, nprocs))/nprocs
       ifin = iini + (maxw-mod(maxw, nprocs))/nprocs - 1
       if (myrank < mod(maxw, nprocs)) ifin = ifin + 1
    endif

    
    !do nw0 = 1, maxw
    do nw0 = iini, ifin
       N0 = (Al0 - mod(Al0, nw0)) / nw0
       do i = 1, N0
          ihead = i*nw0 
          if ( i == 1) then
              k(i) = A(ihead)
          else
              k(i) = A(ihead) - A(ihead - nw0)
          end if
       end do
       kave = sum(k(1:N0)) / real(N0)
       v = sum((k(1:N0) - kave)**2) / real(N0)
       Cn(nw0) = (2.0 * kave - v) / (real(nw0)**2)
       deltas(nw0) = real(nw0)
       kaves(nw0) = kave
    end do
    call MPI_GATHER(Cn(iini), ifin - iini + 1, MPI_DOUBLE_PRECISION,                  &
                    Cn, ifin - iini + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(deltas(iini), ifin - iini + 1, MPI_DOUBLE_PRECISION,              &
                    deltas, ifin - iini + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(kaves(iini), ifin - iini + 1, MPI_DOUBLE_PRECISION,               &
                    kaves, ifin - iini + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (myrank == 0 ) then

      print *, 'size of Cn',size(Cn)
      print *, 'len0', nw0
      print *, 'maxw', maxw

      cost1d%len0 =  maxw
      cost1d%arr = C_loc(Cn(1:maxw))
      cost1d%kavearr = C_loc(kaves(1:maxw))
      cost1d%darr = C_loc(deltas(1:maxw))
    endif
    call MPI_FINALIZE(ierr)
  end function cost1d

end module costfort
