program gettotalintensity
  implicit none
  integer, parameter :: nlines=41352500
  double precision :: intensity(nlines), qx(nlines), qy(nlines), qz(nlines), e(nlines), sig(nlines)
  double precision :: xi, xe, dx
  character(len=80) :: argc
  character(len=80) :: dummy
  integer, parameter :: ndum=9
  integer :: idum, iint, ios

  call getarg(1, argc)
  !print *, trim(argc)
  open(10, file=trim(argc), status='old')
  do idum=1, ndum
     read(10,'(a80)') dummy
     !print *, dummy
  enddo
  do iint=1, nlines
     read(10,*, end=99) qx(iint),qy(iint),qz(iint),e(iint),intensity(iint),sig(iint)
     !print *, qx(iint),qy(iint),qz(iint),e(iint),intensity(iint),sig(iint)
  enddo
  99 print *, "totlal intensity of ", trim(argc)," is ",sum(intensity)
  !print *, "finish reading"
  !print *, "totlal intensity of ", trim(argc)," is ",sum(intensity)
  close(10) 


end program gettotalintensity
