!fortran 90 test program for a kernel density estimation with grobally-optimized bandwidth.
!Kazuyoshi TATSUMI 2021/12/08

module sskernel
  implicit none
  integer, parameter :: xsize=107
  integer, parameter :: tinsize=20
contains

  subroutine ssk
    implicit none 
    double precision x(xsize), tin(tinsize), thist(tinsize+1)
    integer yhist(tinsize)
	data x/4.37, 3.87, 4.00, 4.03, 3.50, 4.08, 2.25, 4.70, 1.73, 4.93, 1.73, 4.62, 3.43, 4.25, 1.68, 3.92, 3.68, 3.10, 4.03, 1.77,&
	4.08, 1.75, 3.20, 1.85, 4.62, 1.97, 4.50, 3.92, 4.35, 2.33, 3.83, 1.88, 4.60, 1.80, 4.73, 1.77, 4.57, 1.85, 3.52, 4.00, 3.70,&
	3.72, 4.25, 3.58, 3.80, 3.77, 3.75, 2.50, 4.50, 4.10, 3.70, 3.80, 3.43, 4.00, 2.27, 4.40, 4.05, 4.25, 3.33, 2.00, 4.33, 2.93,&
	4.58, 1.90, 3.58, 3.73, 3.73, 1.82, 4.63, 3.50, 4.00, 3.67, 1.67, 4.60, 1.67, 4.00, 1.80, 4.42, 1.90, 4.63, 2.93, 3.50, 1.97,&
	4.28, 1.83, 4.13, 1.83, 4.65, 4.20, 3.93, 4.33, 1.83, 4.53, 2.03, 4.18, 4.43, 4.07, 4.13, 3.95, 4.10, 2.72, 4.58, 1.90, 4.50,&
	1.95, 4.83, 4.12/
	double precision dt
	integer i
	dt = (maxval(x) - minval(x))/tinsize
	tin = (/((i*dt+minval(x)), i=1,tinsize)/)
	thist(1:tinsize)=tin(:)
	thist(tinsize+1)=tin(tinsize)+dt
	thist = thist - dt/2
	yhist=hist(x, thist)
	call plothist(yhist)
  end subroutine ssk

  function hist(x, th)
	implicit none
	double precision, intent(in) :: x(xsize), th(tinsize+1)
	integer :: hist(tinsize)
	integer ix, it
	hist(:) = 0
	do ix = 1,xsize
	   do it = 1,xsize
	      if ( x(ix) >= th(it) .and. x(ix) < th(it+1) ) then
			   hist(it) = hist(it) + 1
		  end if
	   end do
	end do
  end function hist
	
  subroutine plothist(yhist)
	integer, intent(in) ::  yhist(tinsize)
	integer ix, j
	do ix = 1, tinsize
      print *, ("*", j=1,yhist(ix))
	end do
  end subroutine plothist
	

end module sskernel


program main
  use sskernel
  call ssk
  stop
  end

