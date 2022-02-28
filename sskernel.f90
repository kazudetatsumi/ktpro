!fortran 90 test program for a kernel density estimation with grobally-optimized bandwidth.
!Kazuyoshi TATSUMI 2021/12/08

module sskernel
  use ISO_C_BINDING
  implicit none
  include 'fftw3.f'
  integer :: xsize
  integer :: tinsize
  double precision, parameter :: pi  = 4 * atan (1.0_8)
contains

  subroutine ssk(optw, xsize0, tinsize0, xdat, tin, yh) bind(C, name="ssk")
	real(c_double), intent(out) :: optw
	integer, intent(in) :: xsize0, tinsize0
	double precision, intent(in) :: xdat(xsize0)
	double precision, intent(in) :: tin(tinsize0)
	double precision, intent(inout) :: yh(tinsize0)
    double precision :: xdatstd(xsize0), xdatstddiff(xsize0-1), xdatstddiffstd(xsize0-1)
    double precision :: thist(tinsize0+1), y_hist(tinsize0), y(tinsize0)
    double precision T, dt_samp, dt, cost, w
    integer :: yhist(tinsize0)
    integer nbin, nsmpl, i, tinsize2
	tinsize=tinsize0
	xsize=xsize0
	!data xdat/4.37, 3.87, 4.00, 4.03, 3.50, 4.08, 2.25, 4.70, 1.73, 4.93, 1.73, 4.62, 3.43, 4.25, 1.68, 3.92, 3.68, 3.10, 4.03, 1.77,&
	!4.08, 1.75, 3.20, 1.85, 4.62, 1.97, 4.50, 3.92, 4.35, 2.33, 3.83, 1.88, 4.60, 1.80, 4.73, 1.77, 4.57, 1.85, 3.52, 4.00, 3.70,&
	!3.72, 4.25, 3.58, 3.80, 3.77, 3.75, 2.50, 4.50, 4.10, 3.70, 3.80, 3.43, 4.00, 2.27, 4.40, 4.05, 4.25, 3.33, 2.00, 4.33, 2.93,&
	!4.58, 1.90, 3.58, 3.73, 3.73, 1.82, 4.63, 3.50, 4.00, 3.67, 1.67, 4.60, 1.67, 4.00, 1.80, 4.42, 1.90, 4.63, 2.93, 3.50, 1.97,&
	!4.28, 1.83, 4.13, 1.83, 4.65, 4.20, 3.93, 4.33, 1.83, 4.53, 2.03, 4.18, 4.43, 4.07, 4.13, 3.95, 4.10, 2.72, 4.58, 1.90, 4.50,&
	!1.95, 4.83, 4.12/
	T=maxval(xdat)-minval(xdat)
	!xdatstd=xdat
	!call quicksort(xdatstd, 1, xsize)
	!xdatstddiff=xdatstd(2:xsize)-xdatstd(1:xsize-1)
	!xdatstddiffstd=xdatstddiff
	!call quicksort(xdatstddiffstd, 1, size(xdat)-1)
	!dt_samp=minval(pack(xdatstddiffstd, xdatstddiffstd > 0.))
	!if (ceiling(T/dt_samp) > 1e3) then
	!	tinsize2 = 1e3
	!else
	!	tinsize2 = ceiling(T/dt_samp)
	!endif
	!print *, 'tinsize2=',tinsize2
	!allocate(tin(tinsize), thist(tinsize+1), y_hist(tinsize), yh(tinsize), yhist(tinsize), y(tinsize))
	dt=T/(tinsize-1)
	!tin = (/(((i-1)*dt+minval(xdat)), i=1,tinsize)/)
	thist(1:tinsize)=tin(:)
	thist(tinsize+1)=tin(tinsize)+dt
	thist = thist - dt/2
	yhist=hist(xdat, thist)
	nsmpl=sum(yhist)
	y_hist=real(yhist)/real(nsmpl)/dt
	call opt(optw, yh, y_hist, xdat, nsmpl, dt)
	write(*, '(f12.5)') optw
  end subroutine ssk

  function hist(x, th)
	double precision, intent(in) :: x(xsize), th(tinsize+1)
	integer :: hist(tinsize)
	integer ix, it
	hist(:) = 0
	do ix = 1, xsize
	   do it = 1, tinsize
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

  subroutine opt(optw, y, y_hist, xdat, nsmpl, dt)
	double precision, intent(in) :: y_hist(tinsize), dt, xdat(xsize)
	integer, intent(in) :: nsmpl
	integer, parameter :: maxiter = 20
	double precision, parameter :: tol = 10e-5
	double precision, parameter :: phi = (5**0.5 + 1) / 2
	double precision :: cost(maxiter), Win(maxiter), dummy(tinsize), yh1(tinsize), yh2(tinsize)
	double precision, intent(out) :: y(tinsize), optw
	double precision :: Winmin, Winmax, a, b, c1, c2, f1, f2
	integer :: kiter
	kiter=1
	cost=0.
	Win=0.
	Winmin = 2*dt
	Winmax = maxval(xdat) - minval(xdat)
	a=ilogexp(Winmin)
	b=ilogexp(Winmax)
	c1=(phi - 1)*a + (2 - phi)*b
	c2=(2 - phi)*a + (phi - 1)*b
	call costfunction(dummy, f1, y_hist, nsmpl, logexp(c1), dt)
	call costfunction(dummy, f2, y_hist, nsmpl, logexp(c2), dt)
	do while ( (abs(a-b) > tol*(abs(c1)+abs(c2))) .and. (kiter <= maxiter) )
	  if (f1 < f2) then
	     b=c2
		 c2=c1
		 c1=(phi-1)*a + (2-phi)*b
		 f2=f1
		 call costfunction(yh1, f1, y_hist, nsmpl, logexp(c1), dt)
		 Win(kiter)=logexp(c1)
		 cost(kiter)=f1
		 optw=logexp(c1)
		 y=yh1/sum(yh1*dt)
	  else
		 a=c1
		 c1=c2
		 c2=(2-phi)*a + (phi-1)*b
		 f1=f2
		 call costfunction(yh2, f2, y_hist, nsmpl, logexp(c2), dt)
		 Win(kiter)=logexp(c2)
		 cost(kiter)=f2
		 optw=logexp(c2)
		 y=yh2/sum(yh2*dt)
	  endif
	  print *, kiter, cost(kiter), optw
	  kiter=kiter+1
	enddo
  end subroutine opt

  subroutine costfunction(yh, cost, y_hist, nsmpl, w, dt)
	integer, intent(in) :: nsmpl
	double precision, intent(in) ::  y_hist(tinsize), w, dt
	double precision, intent(out) :: yh(tinsize), cost
	yh=0.
	yh=fftkernel(y_hist, w/dt)
	cost=sum(yh**2)*dt - 2*sum(yh*y_hist)*dt + 2/(2*pi)**0.5/w/nsmpl
	cost=cost*nsmpl**2
  end subroutine costfunction

  function fftkernel(x, w)
	!include 'fftw3.f'
	double precision, intent(in) :: x(tinsize), w
	double precision :: fftkernel(tinsize), Lmax
	integer :: L, N, i
	integer(8) :: plan
	double complex, allocatable :: input(:), output(:), xinput(:), Xoutput(:), Youtput(:), K(:)
	double precision, allocatable :: f(:), expa(:)
	L=tinsize
	Lmax=L+3*w
	N = 2**ceiling(log(Lmax)/log(2.))
	allocate(input(N), output(N), xinput(N), Xoutput(N), f(N), expa(N), K(N), Youtput(N))
	call dfftw_plan_dft_1d(plan, N, input, output, FFTW_FORWARD, FFTW_ESTIMATE)
	xinput(1:tinsize)=x
	xinput(tinsize+1:N)=0.0
	call dfftw_execute_dft(plan, xinput, Xoutput)
	f(1:N/2+1) = -1.0*real((/(i, i=0, N/2, 1)/))/real(N)
	f(N/2+2:N) = 0.5 - real((/(i, i=1, N/2 - 1, 1)/))/real(N)
	! for preventing the procedure from underflow exception.
	expa=-0.5*(w*2*pi*f)**2
	K = 0.
	where (expa > -708) K=exp(expa)
	call dfftw_destroy_plan(plan)
	call dfftw_plan_dft_1d(plan, N, input, output, FFTW_BACKWARD, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, Xoutput*K, Youtput)
	call dfftw_destroy_plan(plan)
	! Computing a fftw forward followed by a fftw backward results in the original
    ! array scaled by N. So, we should divid the result by N so as to obtain the
	! count dependent cost fucntion.
	fftkernel=real(Youtput(1:L)/N)
  end function fftkernel

  function ilogexp(x)
	double precision, intent(in) :: x
	double precision :: ilogexp
	if ( x < 1e2 ) then
       ilogexp = log(exp(x) - 1)
	else
	   ilogexp = x
	endif 
  end function ilogexp

  function logexp(x)
	double precision, intent(in) :: x
	double precision :: logexp
	if ( x < 1e2 ) then
       logexp = log(exp(x) + 1)
	else
	   logexp = x
	endif 
  end function logexp

! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
  recursive subroutine quicksort(a, first, last)
    implicit none
    real*8  a(*), x, t
    integer first, last
    integer i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

end module sskernel


!program main
!  use sskernel
!  call ssk
!  stop
!  end

