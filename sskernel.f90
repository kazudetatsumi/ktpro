!fortran 90 test program for a kernel density estimation with grobally-optimized bandwidth.
!Kazuyoshi TATSUMI 2021/12/08

module sskernel
  implicit none
  integer, parameter :: xsize=107
  integer, parameter :: tinsize=16
contains

  subroutine ssk
    implicit none 
	!include 'fftw3.f'
    double precision x(xsize), tin(tinsize), thist(tinsize+1), y_hist(tinsize), yh(tinsize), cost, w
    integer yhist(tinsize), N
	data x/4.37, 3.87, 4.00, 4.03, 3.50, 4.08, 2.25, 4.70, 1.73, 4.93, 1.73, 4.62, 3.43, 4.25, 1.68, 3.92, 3.68, 3.10, 4.03, 1.77,&
	4.08, 1.75, 3.20, 1.85, 4.62, 1.97, 4.50, 3.92, 4.35, 2.33, 3.83, 1.88, 4.60, 1.80, 4.73, 1.77, 4.57, 1.85, 3.52, 4.00, 3.70,&
	3.72, 4.25, 3.58, 3.80, 3.77, 3.75, 2.50, 4.50, 4.10, 3.70, 3.80, 3.43, 4.00, 2.27, 4.40, 4.05, 4.25, 3.33, 2.00, 4.33, 2.93,&
	4.58, 1.90, 3.58, 3.73, 3.73, 1.82, 4.63, 3.50, 4.00, 3.67, 1.67, 4.60, 1.67, 4.00, 1.80, 4.42, 1.90, 4.63, 2.93, 3.50, 1.97,&
	4.28, 1.83, 4.13, 1.83, 4.65, 4.20, 3.93, 4.33, 1.83, 4.53, 2.03, 4.18, 4.43, 4.07, 4.13, 3.95, 4.10, 2.72, 4.58, 1.90, 4.50,&
	1.95, 4.83, 4.12/
	double precision dt
	integer i
	dt = (maxval(x) - minval(x))/tinsize
	tin = (/(((i-1)*dt+minval(x)), i=1,tinsize)/)
	thist(1:tinsize)=tin(:)
	thist(tinsize+1)=tin(tinsize)+dt
	thist = thist - dt/2
	yhist=hist(x, thist)
	!call plothist(yhist)
	N=sum(yhist)
	y_hist=real(yhist)/real(N)/dt
	call opt(y_hist, x, N, dt)
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

  subroutine opt(y_hist, x, N, dt)
	double precision, intent(in) :: y_hist(tinsize), dt, x(xsize)
	integer, intent(in) :: N
	integer, parameter :: maxiter = 20
	double precision, parameter :: tol = 10e-5
	double precision, parameter :: phi = (5**0.5 + 1) / 2
	double precision :: cost(maxiter), Win(maxiter), dummy(tinsize), yh1(tinsize), yh2(tinsize)
	double precision :: y(tinsize)
	double precision :: Winmin, Winmax, a, b, c1, c2, f1, f2, optw
	integer :: kiter
	kiter=0
	cost=0.
	Win=0.
	Winmin = 2*dt
	Winmax = maxval(x) - minval(x)
	a=ilogexp(Winmin)
	b=ilogexp(Winmax)
	c1=(phi - 1)*a + (2 - phi)*b
	c2=(2 - phi)*a + (phi - 1)*b
	call costfunction(dummy, f1, y_hist, N, logexp(c1), dt)
	!call costfunction(dummy, f2, y_hist, N, logexp(c2), dt)
	!do while ( (abs(a-b) > tol*(abs(c1)+abs(c2))) .and. (kiter < maxiter) )
	!  if (f1 < f2) then
	!     b=c2
	!	 c2=c1
	!	 c1=(phi-1)*a + (2-phi)*b
	!	 f2=f1
	!	 call costfunction(yh1, f1, y_hist, N, logexp(c1), dt)
	!	 Win(kiter)=logexp(c1)
	!	 cost(kiter)=f1
	!	 optw=logexp(c1)
	!	 y=yh1/sum(yh1*dt)
	!  else
	!	 a=c1
	!	 c1=c2
	!	 c2=(2-phi)*a + (phi-1)*b
	!	 f1=f2
	!	 call costfunction(yh2, f2, y_hist, N, logexp(c2), dt)
	!	 Win(iter)=logexp(c2)
	!	 cost(kiter)=f2
	!	 optw=logexp(c2)
	!	 y=yh2/sum(yh2*dt)
	!  endif
	!  kiter=kiter+1
	!enddo
  end subroutine opt

  subroutine costfunction(yh, cost, y_hist, N, w, dt)
	integer, intent(in) :: N
	double precision, intent(in) ::  y_hist(tinsize), w, dt
	double precision, intent(inout) :: yh(tinsize), cost
	yh=fftkernel(y_hist, w/dt)
	!cost=sum(yh**2)*dt - 2*sum(yh*y_hist)*dt + 2/(2*3.14)**0.5/w/N
	!cost=cost*N**2
  end subroutine costfunction

  function fftkernel(x, w)
	include 'fftw3.f'
	double precision, intent(in) :: x(tinsize), w
	double precision :: fftkernel(tinsize), Lmax
	integer :: L, N, i
	integer(8) :: plan
	double complex, allocatable :: input(:), output(:), xinput(:), Xoutput(:), Youtput(:)
	double precision, allocatable :: f(:), K(:)
	L=tinsize
	Lmax=L+3*w
	N = 2**ceiling(log(Lmax)/log(2.))
	print *, 'N=', N
	allocate(input(N), output(N), xinput(N), Xoutput(N), f(N), K(N), Youtput(N))
	call dfftw_plan_dft_1d(plan, N, input, output, FFTW_FORWARD, FFTW_ESTIMATE)
	xinput(1:tinsize)=x
	xinput(tinsize+1:N)=0.0
	call dfftw_execute_dft(plan, xinput, Xoutput)
	f(1:N/2) = -1.0*real((/(i, i=0, N/2, 1)/))/real(N)
	f(N/2+1:N) = real((/(i, i=1, N/2 - 1, 1)/))/real(N)
	K = exp(-0.5*(w*2*3.14*f)**2)
	call dfftw_destroy_plan(plan)
	call dfftw_plan_dft_1d(plan, N, input, output, FFTW_BACKWARD, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, Xoutput*K, Youtput)
	print *, Youtput(1:20)
	call dfftw_destroy_plan(plan)
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

end module sskernel


program main
  use sskernel
  call ssk
  stop
  end

