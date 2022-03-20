!fortran 90 test program for a kernel density estimation with localy optimized bandwidths.
!Winfun is only Boxcar (hard-coded)
!Kazuyoshi TATSUMI 2022/02/24

module ssvkernel
  use ISO_C_BINDING
  use mpi
  implicit none
  include 'fftw3.f'
  integer :: xsize, tinsize, M
  integer comm, psize, rank, ierr
  character(6) :: WinFunc
  double precision, parameter :: pi  = 4 * atan (1.0_8)
contains

  subroutine ssvk(comm0, M0, winparam, xsize0, tinsize0, WinFuncNo, xdat, tin, optw, yopt) bind(C, name="ssvk")
    integer, intent(in) :: comm0, M0, xsize0, tinsize0, WinFuncNo
    double precision, intent(in) :: winparam, xdat(xsize0), tin(tinsize0)
    double precision, intent(out) :: optw(tinsize0), yopt(tinsize0)
    double precision :: xdatstd(xsize0), xdatstddiff(xsize0-1), xdatstddiffstd(xsize0-1)
    double precision :: thist(tinsize0+1), y_hist(tinsize0), yh(tinsize0), yhist(tinsize0)
    double precision T, dt_samp, dt, cost, Wins(M0), dw, wi, Win, nsmpl
    double precision, dimension(M0, tinsize0) :: cfxw, optws, C_local
	double precision, dimension(M0*tinsize0) :: cfxw1d, C_local1d
	double precision, allocatable :: rcfxw1d(:), rC_local1d(:)
    integer :: minkbwidx(tinsize0)
    !integer nbin, tinsize2
    integer i, kbwidx, winidx, xchidx
	comm = comm0
	call MPI_Comm_size(comm, psize, ierr)
	call MPI_Comm_rank(comm, rank, ierr)
	!print *, "comm, rank, psize", comm, rank, psize
	!if (rank==0) call clock('start')
	if (WinFuncNo==1) then
		WinFunc='Boxcar'
	elseif (WinFuncNo==2) then
		WinFunc='Gauss '
	elseif (WinFuncNo==3) then
		WinFunc='Cauchy'
	endif
	!print *, WinFunc
    tinsize=tinsize0
    xsize=xsize0
    M=M0
    T=maxval(xdat)-minval(xdat)
    !xdatstd=xdat
    !call quicksort(xdatstd, 1, xsize)
    !xdatstddiff=xdatstd(2:xsize)-xdatstd(1:xsize-1)
    !xdatstddiffstd=xdatstddiff
    !call quicksort(xdatstddiffstd, 1, size(xdat)-1)
    !dt_samp=minval(pack(xdatstddiffstd, xdatstddiffstd > 0.))
    !if (ceiling(T/dt_samp) > 1e3) then
    !    tinsize2 = 1e3
    !else
    !    tinsize2 = ceiling(T/dt_samp)
    !endif
    !print *, 'tinsize2=',tinsize2
    !dt=T/(tinsize-1)
    !tin = (/(((xchidx-1)*dt+minval(xdat)), xchidx=1,tinsize)/)
	dt=minval(tin(2:)-tin(1:tinsize-1))
    thist(1:tinsize)=tin(:)
    thist(tinsize+1)=tin(tinsize)+dt
    thist = thist - dt/2
    yhist=hist(xdat, thist)
    nsmpl=sum(yhist)
    y_hist=yhist/dt
    dw=(ilogexp(T)-ilogexp(winparam*dt))/(M-1)
    ! Wins contains all widths to be considered for kernels as well as window functions, 
    ! playing a dual role to put kernel band-widths as well as window-widths.
    Wins=logexparr( (/( (i-1)*dw + ilogexp(winparam*dt), i=1,M)/) )
    !print *, "check ssvkernel param", winparam, M
    !print *, "check Ws", Wins(1:10)
	!if (rank==0) call clock('Wins ')
    !integrand of cost func, for fixed kernel band-widths
	!!mpi
    cfxw=0.
	!cfxw1d=0.
	!print *, 'CHK', M*tinsize/psize
	!print *, 'CHK', M/psize
	allocate(rcfxw1d(M*tinsize/psize))
    !do kbwidx=1, M  ! This loop can be parallelized by using mpi library.
	do kbwidx=1+rank*M/psize, (rank+1)*M/psize
      wi=Wins(kbwidx)
      yh=fftkernel(y_hist, wi/dt)
	  !print *, 'CHKyh', size(yh)
	  !print *, 'CHKrcfxw1d', size(rcfxw1d((kbwidx-1)*tinsize+1:kbwidx*tinsize))
	  !cfxw(kbwidx,:)=yh**2 - 2*yh*y_hist + 2./(2*pi)**0.5/wi*y_hist
	  rcfxw1d((kbwidx-1-rank*M/psize)*tinsize+1:(kbwidx-rank*M/psize)*tinsize)=&
&                                                  yh**2 - 2*yh*y_hist + &
&                                                  2./(2*pi)**0.5/wi*y_hist
    enddo
	!call mpi_barrier(comm, ierr)
	call mpi_allgather(rcfxw1d, M*tinsize/psize, mpi_double_precision,&
&                      cfxw1d, M*tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
    cfxw = transpose(reshape((cfxw1d), (/tinsize, M/)))
	!if (rank==0) call clock('cfxw ')
    !optws is a conversion maxtrix containing an optimum kernel band width for a pair of
    ! a window width and a x channel.
    optws=0.
    allocate(rC_local1d(M*tinsize/psize))
    do winidx=1, M     ! do loop wrt window-widths This loop can be parallelized
      Win=Wins(winidx) ! by using mpi library.
      C_local=0.
	  !mpi
      !do kbwidx=1, M   ! do loop wrt kernel band-widths
	  do kbwidx=1+rank*M/psize, (rank+1)*M/psize
         !C_local(kbwidx, :)=fftkernelWin(cfxw(kbwidx,:), Win/dt)
		 rC_local1d((kbwidx-1-rank*M/psize)*tinsize+1:(kbwidx-rank*M/psize)*tinsize)=&
&                               fftkernelWin(cfxw(kbwidx,:), Win/dt)
      enddo
	  !call mpi_barrier(comm, ierr)
	  call mpi_allgather(rC_local1d, M*tinsize/psize, mpi_double_precision,&
&                        C_local1d, M*tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
      C_local=transpose(reshape((C_local1d), (/tinsize, M/)))
      minkbwidx=minloc(C_local, 1)  
      do xchidx=1, tinsize ! do loop wrt x channels
         optws(winidx, xchidx) = Wins(minkbwidx(xchidx))
      enddo
    enddo
	!if (rank==0) call clock('optss')
    call opt(optw, yopt, y_hist, xdat, nsmpl, tin, dt, Wins, optws)
	!if (rank==0) call clock('optff')
  end subroutine ssvk

  function hist(x, th)
    double precision, intent(in) :: x(xsize), th(tinsize+1)
	double precision :: hist(tinsize)
    integer ix, it
    hist(:) = 0.
    do ix = 1, xsize
       do it = 1, tinsize
          if ( x(ix) >= th(it) .and. x(ix) < th(it+1) ) then
               hist(it) = hist(it) + 1.
          end if
       end do
    end do
  end function hist
    
  !subroutine plothist(yhist)
  !  integer, intent(in) ::  yhist(tinsize)
  !  integer ix, j
  !  do ix = 1, tinsize
  !    print *, ("*", j=1,yhist(ix))
  !  end do
  !end subroutine plothist

  subroutine opt(optw, yopt, y_hist, xdat, nsmpl, tin, dt, Wins, optws)
    double precision, intent(in) :: y_hist(tinsize), dt, xdat(xsize), tin(tinsize)
    double precision, intent(in) :: Wins(M), optws(M, tinsize), nsmpl
    !integer, intent(in) :: nsmpl
    integer, parameter :: maxiter = 30
    double precision, parameter :: tol = 10e-5
    double precision, parameter :: phi = (5**0.5 + 1) / 2
    double precision :: cost(maxiter), gs(maxiter)
    double precision, dimension(tinsize) :: dummy, dummy2, yh1, yh2, optwp1, optwp2, yv1, yv2
    double precision, dimension(tinsize), intent(out) :: yopt, optw
    double precision :: a, b, c1, c2, f1, f2
    integer :: kiter
    kiter=1
    gs=0.
    cost=0.
    a=1e-12
    b=1
    c1=(phi - 1)*a + (2 - phi)*b
    c2=(2 - phi)*a + (phi - 1)*b
    call costfunction(f1, dummy, dummy2, y_hist, nsmpl, tin, dt, optws, Wins, c1)
	!print *, "CHK", f1, c1, rank
    call costfunction(f2, dummy, dummy2, y_hist, nsmpl, tin, dt, optws, Wins, c2)
	!print *, "CHK", f2, c2, rank
    do while ( (abs(a-b) > tol*(abs(c1)+abs(c2))) .and. (kiter <= maxiter) )
      if (f1 < f2) then
         b=c2
         c2=c1
         c1=(phi-1)*a + (2-phi)*b
         f2=f1
         call costfunction(f1, yv1, optwp1, y_hist, nsmpl, tin, dt, optws, Wins, c1)
         yopt=yv1/sum(yv1*dt)
         optw=optwp1
      else
         a=c1
         c1=c2
         c2=(2-phi)*a + (phi-1)*b
         f1=f2
         call costfunction(f2, yv2, optwp2, y_hist, nsmpl, tin, dt, optws, Wins, c2)
         yopt=yv2/sum(yv2*dt)
         optw=optwp2
      endif
      gs(kiter)=c1
      cost(kiter)=f1
	  print *, rank, kiter, cost(kiter), gs(kiter)
      kiter=kiter+1
    enddo
  end subroutine opt

  subroutine costfunction(Cg, yv, optwp, y_hist, nsmpl, tin, dt, optws, Wins, g)
    !integer, intent(in) :: nsmpl
    double precision, dimension(tinsize), intent(in) ::  y_hist, tin
    double precision, intent(in) :: optws(M, tinsize), Wins(M), dt, g, nsmpl
    double precision, intent(out) :: Cg, yv(tinsize), optwp(tinsize)
    double precision, dimension(tinsize) :: optwv, cintegrand, Z
    double precision, allocatable :: y_hist_nz(:), tin_nz(:), roptwp(:), ryv(:)
    double precision :: gammas(M)
    integer :: xchidx, maxidx, wchidx
	!integer comm, psize, rank, ierr
	!if (rank==0) call clock('cost0')
    optwv=0.
    do xchidx=1, tinsize  
      gammas = optws(:, xchidx)/Wins
      if (g > maxval(gammas)) then
        optwv(xchidx)=minval(Wins)
      else
        if (g < minval(gammas)) then
          optwv(xchidx)=maxval(Wins)
        else
          maxidx=maxval(pack([(wchidx, wchidx=1, M)], gammas >= g))
          optwv(xchidx)=g*Wins(maxidx)
        endif
      endif
    enddo
    optwp=0.
	! Nadaraya-Watson kernel regression to smooth optw.
	!if (rank==0) call clock('cost1')
	!!mpi
    !!do xchidx=1, tinsize
	allocate(roptwp(tinsize/psize), ryv(tinsize/psize))
	do xchidx=1+rank*tinsize/psize, (rank+1)*tinsize/psize
	  if (WinFunc == 'Boxcar') Z=Boxcar(tin(xchidx)-tin, optwv/g)
	  if (WinFunc == 'Gauss') Z=vGauss(tin(xchidx)-tin, optwv/g)
	  if (WinFunc == 'Cauchy') Z=Cauchy(tin(xchidx)-tin, optwv/g)
      !!optwp(xchidx)=sum(optwv*Z)/sum(Z)
      roptwp(xchidx-rank*tinsize/psize)=sum(optwv*Z)/sum(Z)
    enddo
	!call mpi_barrier(comm, ierr)
	call mpi_allgather(roptwp, tinsize/psize, mpi_double_precision,&
&                     optwp, tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
	!if (rank==0) call clock('cost2')
	! Balloon estimator only on non-zero bins.
    y_hist_nz=pack(y_hist, y_hist > 0.) 
    tin_nz=pack(tin, y_hist>0)
    yv = 0.
	!! mpi
    !!do xchidx=1, tinsize
	do xchidx=1+rank*tinsize/psize, (rank+1)*tinsize/psize
      ryv(xchidx-rank*tinsize/psize)=sum(y_hist_nz*dt*Gauss(tin(xchidx)-tin_nz, optwp(xchidx)))
    enddo
	!call mpi_barrier(comm, ierr)
	call mpi_allgather(ryv, tinsize/psize, mpi_double_precision,&
&                      yv, tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
    yv=yv*nsmpl/sum(yv*dt)
    cintegrand = yv**2 - 2.*yv*y_hist + 2./(2.*pi)**0.5/optwp*y_hist
    Cg=sum(cintegrand*dt)
	!if (rank==0) call clock('cost3')
  end subroutine costfunction

  function fftkernel(x, w)
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

  function fftkernelWin(x, w)
    double precision, intent(in) :: x(tinsize), w
    double precision :: fftkernelWin(tinsize), Lmax, a
    integer :: L, N, i
    integer(8) :: plan
    double complex, allocatable :: input(:), output(:), xinput(:), Xoutput(:), Youtput(:), K(:), t(:)
    double precision, allocatable :: f(:), expa(:)
    L=tinsize
    Lmax=L+3*w
    N = 2**ceiling(log(Lmax)/log(2.))
    allocate(input(N), output(N), xinput(N), Xoutput(N), f(N), expa(N), K(N), Youtput(N), t(N))
    call dfftw_plan_dft_1d(plan, N, input, output, FFTW_FORWARD, FFTW_ESTIMATE)
    xinput(1:tinsize)=x
    xinput(tinsize+1:N)=0.0
    call dfftw_execute_dft(plan, xinput, Xoutput)
    f(1:N/2+1) = -1.0*real((/(i, i=0, N/2, 1)/))/real(N)
    f(N/2+2:N) = 0.5 - real((/(i, i=1, N/2 - 1, 1)/))/real(N)
    ! for preventing the procedure from underflow exception.
    !expa=-0.5*(w*2*pi*f)**2
    !K = 0.
    !where (expa > -708) K=exp(expa)
	if (WinFunc == 'Boxcar') then
	   t=2*pi*f
	   a=12**0.5*w
	   K(2:)=2*sin(a*t(2:)/2)/(a*t(2:))
	   K(1)=1.
	elseif (WinFunc == 'Gauss') then
	   K = exp(-0.5*(w*t)**2)
	elseif (WinFunc == 'Cauchy') then
	   K = exp(-w*abs(t))
	endif
    call dfftw_destroy_plan(plan)
    call dfftw_plan_dft_1d(plan, N, input, output, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, Xoutput*K, Youtput)
    call dfftw_destroy_plan(plan)
    ! Computing a fftw forward followed by a fftw backward results in the original
    ! array scaled by N. So, we should divid the result by N so as to obtain the
    ! count dependent cost fucntion.
    fftkernelWin=real(Youtput(1:L)/N)
  end function fftkernelWin

  function Gauss(x, w)
    double precision, intent(in) :: x(:), w
    double precision :: Gauss(size(x))
    Gauss =  1. / (2 * pi)**2 / w * exp(-x**2 / 2. / w**2)
  end function Gauss

  function ilogexp(x)
    double precision, intent(in) :: x
    double precision :: ilogexp
    if ( x < 1e2 ) then
       ilogexp = log(exp(x) - 1)
    else
       ilogexp = x
    endif 
  end function ilogexp

  function ilogexparr(x)
    double precision, intent(in) :: x(:)
    double precision :: ilogexparr(size(x))
    where(x<1e2) ilogexparr(:)=log(exp(x(:)) - 1)
    where(x>=1e2) ilogexparr(:)=x(:)
  end function ilogexparr

  function logexp(x)
    double precision, intent(in) :: x
    double precision :: logexp
    if ( x < 1e2 ) then
       logexp = log(exp(x) + 1)
    else
       logexp = x
    endif 
  end function logexp

  function logexparr(x)
    double precision, intent(in) :: x(:)
    double precision :: logexparr(size(x))
    where(x<1e2) logexparr=log(exp(x) + 1)
    where(x>=1e2) logexparr=x
  end function logexparr

  function Boxcar(x,w)
    double precision, intent(in) :: x(tinsize), w(tinsize)
	double precision :: Boxcar(tinsize)
    double precision :: a(tinsize)
	a=12**0.5*w
	Boxcar=1/a
    where(abs(x) > a/2) Boxcar=0.
  end function Boxcar

  function vGauss(x, w)
    double precision, intent(in) :: x(tinsize), w(tinsize)
    double precision :: vGauss(size(x))
    vGauss =  1. / (2 * pi)**2 / w * exp(-x**2 / 2. / w**2)
  end function vGauss

  function Cauchy(x,w)
    double precision, intent(in) :: x(tinsize), w(tinsize)
	double precision :: Cauchy(tinsize)
	Cauchy= 1. / (pi * w * (1. + (x/w)**2))
  end function Cauchy

! quicksort.f -*-f90-*-
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

  subroutine clock(flag)
    implicit none
    !character(8)  :: date
	character(5), intent(in) :: flag
    character(10) :: time
    !character(5)  :: zone
    !integer,dimension(8) :: values
    !call date_and_time(date,time,zone,values)
    !call date_and_time(DATE=date,ZONE=zone)
    call date_and_time(TIME=time)
    !call date_and_time(VALUES=values)
    !print '(a,2x,a,2x,a)', time
	print *, time, flag
  end subroutine clock

end module ssvkernel


!program main
!  use sskernel
!  call ssk
!  stop
!  end
