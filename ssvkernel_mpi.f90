!fortran 89 program for a kernel density estimation with localy optimized bandwidths.
!openmpi version: M and tinsize must be multiples of psize. 
!Kazuyoshi TATSUMI 2022/04/08

module ssvkernel
  use ISO_C_BINDING
  use mpi
  implicit none
  include 'fftw3.f'
  integer :: xsize, tinsize, M
  integer comm, psize, rank, ierr
  character(6) :: WinFunc
  integer  :: nb 
  double precision, parameter :: pi  = 4 * atan (1.0_8)
  double precision :: dt, nsmpl
contains

  subroutine ssvk(comm0, M0, winparam, xsize0, tinsize0, WinFuncNo, nb0, xdat, tin, optw, yopt, yb) bind(C, name="ssvk")
    integer, intent(in) :: comm0, M0, xsize0, tinsize0, WinFuncNo, nb0
    double precision, intent(in) :: winparam, xdat(xsize0), tin(tinsize0)
    double precision, intent(out) :: optw(tinsize0), yopt(tinsize0)
    double precision, intent(out) :: yb(tinsize0,nb0)
    double precision Wins(M0)
    double precision, dimension(M0, tinsize0) :: optws!, C_local
    double precision, allocatable :: cfxw(:,:), y_hist(:)
    integer tidx
    comm = comm0
    call MPI_Comm_size(comm, psize, ierr)
    call MPI_Comm_rank(comm, rank, ierr)
    if (WinFuncNo==1) then
        WinFunc='Boxcar'
    elseif (WinFuncNo==2) then
        WinFunc='Gauss '
    elseif (WinFuncNo==3) then
        WinFunc='Cauchy'
    endif
    tinsize=tinsize0
    xsize=xsize0
    M=M0
    nb=nb0
    if (rank==0) then
        print *, 'y_hist'
    endif
    y_hist=y_histf(xdat, tin)
    if (rank==0) then
        print *, 'Wins'
    endif
    Wins=Winsf(winparam, xdat)
    if (rank==0) then
        print *, 'cfxw'
    endif
    cfxw=cfxwf(Wins, y_hist)
    if (rank==0) then
        print *, 'optws'
    endif
    optws=optwsf(Wins, cfxw)
    if (rank==0) then
        print *, 'opt'
    endif
    call opt(optw, yopt, y_hist, xdat, tin, Wins, optws)
    yb = ybf(tin, xdat, optw) 
    !if (rank==0) then
    !do tidx=1,tinsize
    !   print *, yb(tidx, 1)
    !enddo
    !endif
  end subroutine ssvk

  function Winsf(winparam, xdat)
    double precision, intent(in) :: winparam, xdat(xsize)
    double precision T, Winsf(M), dw
    integer i
    T=maxval(xdat)-minval(xdat)
    dw=(ilogexp(T)-ilogexp(winparam*dt))/(M-1)
    ! Wins contains all widths to be considered for kernels as well as window functions, 
    ! playing a dual role to put kernel band-widths as well as window-widths.
    Winsf=logexparr( (/( (i-1)*dw + ilogexp(winparam*dt), i=1,M)/) )
  end function Winsf

  function y_histf(xdat, tin)
    double precision, intent(in) :: xdat(xsize), tin(tinsize)
    double precision :: thist(tinsize+1)
    double precision, allocatable ::  y_histf(:)
    dt=minval(tin(2:)-tin(1:tinsize-1))
    thist(1:tinsize)=tin(:)
    thist(tinsize+1)=tin(tinsize)+dt
    thist = thist - dt/2
    y_histf=hist(xdat, thist)
    nsmpl=sum(y_histf)
    y_histf=y_histf/dt
  end function y_histf

  function cfxwf(Wins, y_hist)
    double precision, intent(in) :: Wins(M)
    double precision, allocatable, intent(in) :: y_hist(:)
    double precision, allocatable :: cfxwf(:,:)
    double precision :: wi
    double precision, allocatable :: rcfxw(:,:)
    double precision, allocatable :: yh(:)
    integer kbwidx
    cfxwf=0.
    allocate(cfxwf(tinsize, M))
    allocate(rcfxw(tinsize, M/psize))
    allocate(yh(tinsize))
    do kbwidx=1+rank*M/psize, (rank+1)*M/psize
      wi=Wins(kbwidx)
      yh=fftkernel(y_hist, wi/dt)
      rcfxw(:,kbwidx-rank*M/psize)=yh**2 - 2*yh*y_hist + 2./(2*pi)**0.5/wi*y_hist
    enddo
    !call mpi_barrier(comm, ierr)
    call mpi_allgather(rcfxw, M*tinsize/psize, mpi_double_precision,&
                       cfxwf, M*tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
  end function cfxwf

!  function cfxwf(Wins, y_hist)
!    double precision, intent(in) :: Wins(M)
!    double precision, allocatable, intent(in) :: y_hist(:)
!    double precision :: cfxwf(M,tinsize), wi
!    double precision, allocatable :: cfxw1d(:)
!    double precision, allocatable :: rcfxw1d(:)
!    double precision, allocatable :: yh(:)
!    integer kbwidx
!    cfxwf=0.
!    !do kbwidx=1, M  ! This loop can be parallelized by using mpi library.
!    allocate(rcfxw1d(M*tinsize/psize))
!    allocate(yh(tinsize))
!    do kbwidx=1+rank*M/psize, (rank+1)*M/psize
!      wi=Wins(kbwidx)
!      yh=fftkernel(y_hist, wi/dt)
!      rcfxw1d((kbwidx-1-rank*M/psize)*tinsize+1:(kbwidx-rank*M/psize)*tinsize)=&
!&                                                  yh**2 - 2*yh*y_hist + &
!&                                                  2./(2*pi)**0.5/wi*y_hist
!    enddo
!    !call mpi_barrier(comm, ierr)
!    allocate(cfxw1d(M*tinsize))
!    call mpi_allgather(rcfxw1d, M*tinsize/psize, mpi_double_precision,&
!&                      cfxw1d, M*tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
!    cfxwf = transpose(reshape((cfxw1d), (/tinsize, M/)))
!  end function cfxwf

  function optwsf(Wins, cfxw)
    double precision, intent(in) :: Wins(M)
    double precision, allocatable, intent(in) :: cfxw(:,:)
    double precision, dimension(M, tinsize) :: optwsf
    double precision ::  rC_local(tinsize,M/psize), Win
    double precision :: C_local(tinsize,M)
    integer :: minkbwidx(tinsize)
    integer kbwidx, winidx, xchidx
    optwsf=0.
    do winidx=1, M     ! do loop wrt window-widths This loop can be parallelized
      Win=Wins(winidx) ! by using mpi library.
      !mpi
      !do kbwidx=1, M   ! do loop wrt kernel band-widths
      do kbwidx=1+rank*M/psize, (rank+1)*M/psize
         !C_local(kbwidx, :)=fftkernelWin(cfxw(kbwidx,:), Win/dt)
         rC_local(:,kbwidx-rank*M/psize)=fftkernelWin(cfxw(:,kbwidx), Win/dt)
      enddo
      !call mpi_barrier(comm, ierr)
      call mpi_allgather(rC_local, M*tinsize/psize, mpi_double_precision,&
&                        C_local, M*tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
      minkbwidx=minloc(C_local, 2)  
      do xchidx=1, tinsize ! do loop wrt x channels
         optwsf(winidx, xchidx) = Wins(minkbwidx(xchidx))
      enddo
    enddo
  end function optwsf

  function optwsf_alloc(Wins, cfxw)
    double precision, intent(in) :: Wins(M)
    double precision, allocatable, intent(in) :: cfxw(:,:)
    !double precision, dimension(M, tinsize) :: optwsf, C_local
    double precision ::  optwsf_alloc(M, tinsize) 
    double precision, allocatable :: C_local(:,:)
    !double precision ::  rC_local1d(M*tinsize/psize), Win
    double precision, allocatable  ::  rC_local1d(:)
    !double precision :: C_local1d(M*tinsize)
    double precision, allocatable :: C_local1d(:)
    double precision :: Win
    integer :: minkbwidx(tinsize)
    integer kbwidx, winidx, xchidx
    optwsf_alloc=0.
    do winidx=1, M     ! do loop wrt window-widths This loop can be parallelized
      Win=Wins(winidx) ! by using mpi library.
      !mpi
      !do kbwidx=1, M   ! do loop wrt kernel band-widths
      allocate(rC_local1d(M*tinsize/psize))
      do kbwidx=1+rank*M/psize, (rank+1)*M/psize
         !C_local(kbwidx, :)=fftkernelWin(cfxw(kbwidx,:), Win/dt)
         rC_local1d((kbwidx-1-rank*M/psize)*tinsize+1:(kbwidx-rank*M/psize)*tinsize)=&
&                               fftkernelWin(cfxw(:,kbwidx), Win/dt)
      enddo
      !call mpi_barrier(comm, ierr)
      allocate(C_local1d(M*tinsize))
      call mpi_allgather(rC_local1d, M*tinsize/psize, mpi_double_precision,&
&                        C_local1d, M*tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
      deallocate(rC_local1d)
      allocate(C_local(M,tinsize))
      C_local=transpose(reshape((C_local1d), (/tinsize, M/)))
      deallocate(C_local1d)
      minkbwidx=minloc(C_local, 1)  
      deallocate(C_local)
      do xchidx=1, tinsize ! do loop wrt x channels
         optwsf_alloc(winidx, xchidx) = Wins(minkbwidx(xchidx))
      enddo
    enddo
  end function optwsf_alloc

  function hist(x, th)
    double precision, intent(in) :: x(xsize), th(tinsize+1)
    double precision, allocatable :: hist(:)
    integer ix, it
    allocate(hist(tinsize))
    hist(:) = 0.
    do ix = 1, xsize
       do it = 1, tinsize
          if ( x(ix) >= th(it) .and. x(ix) < th(it+1) ) then
               hist(it) = hist(it) + 1.
          end if
       end do
    end do
  end function hist
    
  subroutine opt(optw, yopt, y_hist, xdat, tin, Wins, optws)
    double precision, intent(in) :: xdat(xsize), tin(tinsize)
    double precision, intent(in) :: Wins(M), optws(M, tinsize)
    double precision, allocatable, intent(in) :: y_hist(:)
    integer, parameter :: maxiter = 30
    double precision, parameter :: tol = 10e-5
    double precision, parameter :: phi = (5**0.5 + 1) / 2
    double precision :: cost(maxiter), gs(maxiter)
    double precision, dimension(tinsize) :: dummy, dummy2, optwp1, optwp2, yv1, yv2
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
    call costfunction(f1, dummy, dummy2, y_hist, tin, optws, Wins, c1)
    call costfunction(f2, dummy, dummy2, y_hist, tin, optws, Wins, c2)
    do while ( (abs(a-b) > tol*(abs(c1)+abs(c2))) .and. (kiter <= maxiter) )
      if (f1 < f2) then
         b=c2
         c2=c1
         c1=(phi-1)*a + (2-phi)*b
         f2=f1
         call costfunction(f1, yv1, optwp1, y_hist, tin, optws, Wins, c1)
         yopt=yv1/sum(yv1*dt)
         optw=optwp1
      else
         a=c1
         c1=c2
         c2=(2-phi)*a + (phi-1)*b
         f1=f2
         call costfunction(f2, yv2, optwp2, y_hist, tin, optws, Wins, c2)
         yopt=yv2/sum(yv2*dt)
         optw=optwp2
      endif
      gs(kiter)=c1
      cost(kiter)=f1
      print *, rank, kiter, cost(kiter), gs(kiter)
      kiter=kiter+1
    enddo
  end subroutine opt

  subroutine costfunction(Cg, yv, optwp, y_hist, tin, optws, Wins, g)
    double precision, dimension(tinsize), intent(in) :: tin
    double precision, allocatable, intent(in) :: y_hist(:)
    double precision, intent(in) :: optws(M, tinsize), Wins(M), g
    double precision, intent(out) :: Cg, yv(tinsize), optwp(tinsize)
    double precision, dimension(tinsize) :: optwv, cintegrand, Z
    double precision, allocatable :: y_hist_nz(:), tin_nz(:), roptwp(:), ryv(:)
    double precision :: gammas(M)
    integer :: xchidx, maxidx, wchidx
    y_hist_nz=pack(y_hist, y_hist > 0.) 
    tin_nz=pack(tin, y_hist>0)
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
    ! mpi for Nadaraya-Watson kernel regression to smooth optw
    !         and balloon estimation.
    !do xchidx=1, tinsize
    allocate(roptwp(tinsize/psize), ryv(tinsize/psize))
    do xchidx=1+rank*tinsize/psize, (rank+1)*tinsize/psize
      if (WinFunc == 'Boxcar') Z=Boxcar(tin(xchidx)-tin, optwv/g)
      if (WinFunc == 'Gauss') Z=vGauss(tin(xchidx)-tin, optwv/g)
      if (WinFunc == 'Cauchy') Z=Cauchy(tin(xchidx)-tin, optwv/g)
      roptwp(xchidx-rank*tinsize/psize)=sum(optwv*Z)/sum(Z)
      ryv(xchidx-rank*tinsize/psize)=sum(y_hist_nz*dt*&
                                         Gauss(tin(xchidx)-tin_nz,&
                                         roptwp(xchidx-rank*tinsize/psize)&
                                              )&
                                         )
    enddo
    call mpi_allgather(ryv, tinsize/psize, mpi_double_precision,&
                       yv, tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
    call mpi_allgather(roptwp, tinsize/psize, mpi_double_precision,&
                       optwp, tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
    deallocate(ryv,roptwp)
    yv=yv*nsmpl/sum(yv*dt)
    cintegrand = yv**2 - 2.*yv*y_hist + 2./(2.*pi)**0.5/optwp*y_hist
    Cg=sum(cintegrand*dt)
  end subroutine costfunction
  
  subroutine costfunctionorg(Cg, yv, optwp, y_hist, tin, optws, Wins, g)
    !integer, intent(in) :: nsmpl
    double precision, dimension(tinsize), intent(in) :: tin
    double precision, allocatable, intent(in) :: y_hist(:)
    double precision, intent(in) :: optws(M, tinsize), Wins(M), g
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
    deallocate(roptwp)
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
    deallocate(ryv, tin_nz, y_hist_nz)
    yv=yv*nsmpl/sum(yv*dt)
    cintegrand = yv**2 - 2.*yv*y_hist + 2./(2.*pi)**0.5/optwp*y_hist
    Cg=sum(cintegrand*dt)
    !if (rank==0) call clock('cost3')
  end subroutine costfunctionorg

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

  function ybf(tin, xdat, optw) 
    double precision, intent(in) :: tin(tinsize), xdat(xsize), optw(tinsize)
    double precision :: u(xsize), xb(xsize), thist(tinsize+1)
    double precision, allocatable :: yhistb(:)
    double precision :: ryb(tinsize, nb/psize)
    double precision :: ybf(tinsize, nb)
    double precision, allocatable :: y_histb_nz(:), tinb_nz(:)
    integer :: idx(xsize), sidx, tidx, xchidx, seed(33)
    thist(1:tinsize)=tin(:)
    thist(tinsize+1)=tin(tinsize)+dt
    thist = thist - dt/2
    seed=1000*(rank+1)
    call random_seed(put=seed)
    do sidx=1+rank*nb/psize,(rank+1)*nb/psize
      print *, sidx
      call random_number(u)
      idx=1+floor(u*xsize)
      xb=xdat(idx)
      yhistb=hist(xb, thist)
      y_histb_nz=pack(yhistb, yhistb > 0.) 
      tinb_nz=pack(tin, yhistb > 0.)
      do xchidx=1, tinsize
        ryb(xchidx,sidx-rank*nb/psize)=sum(y_histb_nz*dt*Gauss(tin(xchidx)-tinb_nz, optw(xchidx)))
      enddo
      ryb(:,sidx-rank*nb/psize)=ryb(:,sidx-rank*nb/psize)/sum(ryb(:,sidx-rank*nb/psize)*dt)
    enddo
    call mpi_allgather(ryb, tinsize*nb/psize, mpi_double_precision,&
                       ybf, tinsize*nb/psize, mpi_double_precision, mpi_comm_world, ierr)
  end function ybf

  function ybforg(tin, xdat, optw) 
    double precision, intent(in) :: tin(tinsize), xdat(xsize), optw(tinsize)
    double precision :: u(xsize), xb(xsize), thist(tinsize+1), yvb(tinsize)
    double precision, allocatable :: yhistb(:)
    double precision :: ryvb(tinsize/psize)
    double precision :: ybforg(nb, tinsize)
    double precision, allocatable :: y_histb_nz(:), tinb_nz(:)
    integer :: idx(xsize), sidx, tidx, xchidx
    thist(1:tinsize)=tin(:)
    thist(tinsize+1)=tin(tinsize)+dt
    thist = thist - dt/2
    do sidx=1,nb
      if (rank==0) then
        call random_number(u)
        idx=1+floor(u*xsize)
        xb=xdat(idx)
      endif
      call mpi_bcast(xb(1), xsize, mpi_double_precision, 0, mpi_comm_world, ierr)
      yhistb=hist(xb, thist)
      y_histb_nz=pack(yhistb, yhistb > 0.) 
      tinb_nz=pack(tin, yhistb > 0.)
      yvb = 0.
      do xchidx=1+rank*tinsize/psize, (rank+1)*tinsize/psize
        ryvb(xchidx-rank*tinsize/psize)=sum(y_histb_nz*dt*Gauss(tin(xchidx)-tinb_nz, optw(xchidx)))
      enddo
      call mpi_allgather(ryvb, tinsize/psize, mpi_double_precision,&
                         yvb, tinsize/psize, mpi_double_precision, mpi_comm_world, ierr)
      yvb=yvb/sum(yvb*dt)
      ybforg(sidx,:)=yvb
    enddo

    if (rank==0) then
    do tidx=1,tinsize
       print *, ybforg(1, tidx)
    enddo
    endif
  end function ybforg


end module ssvkernel


!program main
!  use sskernel
!  call ssk
!  stop
!  end

