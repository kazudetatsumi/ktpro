program fftw_c2c
  implicit none
  include 'fftw3.f'

  real(8), parameter :: PI = atan(1.0d0) * 4
  integer, parameter :: N = 16
  integer :: j, k

  integer(8) :: plan
  !complex(8) :: in(0:N-1), out(0:N-1)
  complex(8) :: in(N), out(N)

  ! 1. create a FFT plan first
  call dfftw_plan_dft_1d(plan, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE)

  !do j = 0, N-1
  do j = 1, N
    in(j) = sin(2 * PI * (j-1) / N)
  end do

  ! 2. execute FFT (as many times as you want)
  call dfftw_execute_dft(plan, in, out)

  !do k = 0, N/2
  do k = 1, N/2
    write(*, *) k, k, out(k)
  end do
  !do k = N/2+1, N-1
  do k = N/2+1, N
    write(*, *) k, k - N, out(k)
  end do

  ! 3. destroy the plan finally
  call dfftw_destroy_plan(plan)

  stop
end program fftw_c2c
