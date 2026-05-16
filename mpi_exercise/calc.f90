program main
implicit none
integer(8), parameter :: N = 100000000
real(8), dimension(N) :: x 
real(8),parameter :: alpha = 0.01 
integer(8) :: t1, t2, t_rate, t_max
x(:) = 1.0
call system_clock(t1)
call scal(x, alpha, N)
call system_clock(t2, t_rate, t_max)
write(*,*) (t2-t1)/dble(t_rate) 
end program main


subroutine scal(x, alpha, N) 
implicit none
integer(8), intent(in) :: N
real(8), dimension(N), intent(inout) :: x 
real(8), intent(in) :: alpha
x(:) = x(:) * alpha
end subroutine scal
