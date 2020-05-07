subroutine scal(x, alpha, N) bind(C) 
  implicit none
  integer(8), intent(in) :: N
  real(8), dimension(N), intent(inout) :: x 
  real(8), intent(in) :: alpha
  x(:) = x(:) * alpha 
end subroutine scal
