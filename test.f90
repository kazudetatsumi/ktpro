subroutine test(n, A)
  implicit none
  integer(4), intent(in) :: n
  real(8),    intent(in) :: A(n, n)

  write(*,"(a)")"** write from fortran"
  write(*,"(i8)")n
  write(*,"(1p3e12.4)")A(1,1), A(1,2), A(1,3)
  write(*,"(1p3e12.4)")A(2,1), A(2,2), A(2,3)
  write(*,"(1p3e12.4)")A(3,1), A(3,2), A(3,3)
end subroutine test
