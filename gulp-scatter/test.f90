module TestMod
	implicit none
contains
	subroutine sub(a,m)   
 	  integer ::  m
 	  real(8) :: a(4,*)
 	    print *, "a =",  a(:,1)
 	end subroutine
end module

program main 
  use TestMod
  implicit none
  integer, parameter :: n = 2
  integer :: i, j
  real(8) :: a(n,n)
	do i = 1, n
	  do j = 1, n
	    a(i,j) = i*j
	  end do
	end do
	call sub(a, n*n)
end program
	
