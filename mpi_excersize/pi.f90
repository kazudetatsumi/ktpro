program pi
implicit none
integer, parameter :: SP = kind(1.0)
integer, parameter :: DP = selected_real_kind(2*precision(1.0_SP))
integer, parameter :: n = 10000000
integer :: i
real(DP) :: x, dx, p

print *, DP
dx = 1.0_DP/real(n,DP)
p = 0.0_DP
do i = 1, n
x = real(i, DP) * dx
  p = p + 4.0_DP/(1.0_DP + x ** 2)*dx
end do
print *, p
end program pi
