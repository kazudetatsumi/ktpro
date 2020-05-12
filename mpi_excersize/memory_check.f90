subroutine memory_check(datasize3, datasize2, datasize1, datasize0, A, condition) bind(C)
implicit none
integer, intent(in) :: datasize0, datasize1, datasize2, datasize3
double precision, intent(inout) :: A(datasize0, datasize1, datasize2, datasize3)
integer, intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)
integer i,j,k,l

do i=1, datasize0
do j=1, datasize1
do k=1, datasize2
do l=1, datasize3
   A(i,j,k,l) = i*j*k*l*1.00
enddo
enddo
enddo
enddo
print *, sum(A)

print *, 'Check the memory usage now'
read *
end subroutine memory_check
