program test
implicit none
integer::elem,i
real,dimension(:,:),allocatable::A
elem = 10
allocate(A(elem,2))
FORALL(I = 1:elem) A(I, 1) = i 
FORALL(I = 1:elem) A(I, 2) = 0
print*,A, "<-before"
where (A(:,1)>=3 .and. A(:,1)<=6) A(:,1) = 1 
print*,A
deallocate(A)
print*,mod(-1.1, 1.0)
print*,mod(-0.9, 1.0)
print*,mod(1.6, 1.0)
print*,mod(1.4, 1.0)
end program test 
