module costfort2d
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    integer(c_int) :: len1
    type(c_ptr) :: arr
    type(c_ptr) :: kavearr
    type(c_ptr) :: darr
  end type result

contains

  function cost2d(maxw1, maxw0, Al1, Al0, A) bind(C, name="cost2d")
    !DEC$ ATTRIBUTES DLLEXPORT :: cost2d
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    real(c_double), intent(in) :: A(Al0, Al1)                         
    type(result) :: cost2d                                  
    real(c_double), pointer :: k(:,:)                    
    real(c_double), pointer :: Cn(:,:)                    
    real(c_double), pointer :: kaves(:,:)                    
    real(c_double), pointer :: deltas(:,:)                    
    integer nw0, nw1
    integer N0, N1
    integer i, ihead, j, jhead
    real kave,v
    allocate(k(Al0, Al1))
    allocate(Cn(maxw0, maxw1))
    allocate(kaves(maxw0, maxw1))
    allocate(deltas(maxw0, maxw1))
    Cn(:,:) = 0.0
    kaves(:,:) = 0.0
    deltas(:,:) = 0.0
    do nw0 = 1, maxw0
       N0 = (Al0 - mod(Al0, nw0)) / nw0
       do nw1 = 1, maxw1
          N1 = (Al1 - mod(Al1, nw1)) / nw1
          do i = 1, N0
             ihead = i*nw0 
             do j = 1, N1
                jhead = j*nw1 
                if ( i == 1 .and. j == 1) then
                    k(i, j) = A(ihead, jhead)
                else if ( i /= 1 .and. j == 1) then
                    k(i, j) = A(ihead, jhead) - A(ihead - nw0, jhead)
                else if ( i == 1 .and. j /= 1) then
                    k(i, j) = A(ihead, jhead) - A(ihead, jhead - nw1)
                else 
                    k(i, j) = A(ihead, jhead) - A(ihead - nw0, jhead) - A(ihead, jhead - nw1) + A(ihead - nw0, jhead - nw1)
                end if
             end do
          end do
          kave = sum(k(1:N0, 1:N1)) / real(N0*N1)
          v = sum((k(1:N0, 1:N1) - kave)**2) / real(N0*N1)
          !print *, "cost for ", nw0, nw1, ":", (2.0 * kave - v) / (real(nw0*nw1)**2)
          Cn(nw0, nw1) = (2.0 * kave - v) / (real(nw0*nw1)**2)
          deltas(nw0, nw1) = real(nw0*nw1)
          kaves(nw0, nw1) = kave
       end do
    end do

    print *, "minloc Cn:", minloc(Cn)
    cost2d%len0 =  maxw1
    cost2d%len1 =  maxw0
    cost2d%arr = C_loc(Cn(1:maxw0, 1: maxw1))
    cost2d%kavearr = C_loc(kaves(1:maxw0, 1:maxw1))
    cost2d%darr = C_loc(deltas(1:maxw0, 1:maxw1))
  end function cost2d

!  subroutine delete_array(arr_length, array) bind(C, name="delete_array")
!    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array
!    integer(c_int), intent(in) :: arr_length
!    type(c_ptr), value :: array
!    real(c_double), pointer :: work_array(:)
!
!    call C_F_pointer(array, work_array, [arr_length])
!    deallocate(work_array)
!    array = C_NULL_PTR
!    print *, "Is work_array freed ?, work_array: ",work_array
!  end subroutine delete_array

end module costfort2d
