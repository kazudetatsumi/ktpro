module costfort
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    type(c_ptr) :: arr
    type(c_ptr) :: kavearr
    type(c_ptr) :: darr
  end type result

contains

  function cost1d(maxw, Al0, A) bind(C, name="cost1d")
    !DEC$ ATTRIBUTES DLLEXPORT :: cost1d
    integer(c_int), intent(in) :: maxw
    integer(c_int), intent(in) :: Al0
    real(c_double), intent(in) :: A(Al0)                         
    type(result) :: cost1d                                  
    real(c_double), pointer :: k(:)                    
    real(c_double), pointer :: Cn(:)                    
    real(c_double), pointer :: kaves(:)                    
    real(c_double), pointer :: deltas(:)                    
    integer nw0
    integer N0
    integer i, ihead
    real kave,v
    allocate(k(Al0))
    allocate(Cn(maxw))
    allocate(kaves(maxw))
    allocate(deltas(maxw))
    k(:) = 0.0
    Cn(:) = 0.0
    kaves(:) = 0.0
    deltas(:) = 0.0
    do nw0 = 1, maxw
       N0 = (Al0 - mod(Al0, nw0)) / nw0
       do i = 1, N0
          ihead = i*nw0 
          if ( i == 1) then
              k(i) = A(ihead)
          else
              k(i) = A(ihead) - A(ihead - nw0)
          end if
       end do
       kave = sum(k(1:N0)) / real(N0)
       v = sum((k(1:N0) - kave)**2) / real(N0)
       Cn(nw0) = (2.0 * kave - v) / (real(nw0)**2)
       deltas(nw0) = real(nw0)
       kaves(nw0) = kave
    end do

    cost1d%len0 =  nw0
    cost1d%arr = C_loc(Cn(1:nw0))
    cost1d%kavearr = C_loc(kaves(1:nw0))
    cost1d%darr = C_loc(deltas(1:nw0))
  end function cost1d

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

end module costfort
