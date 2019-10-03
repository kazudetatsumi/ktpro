module histfort
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    type(c_ptr) :: arr
  end type result

contains

  function hist1d(nw0, Al0, A) bind(C, name="hist1d")
    !DEC$ ATTRIBUTES DLLEXPORT :: plus
    integer(c_int), intent(in) :: nw0
    integer(c_int), intent(in) :: Al0
    real(c_double), intent(in) :: A(Al0)                         
    type(result) :: hist1d                                  
    real(c_double), pointer :: work_array(:)                    
    integer N0
    integer i, ihead
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    allocate(work_array(N0))

    do i = 1, N0
       ihead = i*nw0 
       if ( i == 1) then
           work_array(i) = A(ihead)
       else
           work_array(i) = A(ihead) - A(ihead - nw0)
       end if
    end do

    hist1d%len0 =  N0
    hist1d%arr = C_loc(work_array)
  end function hist1d




  subroutine delete_array(arr_length, array) bind(C, name="delete_array")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array
    integer(c_int), intent(in) :: arr_length
    type(c_ptr), value :: array
    real(c_double), pointer :: work_array(:)

    call C_F_pointer(array, work_array, [arr_length])
    deallocate(work_array)
    array = C_NULL_PTR
    print *, "Is work_array freed ?, work_array: ",work_array
  end subroutine delete_array

end module histfort
