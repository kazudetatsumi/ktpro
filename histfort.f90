module histfort
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    type(c_ptr) :: arr
  end type result

contains

  function hist1d(A, Al0, nw0) bind(C, name="hist1d")
    !DEC$ ATTRIBUTES DLLEXPORT :: plus
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: nw0
    real(c_double), intent(in) :: A(Al0)                         
    type(result) :: hist1d                                  
    real(c_double), pointer :: k(:)                    
    integer N0
    integer i, ihead
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    allocate(k(N0))

    do i = 1, N0
       ihead = (i+1)*nw0 - 1
       if ( i == 0) then
           k(i) = A(ihead)
       else
           k(i) = A(ihead) - A(ihead - nw0)
       end if
    end do

    hist1d%len0 =  N0
    hist1d%arr = C_loc(k)
  end function hist1d




  subroutine delete_array(Bl0, darray) bind(C, name="delete_array")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array
    integer(c_int), intent(in) :: Bl0
    type(c_ptr), value :: darray
    real(c_double), pointer :: k(:)

    call C_F_pointer(darray, k, [Bl0])
    print *, "k:",k
    deallocate(k)
    darray = C_NULL_PTR
    print *, "Is work_array freed ?, work_array: ",k
  end subroutine delete_array

end module histfort
