module histfort
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    integer(c_int) :: len1
    type(c_ptr) :: arr
  end type result

contains

  function hist1d(nw0, Al0, A) bind(C, name="hist1d")
    !DEC$ ATTRIBUTES DLLEXPORT :: hist1d
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



  function hist2d(nw1, nw0, Al1, Al0, A) bind(C, name="hist2d")
    !DEC$ ATTRIBUTES DLLEXPORT :: hist2d
    integer(c_int), intent(in) :: nw0
    integer(c_int), intent(in) :: nw1
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    real(c_double), intent(in) :: A(Al0, Al1)                         
    type(result) :: hist2d                                  
    real(c_double), pointer :: work_array(:, :)                    
    integer N0, N1
    integer i, ihead, j, jhead
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    N1 = (Al1 - mod(Al1, nw1)) / nw1
    allocate(work_array(N0, N1))

    do i = 1, N0
       ihead = i*nw0 
       do j = 1, N1
          jhead = j*nw1
          if ( i == 1 .and. j == 1 ) then
             work_array(i, j) = A(ihead, jhead)
          else if ( j == 1 .and. i /= 1 ) then
             work_array(i, j) = A(ihead, jhead) - A(ihead - nw0, jhead)
          else if ( i == 1 .and. j /= 1 ) then
             work_array(i, j) = A(ihead, jhead) - A(ihead, jhead - nw1)
          else
             work_array(i, j) = A(ihead, jhead) - A(ihead - nw0, jhead) - A(ihead, jhead - nw1) + A(ihead - nw0, jhead - nw1)
          end if
       end do
    end do
    print *, "do u understand the invesion of the order of matrix shape?"

    hist2d%len0 =  N1
    hist2d%len1 =  N0
    hist2d%arr = C_loc(work_array)
  end function hist2d


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


  subroutine delete_array2(arr_length0, arr_length1, array) bind(C, name="delete_array2")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array2
    integer(c_int), intent(in) :: arr_length0
    integer(c_int), intent(in) :: arr_length1
    type(c_ptr), value :: array
    real(c_double), pointer :: work_array(:,:)

    call C_F_pointer(array, work_array, [arr_length0, arr_length1])
    deallocate(work_array)
    array = C_NULL_PTR
    print *, "Is work_array freed ?, work_array: ",work_array
  end subroutine delete_array2

end module histfort
