module histfort
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    integer(c_int) :: len1
    integer(c_int) :: len2
    integer(c_int) :: len3
    type(c_ptr) :: arr
  end type result

contains


  function hist4d(nw3, nw2, nw1, nw0, Al3, Al2, Al1, Al0, A) bind(C, name="hist4d")
    !DEC$ ATTRIBUTES DLLEXPORT :: hist3d
    integer(c_int), intent(in) :: nw0
    integer(c_int), intent(in) :: nw1
    integer(c_int), intent(in) :: nw2
    integer(c_int), intent(in) :: nw3
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    integer(c_int), intent(in) :: Al2
    integer(c_int), intent(in) :: Al3
    real(c_double), intent(in) :: A(Al0, Al1, Al2, Al3)                         
    type(result) :: hist4d                                  
    real(c_double), pointer :: work_array(:, :, :, :)                    
    integer N0, N1, N2, N3
    integer i, ihead, j, jhead, h, hhead, l, lhead
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    N1 = (Al1 - mod(Al1, nw1)) / nw1
    N2 = (Al2 - mod(Al2, nw2)) / nw2
    N3 = (Al3 - mod(Al3, nw3)) / nw3
    allocate(work_array(N0, N1, N2, N3))

    do i = 1, N0
       ihead = i*nw0 
       do j = 1, N1
          jhead = j*nw1
          do h = 1, N2
             hhead = h*nw2
             do l = 1, N3
                lhead = l*nw3
                if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead)
                else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead - nw0, jhead, hhead, lhead)
                else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead)
                else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead)
                else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3)
                else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                                          + A(ihead - nw0, jhead - nw1, hhead, lhead)
                else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                                          + A(ihead - nw0, jhead, hhead - nw2, lhead)
                else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead - nw0, jhead, hhead, lhead - nw3)
                else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                                          + A(ihead, jhead - nw1, hhead - nw2, lhead)
                else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead, jhead - nw1, hhead, lhead - nw3)
                else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead, jhead, hhead - nw2, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead, jhead, hhead - nw2, lhead - nw3)

                else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                                          - A(ihead, jhead, hhead - nw2, lhead) &
                                          + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead - nw2, lhead) &
                                          + A(ihead - nw0, jhead - nw1, hhead, lhead) &
                                          - A(ihead - nw0, jhead - nw1, hhead - nw2, lhead)
                else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                                          - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead - nw0, jhead - nw1, hhead, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                                          + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                                          - A(ihead - nw0, jhead - nw1, hhead, lhead - nw3)
                else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                                          - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead - nw0, jhead, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                                          + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                                          - A(ihead - nw0, jhead, hhead - nw2, lhead - nw3)
                else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                                          - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                                          + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                                          - A(ihead, jhead - nw1, hhead - nw2, lhead - nw3)
                else
                   work_array(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                                          - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                                          - A(ihead, jhead, hhead - nw2, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                                          + A(ihead - nw0, jhead - nw1, hhead, lhead) + A(ihead - nw0, jhead, hhead - nw2, lhead) &
                                          + A(ihead - nw0, jhead, hhead, lhead - nw3) + A(ihead, jhead - nw1, hhead - nw2, lhead) &
                                          + A(ihead, jhead - nw1, hhead, lhead - nw3) + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                                          - A(ihead, jhead - nw1, hhead - nw2, lhead - nw3) &
                                          - A(ihead - nw0, jhead, hhead - nw2, lhead - nw3) &
                                          - A(ihead - nw0, jhead - nw1, hhead, lhead - nw3) &
                                          - A(ihead - nw0, jhead - nw1, hhead - nw2, lhead) &
                                          + A(ihead - nw0, jhead - nw1, hhead - nw2, lhead - nw3)
                end if
             end do
          end do
       end do
    end do

    hist4d%len0 =  N3
    hist4d%len1 =  N2
    hist4d%len2 =  N1
    hist4d%len3 =  N0
    hist4d%arr = C_loc(work_array)
  end function hist4d



  subroutine delete_array4(arr_length0, arr_length1, arr_length2, arr_length3, array) bind(C, name="delete_array4")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array4
    integer(c_int), intent(in) :: arr_length0
    integer(c_int), intent(in) :: arr_length1
    integer(c_int), intent(in) :: arr_length2
    integer(c_int), intent(in) :: arr_length3
    type(c_ptr), value :: array
    real(c_double), pointer :: work_array(:,:,:,:)

    call C_F_pointer(array, work_array, [arr_length0, arr_length1, arr_length2, arr_length3])
    deallocate(work_array)
    array = C_NULL_PTR
    !The line below causes a segfault with gfortran v9.4.0.
    !print *, "Is work_array freed ?, work_array: ", work_array
  end subroutine delete_array4

end module histfort
