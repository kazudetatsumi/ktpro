module costfort3d
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    integer(c_int) :: len1
    integer(c_int) :: len2
    type(c_ptr) :: arr
    type(c_ptr) :: kavearr
    type(c_ptr) :: darr
  end type result

contains

  function cost3d(maxw2, maxw1, maxw0, Al2, Al1, Al0, A) bind(C, name="cost3d")
    !DEC$ ATTRIBUTES DLLEXPORT :: cost3d
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: maxw2
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    integer(c_int), intent(in) :: Al2
    real(c_double), intent(in) :: A(Al0, Al1, Al2)                         
    type(result) :: cost3d                                  
    real(c_double), pointer :: k(:,:,:)                    
    real(c_double), pointer :: Cn(:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:)                    
    integer nw0, nw1, nw2
    integer N0, N1, N2
    integer i, ihead, j, jhead, h, hhead
    real kave,v
    allocate(k(Al0, Al1, Al2))
    allocate(Cn(maxw0, maxw1, maxw2))
    allocate(kaves(maxw0, maxw1, maxw2))
    allocate(deltas(maxw0, maxw1, maxw2))
    Cn(:,:,:) = 0.0
    kaves(:,:,:) = 0.0
    deltas(:,:,:) = 0.0
    do nw0 = 1, maxw0
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    do nw1 = 1, maxw1
    N1 = (Al1 - mod(Al1, nw1)) / nw1
    do nw2 = 1, maxw2
    N2 = (Al2 - mod(Al2, nw2)) / nw2
       do i = 1, N0
       ihead = i*nw0 
       do j = 1, N1
       jhead = j*nw1 
       do h = 1, N2
       hhead = h*nw2 
          if ( i == 1 .and. j == 1 .and. h == 1) then
             k(i, j, h) = A(ihead, jhead, hhead)
          else if ( j == 1 .and. i /= 1 .and. h == 1 ) then
             k(i, j, h) = A(ihead, jhead, hhead) - A(ihead - nw0, jhead, hhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 ) then
             k(i, j, h) = A(ihead, jhead, hhead) - A(ihead, jhead - nw1, hhead)
          else if ( i == 1 .and. h /= 1 .and. j == 1 ) then
             k(i, j, h) = A(ihead, jhead, hhead) - A(ihead, jhead, hhead - nw2)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 ) then
             k(i, j, h) = A(ihead, jhead, hhead) - A(ihead - nw0, jhead, hhead) - A(ihead, jhead - nw1, hhead) &
                        + A(ihead - nw0, jhead - nw1, hhead)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 ) then
             k(i, j, h) = A(ihead, jhead, hhead) - A(ihead - nw0, jhead, hhead) - A(ihead, jhead, hhead - nw2) &
                        + A(ihead - nw0, jhead, hhead - nw2)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 ) then
             k(i, j, h) = A(ihead, jhead, hhead) - A(ihead, jhead - nw1, hhead) - A(ihead, jhead, hhead - nw2) &
                        + A(ihead, jhead - nw1, hhead - nw2)
          else
             k(i, j, h) = A(ihead, jhead, hhead) &
                        - A(ihead - nw0, jhead, hhead) - A(ihead, jhead - nw1, hhead) - A(ihead, jhead, hhead - nw2) &
                        + A(ihead, jhead - nw1, hhead - nw2) + A(ihead - nw0, jhead, hhead - nw2) &
                        + A(ihead - nw0, jhead - nw1, hhead) - A(ihead - nw0, jhead - nw1, hhead - nw2)
          end if
       end do
       end do
       end do
       kave = sum(k(1:N0, 1:N1, 1:N2)) / real(N0*N1*N2)
       v = sum((k(1:N0, 1:N1, 1:N2) - kave)**2) / real(N0*N1*N2)
       !print *, "cost for ", nw0, nw1, ":", (2.0 * kave - v) / (real(nw0*nw1)**2)
       Cn(nw0, nw1, nw2) = (2.0 * kave - v) / (real(nw0*nw1*nw2)**2)
       deltas(nw0, nw1, nw2) = real(nw0*nw1*nw2)
       kaves(nw0, nw1, nw2) = kave
    end do
    end do
    end do

    print *, "minloc Cn:", minloc(Cn)
    cost3d%len0 =  maxw2
    cost3d%len1 =  maxw1
    cost3d%len2 =  maxw0
    cost3d%arr = C_loc(Cn(1:maxw0, 1:maxw1, 1:maxw2))
    cost3d%kavearr = C_loc(kaves(1:maxw0, 1:maxw1, 1:maxw2))
    cost3d%darr = C_loc(deltas(1:maxw0, 1:maxw1, 1:maxw2))
  end function cost3d

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

end module costfort3d
