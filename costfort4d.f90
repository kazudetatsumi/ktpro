module costfort4d
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    integer(c_int) :: len1
    integer(c_int) :: len2
    integer(c_int) :: len3
    type(c_ptr) :: arr
    type(c_ptr) :: kavearr
    type(c_ptr) :: darr
  end type result

contains

  !function cost4d(maxw3, maxw2, maxw1, maxw0, Al3, Al2, Al1, Al0, AA, D, CD) bind(C, name="cost4d")
  function cost4d(maxw3, maxw2, maxw1, maxw0, Al3, Al2, Al1, Al0, D, CD) bind(C, name="cost4d")
    !DEC$ ATTRIBUTES DLLEXPORT :: cost4d
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: maxw2
    integer(c_int), intent(in) :: maxw3
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    integer(c_int), intent(in) :: Al2
    integer(c_int), intent(in) :: Al3
    !real(c_double), intent(in) :: AA(Al0, Al1, Al2, Al3)                         
    real(c_double), intent(in) :: D(Al0, Al1, Al2, Al3)                         
    integer(c_int), intent(in) :: CD(Al0, Al1, Al2, Al3)                         
    type(result) :: cost4d                                  
    !real(c_double), pointer :: k(:,:,:,:)                    
    double precision, allocatable :: k(:,:,:,:)                    
    integer, allocatable :: kcond(:,:,:,:)                    
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    integer nw0, nw1, nw2, nw3
    integer N0, N1, N2, N3
    integer i, ihead, j, jhead, h, hhead, l, lhead
    double precision kave, v
    double precision A(Al0, Al1, Al2, Al3)
    !allocate(k(Al0, Al1, Al2, Al3))
    A = cumsum4d(cumsum4d(cumsum4d(cumsum4d(D, 1),2),3),4)
    
    allocate(Cn(maxw0, maxw1, maxw2, maxw3))
    allocate(kaves(maxw0, maxw1, maxw2, maxw3))
    allocate(deltas(maxw0, maxw1, maxw2, maxw3))
    Cn(:,:,:,:) = 0.0
    kaves(:,:,:,:) = 0.0
    deltas(:,:,:,:) = 0.0
    do nw0 = 1, maxw0
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    do nw1 = 1, maxw1
    N1 = (Al1 - mod(Al1, nw1)) / nw1
    do nw2 = 1, maxw2
    N2 = (Al2 - mod(Al2, nw2)) / nw2
    do nw3 = 1, maxw3
    N3 = (Al3 - mod(Al3, nw3)) / nw3
    allocate(k(N0, N1, N2, N3))
    allocate(kcond(N0, N1, N2, N3))
       if (nw0*nw1*nw2*nw3 == 1) then
       k = D
       kcond = CD
       else
       !$omp parallel do
       do i = 1, N0
       ihead = i*nw0 
       do j = 1, N1
       jhead = j*nw1 
       do h = 1, N2
       hhead = h*nw2 
       do l = 1, N3
       lhead = l*nw3 
          !k(i, j, h, l) = sum(D(ihead-nw0+1:ihead, jhead-nw1+1:jhead,  hhead-nw2+1:hhead, lhead-nw3+1:lhead))
          if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead)
          else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead - nw0, jhead, hhead, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead)
          else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead)
          else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead - nw0, jhead, hhead - nw2, lhead)
          else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead, hhead, lhead - nw3)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead, lhead - nw3)
          else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw2, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3)
          else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead - nw2, lhead) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead) - A(ihead - nw0, jhead - nw1, hhead - nw2, lhead)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                           - A(ihead - nw0, jhead - nw1, hhead, lhead - nw3)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                           - A(ihead - nw0, jhead, hhead - nw2, lhead - nw3)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                           - A(ihead, jhead - nw1, hhead - nw2, lhead - nw3)
          else
             k(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
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
       end if
       kave = sum(k) / real(N0*N1*N2*N3)
       v = sum((k - kave)**2) / real(N0*N1*N2*N3)
       print *, "Cn with ", nw0, nw1, nw2, nw3, ":", (2.0 * kave - v) / (real(nw0*nw1*nw2*nw3)**2)
       Cn(nw0, nw1, nw2, nw3) = (2.0 * kave - v) / (real(nw0*nw1*nw2*nw3)**2)
       deltas(nw0, nw1, nw2, nw3) = real(nw0*nw1*nw2*nw3)
       kaves(nw0, nw1, nw2, nw3) = kave
       deallocate(k)
       deallocate(kcond)
    end do
    end do
    end do
    end do

    !deallocate(k)

    print *, "minloc Cn:", minloc(Cn)
    cost4d%len0 =  maxw3
    cost4d%len1 =  maxw2
    cost4d%len2 =  maxw1
    cost4d%len3 =  maxw0
    cost4d%arr = C_loc(Cn(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%kavearr = C_loc(kaves(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%darr = C_loc(deltas(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
  end function cost4d


  subroutine delete_array(arr_length3, arr_length2, arr_length1, arr_length0, carray, karray, darray) bind(C, name="delete_array")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array
    integer(c_int), intent(in) :: arr_length0
    integer(c_int), intent(in) :: arr_length1
    integer(c_int), intent(in) :: arr_length2
    integer(c_int), intent(in) :: arr_length3
    type(c_ptr), value :: carray
    type(c_ptr), value :: karray
    type(c_ptr), value :: darray
    real(c_double), pointer :: Cn(:,:,:,:)
    real(c_double), pointer :: kaves(:,:,:,:)
    real(c_double), pointer :: deltas(:,:,:,:)

    call C_F_pointer(carray, Cn, [arr_length0, arr_length1, arr_length2, arr_length3])
    call C_F_pointer(karray, kaves, [arr_length0, arr_length1, arr_length2, arr_length3])
    call C_F_pointer(darray, deltas, [arr_length0, arr_length1, arr_length2, arr_length3])
    deallocate(Cn)
    deallocate(kaves)
    deallocate(deltas)
    carray = C_NULL_PTR
    karray = C_NULL_PTR
    darray = C_NULL_PTR
  end subroutine delete_array

  function cumsum4d(d, axis)
      integer, intent(in) :: axis
      double precision, intent(in) :: d(:,:,:,:)
      double precision :: cumsum4d(size(d,1), size(d,2), size(d,3), size(d,4))
      integer :: i, j, k, l 
      cumsum4d = d

      if (axis == 4) then
         !$omp parallel do
         do i = 1, size(d,1)
         do j = 1, size(d,2)
         do k = 1, size(d,3)
         do l = 1, size(d,4) - 1
            cumsum4d(i, j, k, l + 1) = cumsum4d(i, j, k, l) + d(i, j, k, l + 1)
         end do
         end do
         end do
         end do
      else if (axis == 3) then
         !$omp parallel do
         do i = 1, size(d,1)
         do j = 1, size(d,2)
         do l = 1, size(d,4)
         do k = 1, size(d,3) - 1
            cumsum4d(i, j, k + 1, l) = cumsum4d(i, j, k, l) + d(i, j, k + 1, l)
         end do
         end do
         end do
         end do
      else if (axis == 2) then
         !$omp parallel do
         do i = 1, size(d,1)
         do k = 1, size(d,3)
         do l = 1, size(d,4)
         do j = 1, size(d,2) - 1
            cumsum4d(i, j + 1, k, l) = cumsum4d(i, j, k, l) + d(i, j + 1, k, l)
         end do
         end do
         end do
         end do
      else if (axis == 1) then
         !$omp parallel do
         do j = 1, size(d,2) 
         do k = 1, size(d,3)
         do l = 1, size(d,4)
         do i = 1, size(d,1) - 1
            cumsum4d(i + 1, j, k, l) = cumsum4d(i, j, k, l) + d(i + 1, j, k, l)
         end do
         end do
         end do
         end do
      end if
   end function cumsum4d       





end module costfort4d
