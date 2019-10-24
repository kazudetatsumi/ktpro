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
  real(c_double), pointer :: cost(:,:,:,:)                    
  real(c_double), pointer :: histaves(:,:,:,:)                    
  real(c_double), pointer :: deltas(:,:,:,:)                    

contains

  function cost4d(maxw3, maxw2, maxw1, maxw0, datasize3, datasize2, datasize1, datasize0, A, data_array, condition) & 
          bind(C, name="cost4d")
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: maxw2
    integer(c_int), intent(in) :: maxw3
    integer(c_int), intent(in) :: datasize0
    integer(c_int), intent(in) :: datasize1
    integer(c_int), intent(in) :: datasize2
    integer(c_int), intent(in) :: datasize3
    real(c_double), intent(in) :: A(datasize0, datasize1, datasize2, datasize3)                         
    real(c_double), intent(in) :: data_array(datasize0, datasize1, datasize2, datasize3)                         
    integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    type(result) :: cost4d                                  
    double precision, allocatable :: k(:,:,:,:)                    
    double precision :: cumdata(datasize0, datasize1, datasize2, datasize3)                    
    integer width_id1, width_id2, width_id3, width_id4
    integer nw(4), histsize(4), datasize(4), maxw(4)
    integer i, ihead, j, jhead, h, hhead, l, lhead
    datasize = (/datasize0, datasize1, datasize2, datasize3/)
    maxw = (/maxw0, maxw1, maxw2, maxw3/)
    allocate(cost(maxw0, maxw1, maxw2, maxw3))
    allocate(histaves(maxw0, maxw1, maxw2, maxw3))
    allocate(deltas(maxw0, maxw1, maxw2, maxw3))
    cost(:,:,:,:) = 0.0
    histaves(:,:,:,:) = 0.0
    deltas(:,:,:,:) = 0.0

    cumdata = cumsum4d(cumsum4d(cumsum4d(cumsum4d(data_array, 1, datasize),2, datasize),3, datasize),4, datasize)

    do width_id1 = 14, maxw(1)
    nw(1) = width_id1
    histsize(1) = (datasize(1) - mod(datasize(1), nw(1))) / nw(1)
    do width_id2 = 14, maxw(2)
    nw(2) = width_id2
    histsize(2) = (datasize(2) - mod(datasize(2), nw(2))) / nw(2)
    do width_id3 = 14, maxw(3)
    nw(3) = width_id3
    histsize(3) = (datasize(3) - mod(datasize(3), nw(3))) / nw(3)
    do width_id4 = 14, maxw(4)
    nw(4) = width_id4
    histsize(4) = (datasize(4) - mod(datasize(4), nw(4))) / nw(4)
       if (width_id1 /= 1 .and. width_id2 /= 1 .and. width_id3 /=1 .and. width_id4 /=1) then
          k = hist4d(cumdata, histsize, nw)
         call stat(k, nw) 
       end if
    end do
    end do
    end do
    end do


    print *, "minloc Cn:", minloc(cost(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%len0 =  maxw3
    cost4d%len1 =  maxw2
    cost4d%len2 =  maxw1
    cost4d%len3 =  maxw0
    cost4d%arr = C_loc(cost(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%kavearr = C_loc(histaves(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%darr = C_loc(deltas(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
  end function cost4d




  subroutine stat(k, nw)
    double precision, allocatable :: k(:,:,:,:)
    integer, intent(in) :: nw(:)
    double precision kave, v
       kave = sum(k) / real(size(k)); v = sum((k - kave)**2) / real(size(k))
       print *, "Cn with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / (real(product(nw))**2)
       histaves(nw(1), nw(2), nw(3), nw(4)) = kave
       deltas(nw(1), nw(2), nw(3), nw(4)) = real(product(nw))
       cost(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / (real(product(nw))**2)
       deallocate(k)
  end subroutine stat

  subroutine delete_array_pointer(L3, L2, L1, L0, carray, karray, darray) bind(C, name="delete_array_pointer")
    integer(c_int), intent(in) :: L3, L2, L1, L0  
    type(c_ptr), value :: carray, karray, darray
    real(c_double), pointer :: Cn(:,:,:,:), kaves(:,:,:,:), deltas(:,:,:,:)
    call C_F_pointer(carray, Cn, [L0, L1, L2, L3])
    call C_F_pointer(karray, kaves, [L0, L1, L2, L3])
    call C_F_pointer(darray, deltas, [L0, L1, L2, L3])
    deallocate(Cn,kaves,deltas)
    carray = C_NULL_PTR; karray = C_NULL_PTR; darray = C_NULL_PTR
  end subroutine delete_array_pointer


  function cumsum4d(data_array, axis, datasize)
      integer, intent(in) :: axis, datasize(4)
      double precision, intent(in) :: data_array(:,:,:,:)
      double precision :: cumsum4d(datasize(1), datasize(2), datasize(3), datasize(4))
      integer :: cum_id
      cumsum4d = data_array
      if (axis == 4) then
         do cum_id = 1, datasize(axis) - 1
            cumsum4d(:, :, :, cum_id + 1) = cumsum4d(:, :, :, cum_id) + data_array(:, :, :, cum_id + 1)
         end do
      else if (axis == 3) then
         do cum_id = 1, datasize(axis) - 1
            cumsum4d(:, :, cum_id + 1, :) = cumsum4d(:, :, cum_id, :) + data_array(:, :, cum_id + 1, :)
         end do
      else if (axis == 2) then
         do cum_id = 1, datasize(axis) - 1
            cumsum4d(:, cum_id + 1, :, :) = cumsum4d(:, cum_id, :, :) + data_array(:, cum_id + 1, :, :)
         end do
      else if (axis == 1) then
         do cum_id = 1, datasize(axis) - 1
            cumsum4d(cum_id + 1, :, :, :) = cumsum4d(cum_id, :, :, :) + data_array(cum_id + 1, :, :, :)
         end do
      end if
   end function cumsum4d       


   function hist4d(A, histsize, nw) 
       double precision, intent(in) :: A(:,:,:,:)
       integer, intent(in) :: histsize(:), nw(:)
       integer :: i, j, h, l, ihead, jhead, hhead, lhead
       double precision :: hist4d(histsize(1), histsize(2), histsize(3), histsize(4))
       do i = 1, histsize(1)
       ihead = i*nw(1)
       do j = 1, histsize(2)
       jhead = j*nw(2)
       do h = 1, histsize(3)
       hhead = h*nw(3)
       !$omp parallel do
       do l = 1, histsize(4)
       lhead = l*nw(4)
          if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead)
          else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead - nw(1), jhead, hhead, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead)
          else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead)
          else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw(4))
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                           + A(ihead - nw(1), jhead - nw(2), hhead, lhead)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                           + A(ihead - nw(1), jhead, hhead - nw(3), lhead)
          else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead - nw(1), jhead, hhead, lhead - nw(4))
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw(2), hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                           + A(ihead, jhead - nw(2), hhead - nw(3), lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw(2), hhead, lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead, jhead - nw(2), hhead, lhead - nw(4))
          else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw(3), lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead, jhead, hhead - nw(3), lhead - nw(4))
          else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                           - A(ihead, jhead, hhead - nw(3), lhead) &
                           + A(ihead, jhead - nw(2), hhead - nw(3), lhead) + A(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                           + A(ihead - nw(1), jhead - nw(2), hhead, lhead) - A(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead - nw(1), jhead - nw(2), hhead, lhead) + A(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                           + A(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                           - A(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4))
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                           - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead - nw(1), jhead, hhead - nw(3), lhead) + A(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                           + A(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                           - A(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4))
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw(2), hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                           - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead, jhead - nw(2), hhead - nw(3), lhead) + A(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                           + A(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                           - A(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4))
          else
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                           - A(ihead, jhead, hhead - nw(3), lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                           + A(ihead - nw(1), jhead - nw(2), hhead, lhead) + A(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                           + A(ihead - nw(1), jhead, hhead, lhead - nw(4)) + A(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                           + A(ihead, jhead - nw(2), hhead, lhead - nw(4)) + A(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                           - A(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4)) &
                           - A(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4)) &
                           - A(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4)) &
                           - A(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead) &
                           + A(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead - nw(4))
          end if
       end do
       end do
       end do
       end do
  end function hist4d


end module costfort4d
