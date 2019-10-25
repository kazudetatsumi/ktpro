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
  integer datasize(4)
  integer maxw(4)

contains

  function cost4d(maxw3, maxw2, maxw1, maxw0, datasize3, datasize2, datasize1, datasize0,& 
      usecond, data_array, condition) bind(C, name="cost4d")
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: maxw2
    integer(c_int), intent(in) :: maxw3
    integer(c_int), intent(in) :: datasize0
    integer(c_int), intent(in) :: datasize1
    integer(c_int), intent(in) :: datasize2
    integer(c_int), intent(in) :: datasize3
    real(c_double), intent(in) :: data_array(datasize0, datasize1, datasize2, datasize3)                         
    integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    logical(c_bool), intent(in) :: usecond
    type(result) :: cost4d                                  
    double precision :: cumdata(datasize0, datasize1, datasize2, datasize3)
    integer :: cum_cond(datasize0, datasize1, datasize2, datasize3)
    double precision, allocatable :: hist_array(:,:,:,:)                    
    integer, allocatable :: hist_cond(:,:,:,:)
    integer nw(4)
    integer histsize(4)
    integer ax_id1, ax_id2, ax_id3, width_id1, width_id2, width_id3, width_id4
    datasize = (/datasize0, datasize1, datasize2, datasize3/)
    cumdata = cumsum4d(cumsum4d(cumsum4d(cumsum4d(data_array, 1), 2), 3), 4)
    if (usecond) cum_cond = cumsum4di(cumsum4di(cumsum4di(cumsum4di(condition, 1), 2), 3), 4)
    maxw = (/maxw0, maxw1, maxw2, maxw3/)
    allocate(cost(maxw(1), maxw(2), maxw(3), maxw(4)))
    allocate(histaves(maxw(1), maxw(2), maxw(3), maxw(4)))
    allocate(deltas(maxw(1), maxw(2), maxw(3), maxw(4)))
    cost(:,:,:,:) = 0.0
    histaves(:,:,:,:) = 0.0
    deltas(:,:,:,:) = 0.0
    !call help0d(data_array, condition, usecond)
    !do ax_id1 = 1, 4
    !  call help1d(data_array, ax_id1, condition, usecond)
    !enddo
    !do ax_id1 = 1, 3
    !do ax_id2 = ax_id1 + 1, 4
    !  call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
    !enddo
    !enddo
    !do ax_id1 = 1, 2
    !do ax_id2 = ax_id1 + 1, 3
    !do ax_id3 = ax_id2 + 1, 4
    !  call help3d(data_array, [ax_id1, ax_id2, ax_id3], condition, usecond)
    !enddo
    !enddo
    !enddo


    !do width_id1 = 1, maxw(1)
    do width_id1 = 2, 11
    nw(1) = width_id1
    histsize(1) = (datasize(1) - mod(datasize(1), nw(1))) / nw(1)
    do width_id2 = 1, maxw(2)
    nw(2) = width_id2
    histsize(2) = (datasize(2) - mod(datasize(2), nw(2))) / nw(2)
    do width_id3 = 1, maxw(3)
    nw(3) = width_id3
    histsize(3) = (datasize(3) - mod(datasize(3), nw(3))) / nw(3)
    do width_id4 = 1, maxw(4)
    nw(4) = width_id4
    histsize(4) = (datasize(4) - mod(datasize(4), nw(4))) / nw(4)
       if (width_id1 /= 1 .and. width_id2 /= 1 .and. width_id3 /=1 .and. width_id4 /=1) then
          hist_array = hist4d(cumdata, histsize, nw)
          if (usecond) then 
              hist_cond = hist4di(cum_cond, histsize, nw)
              call stat1(pack(hist_array, hist_cond == maxval(hist_cond)), nw)
          else
            call stat(hist_array, nw) 
          endif
       end if
    end do
    end do
    end do
    end do


    print *, "minloc cost:", minloc(cost), "with its value:", minval(cost)
    cost4d%len0 =  maxw(4)
    cost4d%len1 =  maxw(3)
    cost4d%len2 =  maxw(2)
    cost4d%len3 =  maxw(1)
    cost4d%arr = C_loc(cost)
    cost4d%kavearr = C_loc(histaves)
    cost4d%darr = C_loc(deltas)
  end function cost4d

  subroutine help0d(data_array, condition, usecond)
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    double precision, allocatable :: hist_array(:,:,:,:)
    integer nw(4)
    nw = 1
    hist_array = data_array
    if (usecond) then
      call stat1(pack(hist_array, condition == maxval(condition)), nw)
    else
      call stat(hist_array, nw)
    end if
  end subroutine help0d
  subroutine help1d(data_array, ax, condition, usecond)
    integer, intent(in) :: ax
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    logical(1) :: Ishistsizeallocated
    double precision, allocatable :: cumdata1d(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)
    integer, allocatable :: histcond(:,:,:,:)
    integer width_id, histsize, nw(4),  cumdata_size(4), cumdata_order(4)
    nw = 1
    Ishistsizeallocated = .true.
    cumdata1d = cumsum4d(data_array, ax)
    if (usecond) then
        cum_cond = cumsum4di(condition, ax)
    endif
    if (ax == 1) then
      Ishistsizeallocated = .false.
    elseif (ax == 2) then
      cumdata_size = (/datasize(2), datasize(1), datasize(3), datasize(4)/) 
      cumdata_order = (/2, 1, 3, 4/)
    elseif (ax == 3) then
      cumdata_size = (/datasize(3), datasize(2), datasize(1), datasize(4)/) 
      cumdata_order = (/3, 2, 1, 4/)
    elseif (ax == 4) then
      cumdata_size = (/datasize(4), datasize(2), datasize(3), datasize(1)/) 
      cumdata_order = (/4, 2, 3, 1/)
    endif
    if (Ishistsizeallocated) then
      cumdata1d = reshape(cumdata1d, cumdata_size, order = cumdata_order)
      if (usecond) then
        cum_cond = reshape(cum_cond, cumdata_size, order = cumdata_order)
      endif
    endif
    do width_id = 2, maxw(ax)
      nw(ax) = width_id
      histsize = (datasize(ax) - mod(datasize(ax), nw(ax))) / nw(ax)
      hist_array = hist1d(cumdata1d, histsize, nw(ax))
      if (usecond) then
        histcond = hist1di(cum_cond, histsize, nw(ax))
        call stat1(pack(hist_array, histcond == maxval(histcond)), nw)
      else
        call stat(hist_array, nw)
      end if
    enddo
  end subroutine help1d

  subroutine help2d(data_array, ax, condition, usecond)
    integer, intent(in) :: ax(:)
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    logical(1) :: Ishistsizeallocated
    double precision, allocatable :: cumdata2d(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)
    integer, allocatable :: histcond(:,:,:,:)
    integer nw(4), histsize(4), cumdata_size_id, width_id1, width_id2
    integer, allocatable ::  cumdata_size(:,:), cumdata_order(:,:)
    nw = 1
    cumdata2d = cumsum4d(cumsum4d(data_array, ax(1)), ax(2))
    if (usecond) cum_cond = cumsum4di(cumsum4di(condition, ax(1)), ax(2))
    if (ax(1) == 1 .and. ax(2) == 2) then
      Ishistsizeallocated = .false.
    else
      Ishistsizeallocated = .true.
      if (ax(1) == 1 .and. ax(2) == 3) then
        allocate(cumdata_size(1,4), cumdata_order(1,4))
        cumdata_size(1,:) = (/datasize(1), datasize(3), datasize(2), datasize(4)/)
        cumdata_order(1,:) = (/1, 3, 2, 4/)
      elseif (ax(1) == 1 .and. ax(2) == 4) then
        allocate(cumdata_size(1,4), cumdata_order(1,4))
        cumdata_size(1,:) = (/datasize(1), datasize(4), datasize(3), datasize(2)/)
        cumdata_order(1,:) = (/1, 4, 3, 2/)
      elseif (ax(1) == 2 .and. ax(2) == 3) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(2), datasize(1), datasize(3), datasize(4)/)
        cumdata_order(1,:) = (/2, 1, 3, 4/)
        cumdata_size(2,:) = (/datasize(2), datasize(3), datasize(1), datasize(4)/)
        cumdata_order(2,:) = (/1, 3, 2, 4/)
      elseif (ax(1) == 2 .and. ax(2) == 4) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(2), datasize(1), datasize(3), datasize(4)/)
        cumdata_order(1,:) = (/2, 1, 3, 4/)
        cumdata_size(2,:) = (/datasize(2), datasize(4), datasize(3), datasize(1)/)
        cumdata_order(2,:) = (/1, 4, 3, 2/)
      elseif (ax(1) == 3 .and. ax(2) == 4) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(3), datasize(2), datasize(1), datasize(4)/)
        cumdata_order(1,:) = (/3, 2, 1, 4/)
        cumdata_size(2,:) = (/datasize(3), datasize(4), datasize(1), datasize(2)/)
        cumdata_order(2,:) = (/1, 4, 3, 2/)
      endif 
    endif
    if (Ishistsizeallocated) then
      do cumdata_size_id = 1, size(cumdata_size,1)
         cumdata2d = reshape(cumdata2d, cumdata_size(cumdata_size_id,1:4), order = cumdata_order(cumdata_size_id,1:4))
         if (usecond) cum_cond = reshape(cum_cond, cumdata_size(cumdata_size_id,1:4), &
                                         order = cumdata_order(cumdata_size_id, 1:4))
      enddo
    endif
    do width_id1 = 2, maxw(ax(1))
      nw(ax(1)) = width_id1
    do width_id2 = 2, maxw(ax(2))
      nw(ax(2)) = width_id2
      histsize = (datasize - mod(datasize, nw)) / nw
      hist_array = hist2d(cumdata2d, [histsize(ax(1)), histsize(ax(2))], [nw(ax(1)), nw(ax(2))])
      if (usecond) then
        histcond = hist2di(cum_cond, [histsize(ax(1)), histsize(ax(2))], [nw(ax(1)), nw(ax(2))])
        call stat1(pack(hist_array, histcond == maxval(histcond)), nw)
      else
        call stat(hist_array, nw)
      endif
    enddo
    enddo
  end subroutine help2d

  subroutine help3d(data_array, ax, condition, usecond)
    integer, intent(in) :: ax(:)
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    logical(1) :: Ishistsizeallocated
    double precision, allocatable :: cumdata3d(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)
    integer, allocatable :: histcond(:,:,:,:)
    integer nw(4), histsize(4), cumdata_size_id, width_id1, width_id2, width_id3
    integer, allocatable ::  cumdata_size(:,:),cumdata_order(:,:)
    nw = 1
    cumdata3d = cumsum4d(cumsum4d(cumsum4d(data_array, ax(1)), ax(2)), ax(3))
    if (usecond) cum_cond = cumsum4di(cumsum4di(cumsum4di(condition, ax(1)), ax(2)), ax(3))
    if (ax(1) == 1 .and. ax(2) == 2 .and. ax(3) == 3) then
      Ishistsizeallocated = .false.
    else
      Ishistsizeallocated = .true.
      if (ax(1) == 1 .and. ax(2) == 2 .and. ax(3) == 4) then
        allocate(cumdata_size(1,4), cumdata_order(1,4))
        cumdata_size(1,:) = (/datasize(1), datasize(2), datasize(4), datasize(3)/)
        cumdata_order(1,:) = (/1, 2, 4, 3/)
      elseif (ax(1) == 1 .and. ax(2) == 3 .and. ax(3) == 4) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(1), datasize(3), datasize(2), datasize(4)/)
        cumdata_order(1,:) = (/1, 3, 2, 4/)
        cumdata_size(2,:) = (/datasize(1), datasize(3), datasize(4), datasize(2)/)
        cumdata_order(2,:) = (/1, 2, 4, 3/)
      elseif (ax(1) == 2 .and. ax(2) == 3 .and. ax(3) == 4) then
        allocate(cumdata_size(3,4), cumdata_order(3,4))
        cumdata_size(1,:) = (/datasize(2), datasize(1), datasize(3), datasize(4)/)
        cumdata_order(1,:) = (/2, 1, 3, 4/)
        cumdata_size(2,:) = (/datasize(2), datasize(3), datasize(1), datasize(4)/)
        cumdata_order(2,:) = (/1, 3, 2, 4/)
        cumdata_size(3,:) = (/datasize(2), datasize(3), datasize(4), datasize(1)/)
        cumdata_order(3,:) = (/1, 2, 4, 3/)
      end if
    endif
    if (Ishistsizeallocated) then
      do cumdata_size_id = 1, size(cumdata_size,1)
        cumdata3d = reshape(cumdata3d, cumdata_size(cumdata_size_id,1:4), order = cumdata_order(cumdata_size_id,1:4))
        if (usecond) cum_cond = reshape(cum_cond, cumdata_size(cumdata_size_id,1:4), order = cumdata_order(cumdata_size_id,1:4))
      enddo
    endif
    do width_id1 = 2, maxw(ax(1))
      nw(ax(1)) = width_id1
    do width_id2 = 2, maxw(ax(2))
      nw(ax(2)) = width_id2
    do width_id3 = 2, maxw(ax(3))
      nw(ax(3)) = width_id3
      histsize = (datasize - mod(datasize, nw)) / nw
      hist_array = hist3d(cumdata3d, [histsize(ax(1)), histsize(ax(2)), histsize(ax(3))], [nw(ax(1)), nw(ax(2)), nw(ax(3))])
      if (usecond) then
        histcond = hist3di(cum_cond, [histsize(ax(1)), histsize(ax(2)), histsize(ax(3))], [nw(ax(1)), nw(ax(2)), nw(ax(3))])
        call stat1(pack(hist_array, histcond == maxval(histcond)), nw)
      else
        call stat(hist_array, nw) 
      endif
    enddo
    enddo
    enddo
  end subroutine help3d

  subroutine stat(k, nw)
    double precision, intent(in) :: k(:,:,:,:)
    integer, intent(in) :: nw(4)
    double precision kave, v
    kave = sum(k) / dble(size(k)); v = sum((k - kave)**2) / dble(size(k))
    print *, "cost with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / (dble(product(nw))**2)
    histaves(nw(1), nw(2), nw(3), nw(4)) = kave
    deltas(nw(1), nw(2), nw(3), nw(4)) = product(nw)
    cost(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / (dble(product(nw))**2)
  end subroutine stat
  subroutine stat1(knonzero, nw)
    double precision, intent(in) :: knonzero(:)
    integer, intent(in) :: nw(4)
    double precision kave, v
    kave = sum(knonzero) / real(size(knonzero)); v = sum((knonzero - kave)**2) / real(size(knonzero))
    !print *, "cost with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / (real(product(nw))**2)
    histaves(nw(1), nw(2), nw(3), nw(4)) = kave
    deltas(nw(1), nw(2), nw(3), nw(4)) = product(nw)
    cost(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / (real(product(nw))**2)
  end subroutine stat1
          
  subroutine delete_array_pointer(L3, L2, L1, L0, carray, karray, darray) bind(C, name="delete_array_pointer")
    integer(c_int), intent(in) :: L3, L2, L1, L0  
    type(c_ptr), value :: carray, karray, darray
    real(c_double), pointer :: cost(:,:,:,:), histaves(:,:,:,:), deltas(:,:,:,:)
    call C_F_pointer(carray, cost, [L0, L1, L2, L3])
    call C_F_pointer(karray, histaves, [L0, L1, L2, L3])
    call C_F_pointer(darray, deltas, [L0, L1, L2, L3])
    deallocate(cost,histaves,deltas)
    carray = C_NULL_PTR; karray = C_NULL_PTR; darray = C_NULL_PTR
  end subroutine delete_array_pointer


  function cumsum4d(data_array, ax)
    integer, intent(in) :: ax
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    double precision :: cumsum4d(datasize(1), datasize(2), datasize(3),datasize(4))
    integer :: cum_id
    cumsum4d = data_array
    if (ax == 4) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(:, :, :, cum_id + 1) = cumsum4d(:, :, :, cum_id) + data_array(:, :, :, cum_id + 1)
       end do
    else if (ax == 3) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(:, :, cum_id + 1, :) = cumsum4d(:, :, cum_id, :) + data_array(:, :, cum_id + 1, :)
       end do
    else if (ax == 2) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(:, cum_id + 1, :, :) = cumsum4d(:, cum_id, :, :) + data_array(:, cum_id + 1, :, :)
       end do
    else if (ax == 1) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(cum_id + 1, :, :, :) = cumsum4d(cum_id, :, :, :) + data_array(cum_id + 1, :, :, :)
       end do
    end if
  end function cumsum4d       
  function cumsum4di(d, ax)
    integer, intent(in) :: ax
    integer, intent(in) :: d(datasize(1), datasize(2), datasize(3), datasize(4))
    integer :: cumsum4di(datasize(1), datasize(2), datasize(3), datasize(4))
    integer :: cum_id
    cumsum4di = d
    if (ax == 4) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(:, :, :, cum_id + 1) = cumsum4di(:, :, :, cum_id) + d(:, :, :, cum_id + 1)
       end do
    else if (ax == 3) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(:, :, cum_id + 1, :) = cumsum4di(:, :, cum_id, :) + d(:, :, cum_id + 1, :)
       end do
    else if (ax == 2) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(:, cum_id + 1, :, :) = cumsum4di(:, cum_id, :, :) + d(:, cum_id + 1, :, :)
       end do
    else if (ax == 1) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(cum_id + 1, :, :, :) = cumsum4di(cum_id, :, :, :) + d(cum_id + 1, :, :, :)
       end do
    end if
  end function cumsum4di

  function hist4d(cumdata, N, nw)
    double precision, intent(in) :: cumdata(:,:,:,:) 
    integer, intent(in) :: N(4), nw(4)
    integer :: i, j, h, l, ihead, jhead, hhead, lhead
    double precision :: hist4d(N(1), N(2), N(3),N(4))
    !$omp parallel do
    do i = 1, N(1)
    ihead = i*nw(1) 
    do j = 1, N(2)
    jhead = j*nw(2) 
    do h = 1, N(3)
    hhead = h*nw(3) 
    do l = 1, N(4)
    lhead = l*nw(4) 
      if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead)
      else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead - nw(1), jhead, hhead, lhead)
      else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead)
      else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead)
      else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4))
      else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead)
      else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead)
      else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4))
      else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead)
      else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4))
      else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4))
      else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead)
      else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4))
      else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4))
      else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4))
      else
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                       + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead - nw(4))
      end if
    end do
    end do
    end do
    end do
  end function hist4d

  function hist4di(cumdata, N, nw)
    integer, intent(in) :: cumdata(:,:,:,:) 
    integer, intent(in) :: N(4), nw(4)
    integer :: i, j, h, l, ihead, jhead, hhead, lhead
    integer :: hist4di(N(1),N(2),N(3),N(4))
    !$omp parallel do
    do i = 1, N(1)
    ihead = i*nw(1) 
    do j = 1, N(2)
    jhead = j*nw(2) 
    do h = 1, N(3)
    hhead = h*nw(3) 
    do l = 1, N(4)
    lhead = l*nw(4) 
       if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead)
       else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead - nw(1), jhead, hhead, lhead)
       else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead)
       else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead)
       else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4))
       else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead)
       else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead)
       else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4))
       else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead)
       else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4))
       else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4))
       else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead)
       else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4))
       else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4))
       else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4))
       else
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead - nw(4))
       end if
    end do
    end do
    end do
    end do
  end function hist4di

  function hist1d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N, nw
    integer :: i, ihead
    double precision :: hist1d(N, size(B,2), size(B,3), size(B,4))
    hist1d(1,:,:,:) = B(nw,:,:,:)
    !$omp parallel do
    do i = 2, N
       ihead = i*nw
       hist1d(i,:,:,:) = B(ihead,:,:,:) - B(ihead - nw,:,:,:)
    end do
  end function hist1d
  function hist1di(C, N, nw)
    integer, intent(in) :: C(:,:,:,:) 
    integer, intent(in) :: N, nw
    integer :: i, ihead
    integer :: hist1di(N, size(C,2), size(C,3), size(C,4))
    hist1di(1,:,:,:) = C(nw,:,:,:)
    !$omp parallel do
    do i = 2, N
       ihead = i*nw
       hist1di(i,:,:,:) = C(ihead,:,:,:) - C(ihead - nw,:,:,:)
    end do
  end function hist1di

  function hist2d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(2), nw(2)
    integer :: i, ihead, j, jhead
    double precision :: hist2d(N(1), N(2), size(B,3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2) 
          if ( i == 1 .and. j == 1) then
              hist2d(i, j, :, :) = B(ihead, jhead,:, :)
          else if ( i /= 1 .and. j == 1) then
              hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :)
          else if ( i == 1 .and. j /= 1) then
              hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead, jhead - nw(2), :, :)
          else 
              hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :) - B(ihead, jhead - nw(2), :, :) &
                                   + B(ihead - nw(1), jhead - nw(2), :, :)
          end if
       end do
    end do
  end function hist2d
  function hist2di(B, N, nw)
    integer, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(2), nw(2)
    integer :: i, ihead, j, jhead
    integer :: hist2di(N(1), N(2), size(B,3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2) 
          if ( i == 1 .and. j == 1) then
              hist2di(i, j, :, :) = B(ihead, jhead,:, :)
          else if ( i /= 1 .and. j == 1) then
              hist2di(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :)
          else if ( i == 1 .and. j /= 1) then
              hist2di(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead, jhead - nw(2), :, :)
          else 
              hist2di(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :) - B(ihead, jhead - nw(2), :, :) &
                                   + B(ihead - nw(1), jhead - nw(2), :, :)
          end if
       end do
    end do
  end function hist2di

  function hist3d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(3), nw(3)
    integer :: i, ihead, j, jhead, h, hhead
    double precision :: hist3d(N(1), N(2), N(3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2)
          do h = 1, N(3)
             hhead = h*nw(3)
             if ( i == 1 .and. j == 1 .and. h == 1) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :)
             else if ( j == 1 .and. i /= 1 .and. h == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :)
             else if ( i == 1 .and. j /= 1 .and. h == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :)
             else if ( i == 1 .and. h /= 1 .and. j == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead, hhead - nw(3), :)
             else if ( i /= 1 .and. j /= 1 .and. h == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                     - B(ihead, jhead - nw(2), hhead, :) + B(ihead - nw(1), jhead - nw(2), hhead, :)
             else if ( i /= 1 .and. j == 1 .and. h /= 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :)
             else if ( i == 1 .and. j /= 1 .and. h /= 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead, jhead - nw(2), hhead - nw(3), :)
             else
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead - nw(2), hhead, :) - B(ihead, jhead, hhead - nw(3), :) &
                                    + B(ihead, jhead - nw(2), hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :) &
                                    + B(ihead - nw(1), jhead - nw(2), hhead, :) - B(ihead - nw(1), jhead - nw(2), hhead - nw(3), :)
             end if
          end do
       end do
    end do
  end function hist3d
  function hist3di(B, N, nw)
    integer, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(3), nw(3)
    integer :: i, ihead, j, jhead, h, hhead
    integer :: hist3di(N(1), N(2), N(3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2)
          do h = 1, N(3)
             hhead = h*nw(3)
             if ( i == 1 .and. j == 1 .and. h == 1) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :)
             else if ( j == 1 .and. i /= 1 .and. h == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :)
             else if ( i == 1 .and. j /= 1 .and. h == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :)
             else if ( i == 1 .and. h /= 1 .and. j == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead, hhead - nw(3), :)
             else if ( i /= 1 .and. j /= 1 .and. h == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                     - B(ihead, jhead - nw(2), hhead, :) + B(ihead - nw(1), jhead - nw(2), hhead, :)
             else if ( i /= 1 .and. j == 1 .and. h /= 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :)
             else if ( i == 1 .and. j /= 1 .and. h /= 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead, jhead - nw(2), hhead - nw(3), :)
             else
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead - nw(2), hhead, :) - B(ihead, jhead, hhead - nw(3), :) &
                                    + B(ihead, jhead - nw(2), hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :) &
                                    + B(ihead - nw(1), jhead - nw(2), hhead, :) - B(ihead - nw(1), jhead - nw(2), hhead - nw(3), :)
             end if
          end do
       end do
    end do
  end function hist3di

end module costfort4d
