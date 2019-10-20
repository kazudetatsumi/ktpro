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

  !function cost4d(maxw3, maxw2, maxw1, maxw0, Al3, Al2, Al1, Al0, A, D, CDA) bind(C, name="cost4d")
  !function cost4d(maxw3, maxw2, maxw1, maxw0, Al3, Al2, Al1, Al0, usecond,  A, D, condition) bind(C, name="cost4d")
  function cost4d(maxw3, maxw2, maxw1, maxw0, Al3, Al2, Al1, Al0, usecond, D, condition) bind(C, name="cost4d")
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: maxw2
    integer(c_int), intent(in) :: maxw3
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    integer(c_int), intent(in) :: Al2
    integer(c_int), intent(in) :: Al3
    !real(c_double), intent(in) :: A(Al0, Al1, Al2, Al3)                         
    real(c_double), intent(in) :: D(Al0, Al1, Al2, Al3)                         
    !integer(c_int), intent(in) :: CDA(Al0, Al1, Al2, Al3)                         
    logical(c_bool), intent(in) :: condition(Al0, Al1, Al2, Al3)                         
    logical(c_bool), intent(in) :: usecond
    type(result) :: cost4d                                  
    double precision :: A(Al0, Al1, Al2, Al3)
    double precision, allocatable :: k(:,:,:,:)                    
    logical, allocatable :: kcond(:,:,:,:)
    !integer, allocatable :: kcond(:,:,:,:)                    
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    integer nw(4),maxw(4)
    integer N(4)
    integer i, j, h, l
    integer :: Al(4)
    !double precision, allocatable :: knonzero(:)
    A = cumsum4d(cumsum4d(cumsum4d(cumsum4d(d, 1),2),3),4)
    Al = (/Al0, Al1, Al2, Al3/)
    maxw = (/maxw0, maxw1, maxw2, maxw3/)
   

    
    allocate(Cn(maxw(1), maxw(2), maxw(3), maxw(4)))
    allocate(kaves(maxw(1), maxw(2), maxw(3), maxw(4)))
    allocate(deltas(maxw(1), maxw(2), maxw(3), maxw(4)))
    Cn(:,:,:,:) = 0.0
    kaves(:,:,:,:) = 0.0
    deltas(:,:,:,:) = 0.0



    call help1d(D, 0, 1, condition, Al, usecond, Cn, kaves, deltas)
    do i=1,4
      call help1d(D, i, maxw(i), condition, Al, usecond, Cn, kaves, deltas)
    enddo
    do i=1,3
    do j=i+1,4
      call help2d(D, [i, j], maxw, condition, Al, usecond, Cn, kaves, deltas)
    enddo
    enddo
    !call help3d(D, 1, maxw1, maxw2, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    !call help3d(D, 2, maxw0, maxw2, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    !call help3d(D, 3, maxw0, maxw1, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    !call help3d(D, 4, maxw0, maxw1, maxw2, Al0, Al1, Al2, Al3, Cn, kaves, deltas)


    do i = 1, maxw(1)
    nw(1) = i
    N(1) = (Al(1) - mod(Al(1), nw(1))) / nw(1)
    do j = 1, maxw(2)
    nw(2) = j
    N(2) = (Al(2) - mod(Al(2), nw(2))) / nw(2)
    do h = 1, maxw(3)
    nw(3) = h
    N(3) = (Al(3) - mod(Al(3), nw(3))) / nw(3)
    do l = 1, maxw(4)
    nw(4) = l
    N(4) = (Al(4) - mod(Al(4), nw(4))) / nw(4)
       if (i /= 1 .and. j /= 1 .and. h /=1 .and. l /=1) then
          k = hist4d(A, N, nw)
       !kcond = hist4di(CDA, N(1), N(2), N(3), N(4), nw(1), nw(2), nw(3), nw(4))
          !! now !! kcond = hist4dcond(condition, N(1), N(2), N(3), N(4), nw(1), nw(2), nw(3), nw(4))
       !knonzero = pack(k, kcond == maxval(kcond))
         call stat(k, Cn, kaves, deltas, nw) 
       end if
    end do
    end do
    end do
    end do


    print *, "minloc Cn:", minloc(Cn)
    cost4d%len0 =  maxw(4)
    cost4d%len1 =  maxw(3)
    cost4d%len2 =  maxw(2)
    cost4d%len3 =  maxw(1)
    cost4d%arr = C_loc(Cn)
    cost4d%kavearr = C_loc(kaves)
    cost4d%darr = C_loc(deltas)
  end function cost4d

  subroutine help1d(D, axis, mxw, condition, Al, usecond,  Cn, kaves, deltas)
    double precision, intent(in) :: D(:,:,:,:)
    integer, intent(in) :: axis, mxw, Al(4)
    logical(1), intent(in) :: condition(Al(1), Al(2), Al(3), Al(4))
    logical(1), allocatable :: cd2(:,:,:,:)
    logical(1), intent(in) :: usecond
    double precision, allocatable :: B(:,:,:,:)
    double precision, allocatable :: k(:,:,:,:)
    integer i, N, nw(4),  bshape(4), border(4)
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    nw = 1
    if (axis == 0) then
        k = D
        if (usecond) then
          call stat1(pack(k, condition), Cn, kaves, deltas, nw)
        else
          call stat(k, Cn, kaves, deltas, nw)
        end if
    else
      B = cumsum4d(d, axis)
      if (axis == 1) then
        bshape = (/Al(1), Al(2), Al(3), Al(4)/) 
        border = (/1, 2, 3, 4/)
      elseif (axis == 2) then
        bshape = (/Al(2), Al(1), Al(3), Al(4)/) 
        border = (/2, 1, 3, 4/)
      elseif (axis == 3) then
        bshape = (/Al(3), Al(2), Al(1), Al(4)/) 
        border = (/3, 2, 1, 4/)
      elseif (axis == 4) then
        bshape = (/Al(4), Al(2), Al(3), Al(1)/) 
        border = (/4, 2, 3, 1/)
      endif
      B = reshape(B, bshape, order = border)
      cd2 = reshape(condition, bshape, order = border)
      do i = 2, mxw
        nw(axis) = i
        N = (Al(axis) - mod(Al(axis), nw(axis))) / nw(axis)
        k = hist1d(B, N, nw(axis))
        if (usecond) then
          call stat1(pack(k, mask1d(cd2, N, nw(axis))),  Cn, kaves, deltas, nw)
        else
          call stat(k, Cn, kaves, deltas, nw)
        end if
      enddo
    end if
 end subroutine help1d

  subroutine help2d(D, ax, maxw,  condition, Al, usecond, Cn, kaves, deltas)
    double precision, intent(in) :: D(:,:,:,:)
    integer, intent(in) :: ax(:), maxw(:), Al(:)
    logical(1), intent(in) :: condition(:,:,:,:)
    logical(1), intent(in) :: usecond
    double precision, allocatable :: B(:,:,:,:)
    double precision, allocatable :: k(:,:,:,:)
    integer i, j, nw(4), N(4)
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    nw = 1
    B = cumsum4d(cumsum4d(d, ax(1)), ax(2))
    if (ax(1) == 1 .and. ax(2) == 3) then
      B = reshape(B, (/Al(1), Al(3), Al(2), Al(4)/), order = (/1, 3, 2, 4/))
    elseif (ax(1) == 1 .and. ax(2) == 4) then
      B = reshape(B, (/Al(1), Al(4), Al(3), Al(2)/), order = (/1, 4, 3, 2/))
    elseif (ax(1) == 2 .and. ax(2) == 3) then
      B = reshape(B, (/Al(2), Al(1), Al(3), Al(4)/), order = (/2, 1, 3, 4/))
      B = reshape(B, (/Al(2), Al(3), Al(1), Al(4)/), order = (/1, 3, 2, 4/))
    elseif (ax(1) == 2 .and. ax(2) == 4) then
      B = reshape(B, (/Al(2), Al(1), Al(3), Al(4)/), order = (/2, 1, 3, 4/))
      B = reshape(B, (/Al(2), Al(4), Al(3), Al(1)/), order = (/1, 4, 3, 2/))
    elseif (ax(1) == 3 .and. ax(2) == 4) then
      B = reshape(B, (/Al(3), Al(2), Al(1), Al(4)/), order = (/3, 2, 1, 4/))
      B = reshape(B, (/Al(3), Al(4), Al(1), Al(2)/), order = (/1, 4, 3, 2/))
    endif 
    do i = 2, maxw(ax(1))
      nw(ax(1)) = i
    do j = 2, maxw(ax(2))
      nw(ax(2)) = j
      N = (Al - mod(Al, nw)) / nw
      k = hist2d(B, [N(ax(1)), N(ax(2))], [nw(ax(1)), nw(ax(2))])
      call stat(k, Cn, kaves, deltas, nw)
    enddo
    enddo
  end subroutine help2d


  subroutine help3d(D, axis, maxw1, maxw2, maxw3, Al, Cn, kaves, deltas)
    double precision, intent(in) :: D(:,:,:,:)
    integer, intent(in) :: axis, maxw1, maxw2, maxw3, Al(4)
    double precision, allocatable :: B(:,:,:,:)
    double precision, allocatable :: k(:,:,:,:)
    integer nw(4), N(4), i, j, l
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    nw = 1
    if (axis == 4) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 1),2),3)
      do i = 2, maxw1
        nw(1) = i
      do j = 2, maxw2
        nw(2) = j
      do l = 2, maxw3
        nw(3) = l
        N = (Al - mod(Al, nw)) / nw
        k = hist3d(B, [N(1), N(2), N(3)], [nw(1), nw(2), nw(3)])
        call stat(k, Cn, kaves, deltas, nw) 
      enddo
      enddo
      enddo
    elseif (axis == 3) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 1),2),4)
      B = reshape(B, (/size(B,1), size(B,2), size(B,4), size(B,3)/), order = (/1, 2, 4, 3/))
      do i = 2, maxw1
        nw(1) = i
      do j = 2, maxw2
        nw(2) = j
      do l = 2, maxw3
        nw(4) = l
        N = (Al - mod(Al, nw)) / nw
        k = hist3d(B, [N(1), N(2), N(4)], [nw(1), nw(2), nw(4)]); call stat(k, Cn, kaves, deltas, nw) 
      enddo
      enddo
      enddo
    elseif (axis == 2) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 1),3),4)
      B = reshape(B, (/size(B,1), size(B,3), size(B,2), size(B,4)/), order = (/1, 3, 2, 4/))
      B = reshape(B, (/size(B,1), size(B,2), size(B,4), size(B,3)/), order = (/1, 2, 4, 3/))
      do i = 2, maxw1
        nw(1) = i
      do j = 2, maxw2
        nw(3) = j
      do l = 2, maxw3
        nw(4) = l
        N = (Al - mod(Al, nw)) / nw
        k = hist3d(B, [N(1), N(3), N(4)], [nw(1), nw(3), nw(4)]); call stat(k, Cn, kaves, deltas, nw) 
      enddo
      enddo
      enddo
    elseif (axis == 1) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 2),3),4)
      B = reshape(B, (/size(B,2), size(B,1), size(B,3), size(B,4)/), order = (/2, 1, 3, 4/))
      B = reshape(B, (/size(B,1), size(B,3), size(B,2), size(B,4)/), order = (/1, 3, 2, 4/))
      B = reshape(B, (/size(B,1), size(B,2), size(B,4), size(B,3)/), order = (/1, 2, 4, 3/))
      do i = 2, maxw1
        nw(2) = i
      do j = 2, maxw2
        nw(3) = j
      do l = 2, maxw3
        nw(4) = l
        N = (Al - mod(Al, nw)) / nw
        k = hist3d(B, [N(2), N(3), N(4)], [nw(2), nw(3), nw(4)]); call stat(k, Cn, kaves, deltas, nw) 
      enddo
      enddo
      enddo
    end if
  end subroutine help3d

  subroutine stat(k, Cn, kaves, deltas, nw)
    double precision, intent(in) :: k(:,:,:,:)
    integer, intent(in) :: nw(4)
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    double precision kave, v
    kave = sum(k) / real(size(k)); v = sum((k - kave)**2) / real(size(k))
    print *, "Cn with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / real(product(nw)**2)
    kaves(nw(1), nw(2), nw(3), nw(4)) = kave
    deltas(nw(1), nw(2), nw(3), nw(4)) = product(nw)
    Cn(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / real(product(nw)**2)
  end subroutine stat
  subroutine stat1(knonzero, Cn, kaves, deltas, nw)
    double precision, intent(in) :: knonzero(:)
    integer, intent(in) :: nw(4)
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    double precision kave, v
    kave = sum(knonzero) / real(size(knonzero)); v = sum((knonzero - kave)**2) / real(size(knonzero))
    print *, "Cn with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / real(product(nw)**2)
    kaves(nw(1), nw(2), nw(3), nw(4)) = kave
    deltas(nw(1), nw(2), nw(3), nw(4)) = product(nw)
    Cn(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / real(product(nw)**2)
  end subroutine stat1
         
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


  function cumsum4d(d, axis)
    integer, intent(in) :: axis
    double precision, intent(in) :: d(:,:,:,:)
    double precision :: cumsum4d(size(d,1), size(d,2), size(d,3), size(d,4))
    integer :: i, j, k, l, N(4)
    cumsum4d = d
    N(1) = size(d,1); N(2) = size(d,2); N(3) = size(d,3); N(4) = size(d,4)


    if (axis == 4) then
       !!$omp parallel do
       !do i = 1, N(1)
       !do j = 1, N(3)
       !do k = 1, N(4)
       do l = 1, N(4) - 1
          !cumsum4d(i, j, k, l + 1) = cumsum4d(i, j, k, l) + d(i, j, k, l + 1)
          cumsum4d(:, :, :, l + 1) = cumsum4d(:, :, :, l) + d(:, :, :, l + 1)
       end do
       !end do
       !end do
       !end do
    else if (axis == 3) then
       !!$omp parallel do
       !do i = 1, N(1)
       !do j = 1, N(2)
       !do l = 1, N(4)
       do k = 1, N(3) - 1
          !cumsum4d(i, j, k + 1, l) = cumsum4d(i, j, k, l) + d(i, j, k + 1, l)
          cumsum4d(:, :, k + 1, :) = cumsum4d(:, :, k, :) + d(:, :, k + 1, :)
       end do
       !end do
       !end do
       !end do
    else if (axis == 2) then
       !!$omp parallel do
       !do i = 1, N(1)
       !do k = 1, N(3)
       !do l = 1, N(4)
       do j = 1, N(2) - 1
          !cumsum4d(i, j + 1, k, l) = cumsum4d(i, j, k, l) + d(i, j + 1, k, l)
          cumsum4d(:, j + 1, :, :) = cumsum4d(:, j, :, :) + d(:, j + 1, :, :)
       end do
       !end do
       !end do
       !end do
    else if (axis == 1) then
       !!$omp parallel do
       !do j = 1, N(2) 
       !do k = 1, N(3)
       !do l = 1, N(4)
       do i = 1, N(1) - 1
          !cumsum4d(i + 1, j, k, l) = cumsum4d(i, j, k, l) + d(i + 1, j, k, l)
          cumsum4d(i + 1, :, :, :) = cumsum4d(i, :, :, :) + d(i + 1, :, :, :)
       end do
       !end do
       !end do
       !end do
    end if
  end function cumsum4d       

  function hist4d(A, N, nw)
    double precision, intent(in) :: A(:,:,:,:) 
    integer, intent(in) :: N(4), nw(4)
    integer :: i, j, h, l, ihead, jhead, hhead, lhead
    double precision :: hist4d(N(1), N(2), N(3),N(4))
    do i = 1, N(1)
    ihead = i*nw(1) 
    do j = 1, N(2)
    jhead = j*nw(2) 
    do h = 1, N(3)
    hhead = h*nw(3) 
    !$omp parallel do
    do l = 1, N(4)
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

  function hist4di(A, N, nw)
    integer, intent(in) :: A(:,:,:,:) 
    integer, intent(in) :: N(4), nw(4)
    integer :: i, j, h, l, ihead, jhead, hhead, lhead
    integer :: hist4di(N(1),N(2),N(3),N(4))
    do i = 1, N(1)
    ihead = i*nw(1) 
    do j = 1, N(2)
    jhead = j*nw(2) 
    do h = 1, N(3)
    hhead = h*nw(3) 
    !$omp parallel do
    do l = 1, N(4)
    lhead = l*nw(4) 
       if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead)
       else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead - nw(1), jhead, hhead, lhead)
       else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead)
       else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead)
       else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw(4))
       else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                        + A(ihead - nw(1), jhead - nw(2), hhead, lhead)
       else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                        + A(ihead - nw(1), jhead, hhead - nw(3), lhead)
       else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                        + A(ihead - nw(1), jhead, hhead, lhead - nw(4))
       else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead, jhead - nw(2), hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                        + A(ihead, jhead - nw(2), hhead - nw(3), lhead)
       else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead, jhead - nw(2), hhead, lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                        + A(ihead, jhead - nw(2), hhead, lhead - nw(4))
       else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead, jhead, hhead - nw(3), lhead) - A(ihead, jhead, hhead, lhead - nw(4)) &
                        + A(ihead, jhead, hhead - nw(3), lhead - nw(4))
       else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                        - A(ihead, jhead, hhead - nw(3), lhead) &
                        + A(ihead, jhead - nw(2), hhead - nw(3), lhead) + A(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + A(ihead - nw(1), jhead - nw(2), hhead, lhead) - A(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead)
       else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead - nw(2), hhead, lhead) &
                        - A(ihead, jhead, hhead, lhead - nw(4)) &
                        + A(ihead - nw(1), jhead - nw(2), hhead, lhead) + A(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + A(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        - A(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4))
       else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead - nw(1), jhead, hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                        - A(ihead, jhead, hhead, lhead - nw(4)) &
                        + A(ihead - nw(1), jhead, hhead - nw(3), lhead) + A(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + A(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - A(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4))
       else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                        - A(ihead, jhead - nw(2), hhead, lhead) - A(ihead, jhead, hhead - nw(3), lhead) &
                        - A(ihead, jhead, hhead, lhead - nw(4)) &
                        + A(ihead, jhead - nw(2), hhead - nw(3), lhead) + A(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        + A(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - A(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4))
       else
          hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
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

  function mask1d(cd, N, nw)
    logical(1), intent(in) :: CD(:,:,:,:) 
    integer, intent(in) :: N, nw
    integer :: i, is, ie
    logical(1) :: mask1d(N, size(CD,2), size(CD,3), size(CD,4))
    do i = 1, N
       is = (i-1)*nw + 1
       ie = i*nw
       mask1d(i,:,:,:) = all(CD(is:ie,:,:,:),1)
    end do
  end function mask1d

  function hist2d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(2), nw(2)
    integer :: i, ihead, j, jhead
    double precision :: hist2d(N(1), N(2), size(B,3), size(B,4))
    do i = 1, N(1)
       ihead = i*nw(1) 
       !$omp parallel do
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

  function hist3d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(3), nw(3)
    integer :: i, ihead, j, jhead, h, hhead
    double precision :: hist3d(N(1), N(2), N(3), size(B,4))
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

end module costfort4d
