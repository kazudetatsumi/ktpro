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
  function cost4d(maxw3, maxw2, maxw1, maxw0, Al3, Al2, Al1, Al0, A, D, condition) bind(C, name="cost4d")
    integer(c_int), intent(in) :: maxw0
    integer(c_int), intent(in) :: maxw1
    integer(c_int), intent(in) :: maxw2
    integer(c_int), intent(in) :: maxw3
    integer(c_int), intent(in) :: Al0
    integer(c_int), intent(in) :: Al1
    integer(c_int), intent(in) :: Al2
    integer(c_int), intent(in) :: Al3
    real(c_double), intent(in) :: A(Al0, Al1, Al2, Al3)                         
    real(c_double), intent(in) :: D(Al0, Al1, Al2, Al3)                         
    !integer(c_int), intent(in) :: CDA(Al0, Al1, Al2, Al3)                         
    logical(c_bool), intent(in) :: condition(Al0, Al1, Al2, Al3)                         
    type(result) :: cost4d                                  
    double precision, allocatable :: k(:,:,:,:)                    
    logical, allocatable :: kcond(:,:,:,:)
    !integer, allocatable :: kcond(:,:,:,:)                    
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    integer nw0, nw1, nw2, nw3
    integer N0, N1, N2, N3
    integer i, ihead, j, jhead, h, hhead, l, lhead
    !double precision, allocatable :: knonzero(:)

    
    allocate(Cn(maxw0, maxw1, maxw2, maxw3))
    allocate(kaves(maxw0, maxw1, maxw2, maxw3))
    allocate(deltas(maxw0, maxw1, maxw2, maxw3))
    Cn(:,:,:,:) = 0.0
    kaves(:,:,:,:) = 0.0
    deltas(:,:,:,:) = 0.0



    call help1d(D, 0, 1, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help1d(D, 1, maxw0, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help1d(D, 2, maxw1, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help1d(D, 3, maxw2, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help1d(D, 4, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help2d(D, 1, 2, maxw0, maxw1, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help2d(D, 1, 3, maxw0, maxw2, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help2d(D, 1, 4, maxw0, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help2d(D, 2, 3, maxw1, maxw2, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help2d(D, 2, 4, maxw1, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help2d(D, 3, 4, maxw2, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help3d(D, 1, maxw1, maxw2, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help3d(D, 2, maxw0, maxw2, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help3d(D, 3, maxw0, maxw1, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    call help3d(D, 4, maxw0, maxw1, maxw2, Al0, Al1, Al2, Al3, Cn, kaves, deltas)


    do nw0 = 1, maxw0
    N0 = (Al0 - mod(Al0, nw0)) / nw0
    do nw1 = 1, maxw1
    N1 = (Al1 - mod(Al1, nw1)) / nw1
    do nw2 = 1, maxw2
    N2 = (Al2 - mod(Al2, nw2)) / nw2
    do nw3 = 1, maxw3
    N3 = (Al3 - mod(Al3, nw3)) / nw3
       if (nw0 /= 1 .and. nw1 /= 1 .and. nw2 /=1 .and. nw3 /=1) then
          k = hist4d(A, N0, N1, N2, N3, nw0, nw1, nw2, nw3)
       !kcond = hist4di(CDA, N0, N1, N2, N3, nw0, nw1, nw2, nw3)
          !! now !! kcond = hist4dcond(condition, N0, N1, N2, N3, nw0, nw1, nw2, nw3)
       !knonzero = pack(k, kcond == maxval(kcond))
         call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
       end if
    end do
    end do
    end do
    end do


    print *, "minloc Cn:", minloc(Cn(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%len0 =  maxw3
    cost4d%len1 =  maxw2
    cost4d%len2 =  maxw1
    cost4d%len3 =  maxw0
    cost4d%arr = C_loc(Cn(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%kavearr = C_loc(kaves(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
    cost4d%darr = C_loc(deltas(1:maxw0, 1:maxw1, 1:maxw2, 1:maxw3))
  end function cost4d

  subroutine help1d(D, axis, maxw, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    double precision, intent(in) :: D(:,:,:,:)
    integer, intent(in) :: axis, maxw, Al0, Al1, Al2, Al3
    double precision, allocatable :: B(:,:,:,:)
    double precision, allocatable :: k(:,:,:,:)
    integer nw0, nw1, nw2, nw3, N0, N1, N2, N3
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    if (axis == 0) then
        nw0 = 1; nw1 = 1; nw2 = 1; nw3 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k=d; call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
    elseif (axis == 1) then
      B = cumsum4d(d, axis)
      do nw0 = 2, maxw
        nw1 = 1; nw2 = 1; nw3 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        k = hist1d(B, N0, nw0);  call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
      enddo
    elseif (axis == 2) then
      B = cumsum4d(d, axis)
      B = reshape(B, (/size(B,2), size(B,1), size(B,3), size(B,4)/), order = (/2, 1, 3, 4/))
      do nw1 = 2, maxw
        nw0 = 1; nw2 = 1; nw3 = 1
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        k = hist1d(B, N1, nw1);  call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
    elseif (axis == 3) then
      B = cumsum4d(d, axis)
      B = reshape(B, (/size(B,3), size(B,2), size(B,1), size(B,4)/), order = (/3, 2, 1, 4/))
      do nw2 = 2, maxw
        nw0 = 1; nw1 = 1; nw3 = 1
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        k = hist1d(B, N2, nw2); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
    elseif (axis == 4) then
      B = cumsum4d(d, axis)
      B = reshape(B, (/size(B,4), size(B,2), size(B,3), size(B,1)/), order = (/4, 2, 3, 1/))
      do nw3 = 2, maxw
        nw0 = 1; nw1 = 1; nw2 = 1
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist1d(B, N3, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
      enddo
    endif
  end subroutine help1d
    

  subroutine help2d(D, axis1, axis2, maxw1, maxw2, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    double precision, intent(in) :: D(:,:,:,:)
    integer, intent(in) :: axis1, axis2, maxw1, maxw2, Al0, Al1, Al2, Al3
    double precision, allocatable :: B(:,:,:,:)
    double precision, allocatable :: k(:,:,:,:)
    integer nw0, nw1, nw2, nw3, N0, N1, N2, N3
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    if (axis1 == 1 .and. axis2 == 2) then
      B = cumsum4d(cumsum4d(d, axis1), axis2)
      do nw0 = 2, maxw1
      do nw1 = 2, maxw2
        nw2 = 1; nw3 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        k = hist2d(B, N0, N1, nw0, nw1); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
      enddo
    elseif (axis1 == 1 .and. axis2 == 3) then
      B = cumsum4d(cumsum4d(d, axis1), axis2)
      B = reshape(B, (/size(B,1), size(B,3), size(B,2), size(B,4)/), order = (/1, 3, 2, 4/))
      do nw0 = 2, maxw1
      do nw2 = 2, maxw2
        nw1 = 1; nw3 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        k = hist2d(B, N0, N2, nw0, nw2); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
      enddo
    elseif (axis1 == 1 .and. axis2 == 4) then
      B = cumsum4d(cumsum4d(d, axis1), axis2)
      B = reshape(B, (/size(B,1), size(B,4), size(B,3), size(B,2)/), order = (/1, 4, 3, 2/))
      do nw0 = 2, maxw1
      do nw3 = 2, maxw2
        nw1 = 1; nw2 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist2d(B, N0, N3, nw0, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
      enddo
    elseif (axis1 == 2 .and. axis2 == 3) then
      B = cumsum4d(cumsum4d(d, axis1), axis2)
      B = reshape(B, (/size(B,2), size(B,1), size(B,3), size(B,4)/), order = (/2, 1, 3, 4/))
      B = reshape(B, (/size(B,1), size(B,3), size(B,2), size(B,4)/), order = (/1, 3, 2, 4/))
      do nw1 = 2, maxw1
      do nw2 = 2, maxw2
        nw0 = 1; nw3 = 1
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        k = hist2d(B, N1, N2, nw1, nw2); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
      enddo
    elseif (axis1 == 2 .and. axis2 == 4) then
      B = cumsum4d(cumsum4d(d, axis1), axis2)
      B = reshape(B, (/size(B,2), size(B,1), size(B,3), size(B,4)/), order = (/2, 1, 3, 4/))
      B = reshape(B, (/size(B,1), size(B,4), size(B,3), size(B,2)/), order = (/1, 4, 3, 2/))
      do nw1 = 2, maxw1
      do nw3 = 2, maxw2
        nw0 = 1; nw2 = 1
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist2d(B, N1, N3, nw1, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
      enddo
    elseif (axis1 == 3 .and. axis2 == 4) then
      B = cumsum4d(cumsum4d(d, axis1), axis2)
      B = reshape(B, (/size(B,3), size(B,2), size(B,1), size(B,4)/), order = (/3, 2, 1, 4/))
      B = reshape(B, (/size(B,1), size(B,4), size(B,3), size(B,2)/), order = (/1, 4, 3, 2/))
      do nw2 = 2, maxw1
      do nw3 = 2, maxw2
        nw0 = 1; nw1 = 1
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist2d(B, N2, N3, nw2, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
      enddo
      enddo
      deallocate(B)
    endif 
  end subroutine help2d

  subroutine help3d(D, axis, maxw1, maxw2, maxw3, Al0, Al1, Al2, Al3, Cn, kaves, deltas)
    double precision, intent(in) :: D(:,:,:,:)
    integer, intent(in) :: axis, maxw1, maxw2, maxw3, Al0, Al1, Al2, Al3
    double precision, allocatable :: B(:,:,:,:)
    double precision, allocatable :: k(:,:,:,:)
    integer nw0, nw1, nw2, nw3, N0, N1, N2, N3
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    if (axis == 4) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 1),2),3)
      do nw0 = 2, maxw1
      do nw1 = 2, maxw2
      do nw2 = 2, maxw3
        nw3 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        k = hist3d(B, N0, N1, N2, nw0, nw1, nw2); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
      enddo
      enddo
      enddo
    elseif (axis == 3) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 1),2),4)
      B = reshape(B, (/size(B,1), size(B,2), size(B,4), size(B,3)/), order = (/1, 2, 4, 3/))
      do nw0 = 2, maxw1
      do nw1 = 2, maxw2
      do nw3 = 2, maxw3
        nw2 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist3d(B, N0, N1, N3, nw0, nw1, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
      enddo
      enddo
      enddo
    elseif (axis == 2) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 1),3),4)
      B = reshape(B, (/size(B,1), size(B,3), size(B,2), size(B,4)/), order = (/1, 3, 2, 4/))
      B = reshape(B, (/size(B,1), size(B,2), size(B,4), size(B,3)/), order = (/1, 2, 4, 3/))
      do nw0 = 2, maxw1
      do nw2 = 2, maxw2
      do nw3 = 2, maxw3
        nw1 = 1
        N0 = (Al0 - mod(Al0, nw0)) / nw0
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist3d(B, N0, N2, N3, nw0, nw2, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
      enddo
      enddo
      enddo
    elseif (axis == 1) then
      B = cumsum4d(cumsum4d(cumsum4d(d, 2),3),4)
      B = reshape(B, (/size(B,2), size(B,1), size(B,3), size(B,4)/), order = (/2, 1, 3, 4/))
      B = reshape(B, (/size(B,1), size(B,3), size(B,2), size(B,4)/), order = (/1, 3, 2, 4/))
      B = reshape(B, (/size(B,1), size(B,2), size(B,4), size(B,3)/), order = (/1, 2, 4, 3/))
      do nw1 = 2, maxw1
      do nw2 = 2, maxw2
      do nw3 = 2, maxw3
        nw0 = 1
        N1 = (Al1 - mod(Al1, nw1)) / nw1
        N2 = (Al2 - mod(Al2, nw2)) / nw2
        N3 = (Al3 - mod(Al3, nw3)) / nw3
        k = hist3d(B, N1, N2, N3, nw1, nw2, nw3); call stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3) 
      enddo
      enddo
      enddo
    end if
  end subroutine help3d

  subroutine stat(k, Cn, kaves, deltas, nw0, nw1, nw2, nw3)
    double precision, allocatable :: k(:,:,:,:)
    integer, intent(in) :: nw0, nw1, nw2, nw3
    real(c_double), pointer :: Cn(:,:,:,:)                    
    real(c_double), pointer :: kaves(:,:,:,:)                    
    real(c_double), pointer :: deltas(:,:,:,:)                    
    double precision kave, v
       kave = sum(k) / real(size(k)); v = sum((k - kave)**2) / real(size(k))
       print *, "Cn with ", nw0, nw1, nw2, nw3, ":", (2.0 * kave - v) / (real(nw0*nw1*nw2*nw3)**2)
       kaves(nw0, nw1, nw2, nw3) = kave
       deltas(nw0, nw1, nw2, nw3) = nw0*nw1*nw2*nw3
       Cn(nw0, nw1, nw2, nw3) = (2.0 * kave - v) / (real(nw0*nw1*nw2*nw3)**2)
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


  function cumsum4d(d, axis)
      integer, intent(in) :: axis
      double precision, intent(in) :: d(:,:,:,:)
      double precision :: cumsum4d(size(d,1), size(d,2), size(d,3), size(d,4))
      integer :: i, j, k, l, N1, N2, N3, N4
      cumsum4d = d
      N1 = size(d,1); N2 = size(d,2); N3 = size(d,3); N4 = size(d,4)


      if (axis == 4) then
         !!$omp parallel do
         !do i = 1, N1
         !do j = 1, N2
         !do k = 1, N3
         do l = 1, N4 - 1
            !cumsum4d(i, j, k, l + 1) = cumsum4d(i, j, k, l) + d(i, j, k, l + 1)
            cumsum4d(:, :, :, l + 1) = cumsum4d(:, :, :, l) + d(:, :, :, l + 1)
         end do
         !end do
         !end do
         !end do
      else if (axis == 3) then
         !!$omp parallel do
         !do i = 1, N1
         !do j = 1, N2
         !do l = 1, N4
         do k = 1, N3 - 1
            !cumsum4d(i, j, k + 1, l) = cumsum4d(i, j, k, l) + d(i, j, k + 1, l)
            cumsum4d(:, :, k + 1, :) = cumsum4d(:, :, k, :) + d(:, :, k + 1, :)
         end do
         !end do
         !end do
         !end do
      else if (axis == 2) then
         !!$omp parallel do
         !do i = 1, N1
         !do k = 1, N3
         !do l = 1, N4
         do j = 1, N2 - 1
            !cumsum4d(i, j + 1, k, l) = cumsum4d(i, j, k, l) + d(i, j + 1, k, l)
            cumsum4d(:, j + 1, :, :) = cumsum4d(:, j, :, :) + d(:, j + 1, :, :)
         end do
         !end do
         !end do
         !end do
      else if (axis == 1) then
         !!$omp parallel do
         !do j = 1, N2 
         !do k = 1, N3
         !do l = 1, N4
         do i = 1, N1 - 1
            !cumsum4d(i + 1, j, k, l) = cumsum4d(i, j, k, l) + d(i + 1, j, k, l)
            cumsum4d(i + 1, :, :, :) = cumsum4d(i, :, :, :) + d(i + 1, :, :, :)
         end do
         !end do
         !end do
         !end do
      end if
   end function cumsum4d       

   function hist4d(A, N0, N1, N2, N3, nw0, nw1, nw2, nw3)
       double precision, intent(in) :: A(:,:,:,:) 
       integer, intent(in) :: N0, N1, N2, N3, nw0, nw1, nw2, nw3
       integer :: i, j, h, l, ihead, jhead, hhead, lhead
       double precision :: hist4d(N0,N1,N2,N3)
       do i = 1, N0
       ihead = i*nw0 
       do j = 1, N1
       jhead = j*nw1 
       do h = 1, N2
       hhead = h*nw2 
       !$omp parallel do
       do l = 1, N3
       lhead = l*nw3 
          if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead)
          else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead - nw0, jhead, hhead, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead)
          else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead)
          else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead - nw0, jhead, hhead - nw2, lhead)
          else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead, hhead, lhead - nw3)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead, lhead - nw3)
          else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw2, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3)
          else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead - nw2, lhead) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead) - A(ihead - nw0, jhead - nw1, hhead - nw2, lhead)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                           - A(ihead - nw0, jhead - nw1, hhead, lhead - nw3)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                           - A(ihead - nw0, jhead, hhead - nw2, lhead - nw3)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                           - A(ihead, jhead - nw1, hhead - nw2, lhead - nw3)
          else
             hist4d(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
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
   end function hist4d

   function hist4di(A, N0, N1, N2, N3, nw0, nw1, nw2, nw3)
       integer, intent(in) :: A(:,:,:,:) 
       integer, intent(in) :: N0, N1, N2, N3, nw0, nw1, nw2, nw3
       integer :: i, j, h, l, ihead, jhead, hhead, lhead
       integer :: hist4di(N0,N1,N2,N3)
       do i = 1, N0
       ihead = i*nw0 
       do j = 1, N1
       jhead = j*nw1 
       do h = 1, N2
       hhead = h*nw2 
       !$omp parallel do
       do l = 1, N3
       lhead = l*nw3 
          if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead)
          else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead - nw0, jhead, hhead, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead)
          else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead)
          else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead - nw0, jhead, hhead - nw2, lhead)
          else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead, hhead, lhead - nw3)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead)
          else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead, lhead - nw3)
          else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw2, lhead) - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3)
          else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           - A(ihead, jhead, hhead - nw2, lhead) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead - nw2, lhead) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead) - A(ihead - nw0, jhead - nw1, hhead - nw2, lhead)
          else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead - nw1, hhead, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead - nw1, hhead, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                           - A(ihead - nw0, jhead - nw1, hhead, lhead - nw3)
          else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead - nw0, jhead, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead - nw0, jhead, hhead - nw2, lhead) + A(ihead - nw0, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                           - A(ihead - nw0, jhead, hhead - nw2, lhead - nw3)
          else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
                           - A(ihead, jhead - nw1, hhead, lhead) - A(ihead, jhead, hhead - nw2, lhead) &
                           - A(ihead, jhead, hhead, lhead - nw3) &
                           + A(ihead, jhead - nw1, hhead - nw2, lhead) + A(ihead, jhead - nw1, hhead, lhead - nw3) &
                           + A(ihead, jhead, hhead - nw2, lhead - nw3) &
                           - A(ihead, jhead - nw1, hhead - nw2, lhead - nw3)
          else
             hist4di(i, j, h, l) = A(ihead, jhead, hhead, lhead) &
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

   function hist2d(B, N0, N1, nw0, nw1)
       double precision, intent(in) :: B(:,:,:,:) 
       integer, intent(in) :: N0, N1, nw0, nw1
       integer :: i, ihead, j, jhead
       double precision :: hist2d(N0, N1, size(B,3), size(B,4))
       do i = 1, N0
          ihead = i*nw0 
          !$omp parallel do
          do j = 1, N1
             jhead = j*nw1 
             if ( i == 1 .and. j == 1) then
                 hist2d(i, j, :, :) = B(ihead, jhead,:, :)
             else if ( i /= 1 .and. j == 1) then
                 hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw0, jhead, :, :)
             else if ( i == 1 .and. j /= 1) then
                 hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead, jhead - nw1, :, :)
             else 
                 hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw0, jhead, :, :) - B(ihead, jhead - nw1, :, :) &
                                      + B(ihead - nw0, jhead - nw1, :, :)
             end if
          end do
       end do
   end function hist2d

   function hist3d(B, N0, N1, N2, nw0, nw1, nw2)
       double precision, intent(in) :: B(:,:,:,:) 
       integer, intent(in) :: N0, N1, N2, nw0, nw1, nw2
       integer :: i, ihead, j, jhead, h, hhead
       double precision :: hist3d(N0, N1, N2, size(B,4))
       do i = 1, N0
          ihead = i*nw0 
          do j = 1, N1
             jhead = j*nw1
             do h = 1, N2
                hhead = h*nw2
                if ( i == 1 .and. j == 1 .and. h == 1) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :)
                else if ( j == 1 .and. i /= 1 .and. h == 1 ) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw0, jhead, hhead, :)
                else if ( i == 1 .and. j /= 1 .and. h == 1 ) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw1, hhead, :)
                else if ( i == 1 .and. h /= 1 .and. j == 1 ) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead, hhead - nw2, :)
                else if ( i /= 1 .and. j /= 1 .and. h == 1 ) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw0, jhead, hhead, :) &
                                        - B(ihead, jhead - nw1, hhead, :) + B(ihead - nw0, jhead - nw1, hhead, :)
                else if ( i /= 1 .and. j == 1 .and. h /= 1 ) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw0, jhead, hhead, :) &
                                       - B(ihead, jhead, hhead - nw2, :) + B(ihead - nw0, jhead, hhead - nw2, :)
                else if ( i == 1 .and. j /= 1 .and. h /= 1 ) then
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw1, hhead, :) &
                                       - B(ihead, jhead, hhead - nw2, :) + B(ihead, jhead - nw1, hhead - nw2, :)
                else
                   hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw0, jhead, hhead, :) &
                                       - B(ihead, jhead - nw1, hhead, :) - B(ihead, jhead, hhead - nw2, :) &
                                       + B(ihead, jhead - nw1, hhead - nw2, :) + B(ihead - nw0, jhead, hhead - nw2, :) &
                                       + B(ihead - nw0, jhead - nw1, hhead, :) - B(ihead - nw0, jhead - nw1, hhead - nw2, :)
                end if
             end do
          end do
       end do
   end function hist3d

end module costfort4d
