module costfort1din4d
  use ISO_C_binding
  implicit none
  type, bind(C) :: result1d
    integer(c_int) :: len3
    type(c_ptr) :: arr
    type(c_ptr) :: kavearr
    type(c_ptr) :: darr
  end type result1d

contains

  function cost1d(maxw, Bl3, Bl2, Bl1, Bl0, B, CDB) bind(C, name="cost1d")
    !DEC$ ATTRIBUTES DLLEXPORT :: cost1d
    integer(c_int), intent(in) :: maxw
    integer(c_int), intent(in) :: Bl0
    integer(c_int), intent(in) :: Bl1
    integer(c_int), intent(in) :: Bl2
    integer(c_int), intent(in) :: Bl3
    real(c_double), intent(in) :: B(Bl0, Bl1, Bl2, Bl3)                         
    integer(c_int), intent(in) :: CDB(Bl0, Bl1, Bl2, Bl3)                         
    type(result1d) :: cost1d                                  
    real(c_double), pointer :: Cn(:)                    
    real(c_double), pointer :: kaves(:)                    
    real(c_double), pointer :: deltas(:)                    
    double precision, allocatable :: k(:,:,:,:)
    double precision, allocatable :: kcond(:,:,:,:)
    double precision, allocatable :: knonzero(:)
    integer nw3
    integer N3
    integer i, ihead
    real kave,v
    allocate(Cn(maxw))
    allocate(kaves(maxw))
    allocate(deltas(maxw))
    Cn(:) = 0.0
    kaves(:) = 0.0
    deltas(:) = 0.0
    do nw3 = 1, maxw
    N3 = (Bl3 - mod(Bl3, nw3)) / nw3
    allocate(k(Bl0, Bl1, Bl2, N3))
    allocate(kcond(Bl0, Bl1, Bl2, N3))
       k = hist1d(B, N3, nw3)
       kcond = hist1di(CDB, N3, nw3)
       knonzero = pack(k, kcond == maxval(kcond))
       kave = sum(knonzero) / real(size(knonzero))
       v = sum((knonzero - kave)**2) / real(size(knonzero))
       Cn(nw3) = (2.0 * kave - v) / (real(nw3)**2)
       deltas(nw3) = real(nw3)
       kaves(nw3) = kave
       deallocate(k,kcond,knonzero)
    end do

    cost1d%len3 =  maxw
    cost1d%arr = C_loc(Cn(1:maxw))
    cost1d%kavearr = C_loc(kaves(1:maxw))
    cost1d%darr = C_loc(deltas(1:maxw))
  end function cost1d



  function hist1d(B, N, nw)
     double precision, intent(in) :: B(:,:,:,:)
     integer, intent(in) :: N, nw
     integer :: i, ihead
     double precision :: hist1d(size(B,1), size(B,2), size(B,3), N)
     hist1d(:,:,:,1) = B(:,:,:,nw)
     do i = 2, N
        ihead = i*nw
        hist1d(:,:,:,i) = B(:,:,:,ihead) - B(:,:,:,ihead - nw)
     end do
  end function hist1d

  function hist1di(CDB, N, nw)
     integer, intent(in) :: CDB(:,:,:,:)
     integer, intent(in) :: N, nw
     integer :: i, ihead
     integer :: hist1di(size(CDB,1), size(CDB,2), size(CDB,3), N)
     hist1di(:,:,:,1) = CDB(:,:,:,nw)
     do i = 2, N
        ihead = i*nw
        hist1di(:,:,:,i) = CDB(:,:,:,ihead) - CDB(:,:,:,ihead - nw)
     end do
  end function hist1di


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

end module costfort1din4d
