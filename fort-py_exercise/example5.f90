module example5
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len0
    integer(c_int) :: len1
    type(c_ptr) :: arr
  end type result

contains

  !function extract_plus(arr_length, array) bind(C, name="extract_plus")
  !  !DEC$ ATTRIBUTES DLLEXPORT :: extruct_plus
  !  integer(c_int), intent(in) :: arr_length
  !!!  integer(c_int), intent(in) :: array(arr_length)
!  !!  type(result) :: extract_plus
!    integer(c_int) :: plus_count
  !  integer(c_int), pointer :: work_array(:)
!
!    plus_count = count(array(:)>=0)
!    allocate(work_array(plus_count))
!    work_array(:) = pack(array(:), mask=(array(:)>=0))
!!
!    extract_plus%len = plus_count
!    extract_plus%arr = C_loc(work_array)
!  end function extract_plus


  function plus(arr_length0, arr_length1,  array) bind(C, name="plus1")
    !DEC$ ATTRIBUTES DLLEXPORT :: plus
    integer(c_int), intent(in) :: arr_length0
    integer(c_int), intent(in) :: arr_length1
    !integer(c_int), intent(in) :: array(arr_length)
    real(c_double), intent(in) :: array(arr_length1, arr_length0) ! u should understand the inverted order of the array sizes
    type(result) :: plus                                          ! fortran store the argv array with the inverted size order, i.e.,
    !integer(c_int), pointer :: work_array(:)                     ! python array A[dim1, dim2, dim3] -> fortran "array" should be transpose(A)[dim3, dim2, dim1] 
    real(c_double), pointer :: work_array(:,:) 

    allocate(work_array(3,2))
    work_array(1:3,1:2)=array(1:3,1:2)                           
    !plus%len0 =  arr_length0
    !plus%len1 =  arr_length1
    plus%len0 =  2
    plus%len1 =  3
    plus%arr = C_loc(work_array)
  end function plus




  subroutine delete_array(arr_length0, arr_length1, darray) bind(C, name="delete_array")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array
    integer(c_int), intent(in) :: arr_length0
    integer(c_int), intent(in) :: arr_length1
    type(c_ptr), value :: darray
    real(c_double), pointer :: work_array(:,:)

    call C_F_pointer(darray, work_array, [arr_length1, arr_length0])
    print *, "work_array:",work_array
    deallocate(work_array)
    darray = C_NULL_PTR
    print *, "Is work_array freed ?, work_array: ",work_array
  end subroutine delete_array

end module example5
