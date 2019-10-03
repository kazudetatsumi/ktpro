module example5
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: len
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


  function plus(arr_length, array) bind(C, name="plus1")
    !DEC$ ATTRIBUTES DLLEXPORT :: plus
    integer(c_int), intent(in) :: arr_length
    !integer(c_int), intent(in) :: array(arr_length)
    real(c_double), intent(in) :: array(arr_length)
    type(result) :: plus
    !integer(c_int), pointer :: work_array(:)
    real(c_double), pointer :: work_array(:)

    allocate(work_array(arr_length))
    work_array(:) = array(:) + 1.0
    plus%len = arr_length
    plus%arr = C_loc(work_array)
  end function plus




  subroutine delete_array(arr_length, array) bind(C, name="delete_array")
    !DEC$ ATTRIBUTES DLLEXPORT :: delete_array
    integer(c_int), intent(in) :: arr_length
    type(c_ptr), value :: array
    real(c_double), pointer :: work_array(:)

    call C_F_pointer(array, work_array, [arr_length])
    deallocate(work_array)
    array = C_NULL_PTR
  end subroutine delete_array

end module example5
