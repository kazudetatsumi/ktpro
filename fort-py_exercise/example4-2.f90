module example4
  use ISO_C_binding
  implicit none

contains

  function sum_all_func(arr_length, array) bind(C, name="sum_all_func")
    !DEC$ ATTRIBUTES DLLEXPORT :: sum_all_func
    integer(c_int), intent(in) :: arr_length
    integer(c_int), intent(in) :: array(arr_length)
    integer(c_int) :: sum_all_func

    sum_all_func = sum(array(:))

  end function sum_all_func


end module example4
