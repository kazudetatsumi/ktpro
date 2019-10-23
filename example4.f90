module example4
  use ISO_C_binding
  implicit none

contains

  subroutine sum_all_sub(arr_length, array, result) bind(C, name="sum_all_sub")
    !DEC$ ATTRIBUTES DLLEXPORT :: sum_all_sub
    integer(c_int), intent(in) :: arr_length
    integer(c_int), intent(in) :: array(arr_length)
    integer(c_int), intent(out) :: result
    
    result = sum(array(:)) + 10000
  end subroutine sum_all_sub


  function sum_all_func(arr_length, array) bind(C, name="sum_all_func")
    !DEC$ ATTRIBUTES DLLEXPORT :: sum_all_func
    integer(c_int), intent(in) :: arr_length
    integer(c_int), intent(in) :: array(arr_length)
    integer(c_int) :: sum_all_func

    sum_all_func = sum(array(:)) 

  end function sum_all_func

end module example4

