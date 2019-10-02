module example2
  use ISO_C_binding
  implicit none

  ! scalar
  integer(C_int) :: int_val 

  ! array
  integer(C_int), parameter :: array_length = 10
  integer(C_int) :: array(array_length)

contains


  function add(val_a, val_b) bind(C, name="add")
    !DEC$ ATTRIBUTES DLLEXPORT :: add
    integer(c_int), intent(in) :: val_a, val_b
    integer(c_int) :: add

    add = val_a + val_b

  end function add
 
end module example2
