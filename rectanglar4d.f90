!Kazuyoshi TATSUMI 2020/04/06

module costfort4d
  use ISO_C_binding
  implicit none
  type, bind(C) :: result
    integer(c_int) :: lb_0
    integer(c_int) :: ub_0
    integer(c_int) :: lb_1
    integer(c_int) :: ub_1
    integer(c_int) :: lb_2
    integer(c_int) :: ub_2
  end type result
  integer datasize(4)

contains

  function rectanglar(datasize3, datasize2, datasize1, datasize0,& 
      condition) bind(C, name="rectanglar")
    integer(c_int), intent(in) :: datasize0, datasize1, datasize2, datasize3
    integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    integer lb1, lb2, lb3, lb4, ub1, ub2, ub3, ub4, v, maxv, maxlb4
    integer local_maxv(datasize3)
    integer argmaxv_lb(4,datasize3)
    integer argmaxv_ub(4,datasize3)
    integer lb(4), ub(4), max_lb(4), min_ub(4)
    type(result) :: rectanglar                              
    datasize = (/datasize0, datasize1, datasize2, datasize3/)


    max_lb = 1
    min_ub = datasize
    local_maxv = 0
    argmaxv_lb = 0
    argmaxv_ub = 0


    max_lb(1) = 5
    do lb2 = 1, datasize(2)
      if ( sum(condition(5:35, lb2, :, :)) > 0.1 ) then
           max_lb(2) = lb2
           exit
      end if
    end do
    do lb3 = 1, datasize(3)
      if ( sum(condition(5:35, :, lb3, :)) > 0.1 ) then
           max_lb(3) = lb3
           exit
      end if
    end do
    do lb4 = 1, datasize(4)
      if ( sum(condition(5:35, :, :, lb4)) > 0.1 ) then
           max_lb(4) = lb4
           exit
      end if
    end do
    print *, "max_lb(:)", max_lb(:)

    min_ub(1) = 35
    do ub2 = datasize(2), 1, -1
      if ( sum(condition(5:35, ub2, :, :)) > 0.1 ) then
           min_ub(2) = ub2
           exit
      end if
    end do
    do ub3 = datasize(3), 1, -1
      if ( sum(condition(5:35, :, ub3, :)) > 0.1 ) then
           min_ub(3) = ub3 
           exit
      end if
    end do
    do ub4 = datasize(4), 1, -1
      if ( sum(condition(5:35, :, :, ub4)) > 0.1 ) then
           min_ub(4) = ub4
           exit
      end if
    end do
    print *, "min_ub(:)", min_ub(:)
       
       

    lb1 = 5
    ub1 = 35
    do lb2 = max_lb(2), min_ub(2) - 10, 5
    do ub2 = lb2 + 10, min_ub(2), 1
    do lb3 = max_lb(3), min_ub(3) - 10, 5
    do ub3 = lb3 + 10, min_ub(3), 1
    !$omp parallel do
    do lb4 = max_lb(4), min_ub(4) - 10, 5
    do ub4 = lb4 + 10, min_ub(4), 5
    if ( real(sum(condition(lb1:ub1, lb2, lb3:ub3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub3-lb3)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, ub2, lb3:ub3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub3-lb3)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, lb3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, ub3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, lb3:ub3, lb4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub3-lb3) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, lb3:ub3, ub4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub3-lb3) ) then

    v = (ub1 - lb1)*(ub2 - lb2)*(ub3 - lb3)*(ub4 - lb4)
    if (real(sum(condition(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4))) >= 0.9 * v .and. v > local_maxv(lb4) ) then
    local_maxv(lb4) = v
    argmaxv_lb(:, lb4) = (/lb1, lb2, lb3, lb4/)
    argmaxv_ub(:, lb4) = (/ub1, ub2, ub3, ub4/)
    end if
    end if
    end if
    end if
    end if
    end if
    end if
    end do
    end do

    end do
    end do
    end do
    end do
    maxlb4 = maxloc(local_maxv, 1)
    print *, maxlb4, local_maxv(maxlb4)
    print *, argmaxv_lb(:,maxlb4)
    print *, argmaxv_ub(:,maxlb4)

    rectanglar%lb_0 = argmaxv_lb(1, maxlb4)
    rectanglar%lb_1 = argmaxv_lb(2, maxlb4)
    rectanglar%lb_2 = argmaxv_lb(3, maxlb4)
    rectanglar%ub_0 = argmaxv_ub(1, maxlb4)
    rectanglar%ub_1 = argmaxv_ub(2, maxlb4)
    rectanglar%ub_2 = argmaxv_ub(3, maxlb4)

  end function rectanglar

end module costfort4d
