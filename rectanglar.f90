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

contains

  function rectanglar(datasize3, datasize2, datasize1, datasize0,& 
  !function rectanglar(datasize2, datasize1, datasize0,& 
      condition) bind(C, name="rectanglar")
    integer(c_int), intent(in) :: datasize0, datasize1, datasize2, datasize3
    !integer(c_int), intent(in) :: datasize0, datasize1, datasize2
    integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    !integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2)                         
    !integer lb0, lb1, lb2, ub0, ub1, ub2, v, maxv, max_lb0, max_lb1, max_lb2, min_ub0, min_ub1, min_ub2
    integer lb0, lb1, lb2, lb3, ub0, ub1, ub2, ub3, v, maxv, max_lb0, max_lb1, max_lb2, max_lb3,  min_ub0, min_ub1, min_ub2, min_ub3
    integer local_maxv(datasize3)
    !integer argmaxv_lb0, argmaxv_lb1, argmaxv_lb2, argmaxv_ub0, argmaxv_ub1, argmaxv_ub2
    !integer argmaxv_lb0, argmaxv_lb1, argmaxv_lb2, argmaxv_lb3, argmaxv_ub0, argmaxv_ub1, argmaxv_ub2, argmaxv_ub3
    integer argmaxv_lb0(datasize3), argmaxv_lb1(datasize3), argmaxv_lb2(datasize3), argmaxv_lb3(datasize3), argmaxv_ub0(datasize3)
    integer argmaxv_ub1(datasize3), argmaxv_ub2(datasize3), argmaxv_ub3(datasize3)
    integer maxlb3
    type(result) :: rectanglar                              

    print *, datasize0, datasize1, datasize2, datasize3
    print *, shape(condition)
    print *, "maxval of condition  in fortran",maxval(condition)
    print *, "sum of condition in fortran", sum(condition(:,:,:,:))
    !print *, "sum of condition in fortran", sum(condition(:,:,:))
    print *, "condition(34, 20, 20, :):", condition(34, 20, 20, :)
    !print *, "condition(20, 20, :):", condition(20, 20, :)


    max_lb0 = 1
    max_lb1 = 1
    max_lb2 = 1
    max_lb3 = 1
    min_ub0 = datasize0
    min_ub1 = datasize1
    min_ub2 = datasize2
    min_ub3 = datasize3
    !maxv = 1
    local_maxv = 0
    argmaxv_lb0 = 0
    argmaxv_lb1 = 0
    argmaxv_lb2 = 0
    argmaxv_lb3 = 0
    argmaxv_ub0 = 0
    argmaxv_ub1 = 0
    argmaxv_ub2 = 0
    argmaxv_ub3 = 0


    do lb0 = 1, datasize0
      !if ( sum(condition(lb0, :, :)) > 0.1 ) then
      if ( sum(condition(lb0, :, :, :)) > 0.1 ) then
           max_lb0 = lb0
           exit
      end if
    end do
    print *, "max_lb0", max_lb0
    do lb1 = 1, datasize1
      !if ( sum(condition(:, lb1, :)) > 0.1 ) then
      if ( sum(condition(:, lb1, :, :)) > 0.1 ) then
           max_lb1 = lb1
           exit
      end if
    end do
    print *, "max_lb1", max_lb1
    do lb2 = 1, datasize2
      !if ( sum(condition(:, :, lb2)) > 0.1 ) then
      if ( sum(condition(:, :, lb2, :)) > 0.1 ) then
           max_lb2 = lb2
           exit
      end if
    end do
    print *, "max_lb2", max_lb2
    do lb3 = 1, datasize3
      if ( sum(condition(:, :, :, lb3)) > 0.1 ) then
           max_lb3 = lb3
           exit
      end if
    end do
    print *, "max_lb3", max_lb3


    do ub0 = datasize0, 1, -1
      !if ( sum(condition(ub0, :, :)) > 0.1 ) then
      if ( sum(condition(ub0, :, :, :)) > 0.1 ) then
           min_ub0 = ub0
           exit
      end if
    end do
    print *, "min_ub0", min_ub0
    do ub1 = datasize1, 1, -1
      !if ( sum(condition(:, ub1, :)) > 0.1 ) then
      if ( sum(condition(:, ub1, :, :)) > 0.1 ) then
           min_ub1 = ub1
           exit
      end if
    end do
    print *, "min_ub1", min_ub1
    do ub2 = datasize2, 1, -1
      !if ( sum(condition(:, :, ub2)) > 0.1 ) then
      if ( sum(condition(:, :, ub2, :)) > 0.1 ) then
           min_ub2 = ub2 
           exit
      end if
    end do
    print *, "min_ub2", min_ub2
    do ub3 = datasize3, 1, -1
      if ( sum(condition(:, :, :, ub3)) > 0.1 ) then
           min_ub3 = ub3 
           exit
      end if
    end do
    print *, "min_ub3", min_ub3
       
       

       

    do lb0 = max_lb0, min_ub0, 5
    do ub0 = lb0, min_ub0, 5
    do lb1 = max_lb1, min_ub1, 5
    do ub1 = lb1, min_ub1, 5
    do lb2 = max_lb2, min_ub2, 5
    do ub2 = lb2, min_ub2, 5
    !$omp parallel do
    do lb3 = max_lb3, min_ub3, 5
    do ub3 = lb3, min_ub3, 5

      !v = (ub0 - lb0)*(ub1 - lb1)*(ub2 - lb2)
      v = (ub0 - lb0)*(ub1 - lb1)*(ub2 - lb2)*(ub3 - lb3)
      !if (real(sum(condition(35, lb0:ub0, lb1:ub1, lb2:ub2))) >= 0.9 * v .and. v > maxv ) then
      !if (real(sum(condition(lb0:ub0, lb1:ub1, lb2:ub2))) >= 0.9 * v .and. v > maxv ) then
      !if (real(sum(condition(lb0:ub0, lb1:ub1, lb2:ub2, lb3:ub3))) >= 0.9 * v .and. v > maxv ) then
      if (real(sum(condition(lb0:ub0, lb1:ub1, lb2:ub2, lb3:ub3))) >= 0.9 * v .and. v > local_maxv(lb3) ) then
          !maxv = v
          local_maxv(lb3) = v
          argmaxv_lb0(lb3) = lb0
          argmaxv_ub0(lb3) = ub0
          argmaxv_lb1(lb3) = lb1
          argmaxv_ub1(lb3) = ub1
          argmaxv_lb2(lb3) = lb2
          argmaxv_ub2(lb3) = ub2
          argmaxv_lb3(lb3) = lb3
          argmaxv_ub3(lb3) = ub3
          !print *, maxv, argmaxv_lb0, argmaxv_ub0, argmaxv_lb1, argmaxv_ub1, argmaxv_lb2, argmaxv_ub2, argmaxv_lb3, argmaxv_ub3
      end if
    end do
    end do
    end do
    end do

    end do
    end do
    end do
    end do
    !print *, "maxv=", maxv, " with bondaries: ", argmaxv_lb0, argmaxv_ub0, argmaxv_lb1, argmaxv_ub1, argmaxv_lb2, argmaxv_ub2
    !!do lb3=1,datasize3
    !      print *, lb3, local_maxv(lb3), argmaxv_lb0(lb3), argmaxv_ub0(lb3), argmaxv_lb1(lb3), argmaxv_ub1(lb3), argmaxv_lb2(lb3), &
    !                argmaxv_ub2(lb3), argmaxv_lb3(lb3), argmaxv_ub3(lb3)
    !end do
    !alb3 = maxloc(local_maxv)
    !print *, alb3, local_maxv(alb3), argmaxv_lb0(alb3), argmaxv_ub0(alb3), argmaxv_lb1(alb3), argmaxv_ub1(alb3), argmaxv_lb2(alb3), &
    !         argmaxv_ub2(alb3), argmaxv_lb3(alb3), argmaxv_ub3(alb3)
    maxlb3 = maxloc(local_maxv, 1)
    print *, maxlb3, local_maxv(maxlb3), argmaxv_lb0(maxlb3), argmaxv_ub0(maxlb3), argmaxv_lb1(maxlb3), argmaxv_ub1(maxlb3),  &
             argmaxv_lb2(maxlb3), argmaxv_ub2(maxlb3), argmaxv_lb3(maxlb3), argmaxv_ub3(maxlb3)

    rectanglar%lb_0 = argmaxv_lb0(1)
    rectanglar%lb_1 = argmaxv_lb1(1)
    rectanglar%lb_2 = argmaxv_lb2(1)
    rectanglar%ub_0 = argmaxv_ub0(1)
    rectanglar%ub_1 = argmaxv_ub1(1)
    rectanglar%ub_2 = argmaxv_ub2(1)

  end function rectanglar

end module costfort4d
