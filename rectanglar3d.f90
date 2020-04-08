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

  function rectanglar(datasize3, datasize2, datasize1, datasize0,  & 
      condition) bind(C, name="rectanglar")
    integer(c_int), intent(in) :: datasize0, datasize1, datasize2, datasize3
    integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    integer lb1, lb2, lb3, lb4, ub1, ub2, ub3, ub4, v, maxv, maxlb4, dq, dqnew
    integer local_maxv(datasize3)
    integer argmaxv_lb(4)
    integer argmaxv_ub(4)
    integer lb(4), ub(4), max_lb(4), min_ub(4), one(4)
    integer argmaxv_lub_ei(4, 2), argmaxv_lub_ee(4, 2)
    integer ei, ee, i
    type(result) :: rectanglar                              
    datasize = (/datasize0, datasize1, datasize2, datasize3/)
    print *, "datasize", datasize
    ei = 5
    ee = 35
    one = 1
    argmaxv_lub_ei =  argmaxv_lub(ei, ei, ee, datasize, condition)
    argmaxv_lub_ee =  argmaxv_lub(ee, ei, ee, datasize, condition)
    do i = 1, 4
       argmaxv_lb(i) = max(argmaxv_lub_ei(i, 1), argmaxv_lub_ee(i, 1))
       argmaxv_ub(i) = min(argmaxv_lub_ei(i, 2), argmaxv_lub_ee(i, 2))
    end do
    print *, "argmaxv_lb:", argmaxv_lb
    print *, "argmaxv_ub:", argmaxv_ub


    !call bondaries(ei, ee, max_lb, (/1, 1, 1, 1/), datasize, 1, datasize, condition)
    !call bondaries(ei, ee, min_ub, datasize, (/1, 1, 1, 1/), -1, datasize, condition)
    !do dq = 5, 1, -1
    !  call argmaxvlight(ee, dq, max_lb, min_ub, argmaxv_lb, argmaxv_ub, datasize, condition)
    !  max_lb = argmaxv_lb 
    !  min_ub = argmaxv_ub 
    !  do i = 2, 4
    !     if ( max_lb(i) - dq > 0 ) then 
    !         max_lb(i) = max_lb(i) - dq
    !     else
    !         max_lb(i) = 1
    !     end if
    !     if ( min_ub(i) + dq <= datasize(i) ) then 
    !         min_ub(i) = min_ub(i) + dq
    !     else
    !         min_ub(i) = datasize(i)
    !     end if
    !  end do
    !end do
    !print *, "argmaxv_lb:", argmaxv_lb
    !print *, "argmaxv_ub:", argmaxv_ub

    !call bondaries(ei, ee, max_lb, (/1, 1, 1, 1/), datasize, 1, datasize, condition)
    !call bondaries(ei, ee, min_ub, datasize, (/1, 1, 1, 1/), -1, datasize, condition)
    !do dq = 5, 1, -1
    !  call argmaxvlight(ei, dq, max_lb, min_ub, argmaxv_lb, argmaxv_ub, datasize, condition)
    !  max_lb = argmaxv_lb 
    !  min_ub = argmaxv_ub 
    !  do i = 2, 4
    !     if ( max_lb(i) - dq > 0 ) then 
    !         max_lb(i) = max_lb(i) - dq
    !     else
    !         max_lb(i) = 1
    !     end if
    !     if ( min_ub(i) + dq <= datasize(i) ) then 
    !         min_ub(i) = min_ub(i) + dq
    !     else
    !         min_ub(i) = datasize(i)
    !     end if
    !  end do
    !end do
    !print *, "argmaxv_lb:", argmaxv_lb
    !print *, "argmaxv_ub:", argmaxv_ub

       

    rectanglar%lb_0 = argmaxv_lb(1)
    rectanglar%lb_1 = argmaxv_lb(2)
    rectanglar%lb_2 = argmaxv_lb(3)
    rectanglar%ub_0 = argmaxv_ub(1)
    rectanglar%ub_1 = argmaxv_ub(2)
    rectanglar%ub_2 = argmaxv_ub(3)

  end function rectanglar


  subroutine bondaries(ei, ee, qb, iniq, endq, dq,  datasize, condition)
    integer, intent(out) :: qb(4)
    integer, intent(in)  :: iniq(4), endq(4), dq
    integer, intent(in)  :: datasize(4), condition(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in)  :: ei, ee
    integer  qid
    qb = iniq
    do qid = iniq(2), endq(2), dq
      if ( sum(condition(ei:ee, qid, :, :)) > 0.1 ) then
           qb(2) = qid
           exit
      end if
    end do
    do qid = iniq(3), endq(3), dq
      if ( sum(condition(ei:ee, :, qid, :)) > 0.1 ) then
           qb(3) = qid
           exit
      end if
    end do
    do qid = iniq(4), endq(4), dq
      if ( sum(condition(ei:ee, :, :, qid)) > 0.1 ) then
           qb(4) = qid
           exit
      end if
    end do
  end subroutine bondaries


  subroutine argmaxv(ei, ee, max_lb, min_ub, argmaxv_lb, argmaxv_ub, datasize, condition)
    integer, intent(in)  :: ei, ee
    integer, intent(in)  :: max_lb(4), min_ub(4)
    integer, intent(in)  :: datasize(4), condition(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(out) :: argmaxv_lb(4), argmaxv_ub(4)
    integer  lb1, lb2, lb3, lb4, ub1, ub2, ub3, ub4, v, maxlb4(3)
    integer  local_maxv(datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_lb(4, datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_ub(4, datasize(2), datasize(3), datasize(4))
    local_maxv = 0
    argmaxv_lb = 0
    argmaxv_ub = 0


    lb1 = ei
    ub1 = ee
    do lb2 = max_lb(2), min_ub(2) - 15, 1
    do ub2 = lb2 + 15, min_ub(2), 1
    !$omp parallel do
    do lb3 = max_lb(3), min_ub(3) - 35, 1
    do ub3 = lb3 + 35, min_ub(3), 1
    do lb4 = max_lb(4), min_ub(4) - 35, 1
    do ub4 = lb4 + 35, min_ub(4), 1
    if ( real(sum(condition(lb1:ub1, lb2, lb3:ub3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub3-lb3)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, ub2, lb3:ub3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub3-lb3)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, lb3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, ub3, lb4:ub4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub4-lb4) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, lb3:ub3, lb4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub3-lb3) ) then  
    if ( real(sum(condition(lb1:ub1, lb2:ub2, lb3:ub3, ub4))) >= 0.95 * (ub1-lb1)*(ub2-lb2)*(ub3-lb3) ) then

    v = (ub1 - lb1)*(ub2 - lb2)*(ub3 - lb3)*(ub4 - lb4)
    if (real(sum(condition(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4))) >= 0.9 * v .and. v > local_maxv(lb2, lb3, lb4) ) then
    local_maxv(lb2, lb3, lb4) = v
    local_argmaxv_lb(:, lb2, lb3, lb4) = (/lb1, lb2, lb3, lb4/)
    local_argmaxv_ub(:, lb2, lb3, lb4) = (/ub1, ub2, ub3, ub4/)
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
    maxlb4 = maxloc(local_maxv)
    print *, maxlb4, local_maxv(maxlb4(1), maxlb4(2), maxlb4(3))
    print *, local_argmaxv_lb(:,maxlb4(1), maxlb4(2), maxlb4(3))
    print *, local_argmaxv_ub(:,maxlb4(1), maxlb4(2), maxlb4(3))

  end subroutine argmaxv


  subroutine argmaxvlight(eidx, dq, max_lb, min_ub, argmaxv_lb, argmaxv_ub, datasize, condition)
    integer, intent(in)  :: eidx, dq
    integer, intent(in)  :: max_lb(4), min_ub(4)
    integer, intent(in)  :: datasize(4), condition(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(inout) :: argmaxv_lb(4), argmaxv_ub(4)
    integer  lb2, lb3, lb4, ub2, ub3, ub4, v, maxlb4(3)
    integer  local_maxv(datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_lb(4, datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_ub(4, datasize(2), datasize(3), datasize(4))
    local_maxv = 0


    do lb2 = max_lb(2), min_ub(2) - 10, dq
    do ub2 = lb2 + 10, min_ub(2), dq
    !$omp parallel do
    do lb3 = max_lb(3), min_ub(3) - 30, dq
    do ub3 = lb3 + 30, min_ub(3), dq
    do lb4 = max_lb(4), min_ub(4) - 30, dq
    do ub4 = lb4 + 30, min_ub(4), dq
    if ( real(sum(condition(eidx, lb2, lb3:ub3, lb4:ub4))) >= 0.95 * (ub3-lb3)*(ub4-lb4) ) then  
    if ( real(sum(condition(eidx, ub2, lb3:ub3, lb4:ub4))) >= 0.95 * (ub3-lb3)*(ub4-lb4) ) then  
    if ( real(sum(condition(eidx, lb2:ub2, lb3, lb4:ub4))) >= 0.95 * (ub2-lb2)*(ub4-lb4) ) then  
    if ( real(sum(condition(eidx, lb2:ub2, ub3, lb4:ub4))) >= 0.95 * (ub2-lb2)*(ub4-lb4) ) then  
    if ( real(sum(condition(eidx, lb2:ub2, lb3:ub3, lb4))) >= 0.95 * (ub2-lb2)*(ub3-lb3) ) then  
    if ( real(sum(condition(eidx, lb2:ub2, lb3:ub3, ub4))) >= 0.95 * (ub2-lb2)*(ub3-lb3) ) then

    v = (ub2 - lb2)*(ub3 - lb3)*(ub4 - lb4)
    if (real(sum(condition(eidx, lb2:ub2, lb3:ub3, lb4:ub4))) >= 0.9 * v .and. v > local_maxv(lb2, lb3, lb4) ) then
    local_maxv(lb2, lb3, lb4) = v
    local_argmaxv_lb(:, lb2, lb3, lb4) = (/eidx, lb2, lb3, lb4/)
    local_argmaxv_ub(:, lb2, lb3, lb4) = (/eidx, ub2, ub3, ub4/)
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
    if ( maxval(local_maxv) > 0 ) then
      maxlb4 = maxloc(local_maxv)
      argmaxv_lb = local_argmaxv_lb(:, maxlb4(1), maxlb4(2), maxlb4(3))
      argmaxv_ub = local_argmaxv_ub(:, maxlb4(1), maxlb4(2), maxlb4(3))
    end if
    print *, argmaxv_lb
    print *, argmaxv_ub

  end subroutine argmaxvlight

  function argmaxv_lub(eidx, ei, ee, datasize, condition)
    integer, intent(in)    :: eidx, ei, ee
    integer, intent(in)    :: datasize(4), condition(datasize(1), datasize(2), datasize(3), datasize(4))
    integer :: max_lb(4), min_ub(4)
    integer :: dq, i
    integer :: argmaxv_lub(4, 2)

    call bondaries(ei, ee, max_lb, (/1, 1, 1, 1/), datasize, 1, datasize, condition)
    call bondaries(ei, ee, min_ub, datasize, (/1, 1, 1, 1/), -1, datasize, condition)
    do dq = 5, 1, -1
      call argmaxvlight(eidx, dq, max_lb, min_ub, argmaxv_lub(:, 1), argmaxv_lub(:, 2), datasize, condition)
      max_lb = argmaxv_lub(:, 1) 
      min_ub = argmaxv_lub(:, 2) 
      do i = 2, 4
         if ( max_lb(i) - dq > 0 ) then 
             max_lb(i) = max_lb(i) - dq
         else
             max_lb(i) = 1
         end if
         if ( min_ub(i) + dq <= datasize(i) ) then 
             min_ub(i) = min_ub(i) + dq
         else
             min_ub(i) = datasize(i)
         end if
      end do
    end do
  end function argmaxv_lub

end module costfort4d
