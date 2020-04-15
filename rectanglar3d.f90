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
    integer(c_int) :: lb_3
    integer(c_int) :: ub_3
  end type result
  integer datasize(4), ei, ee
  real a_v, a_b
  integer, allocatable ::  cond(:,:,:,:)

contains

  function rectanglar(uselight, e_i, e_e, av, ab, datasize3, datasize2, datasize1, datasize0, condition) bind(C, name="rectanglar")
    logical(c_bool), intent(in) :: uselight
    integer(c_int), intent(in) :: e_i, e_e, datasize0, datasize1, datasize2, datasize3
    integer(c_int), intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    real(c_double), intent(in) :: av, ab

    integer argmaxv_lb(4)
    integer argmaxv_ub(4)
    integer argmaxv_lub_ei(4, 2), argmaxv_lub_ee(4, 2)
    integer didx
    type(result) :: rectanglar                              
    datasize = (/datasize0, datasize1, datasize2, datasize3/)
    cond = condition
    print *, "datasize:", datasize
    ! the lower and upper energy boundaries
    ei = e_i + 1
    ee = e_e + 1
    ! parameters softening the conditions to select volume and boundaries
    ! a_v = 0.9, ab = 0.999
    a_v = av
    a_b = ab
    print *, "a_v, a_b", a_v, a_b
    print *, "uselight:", uselight
    if ( uselight ) then
       !!search upper and lower boundaries on three q axis for two partial matrices condition(ei, :, :, :) and condition(ee, :, :, :)
       !!and select the smaller (larger) upper (lower) boundaries.
       argmaxv_lub_ei =  argmaxv3d_lub(ei)
       argmaxv_lub_ee =  argmaxv3d_lub(ee)
       do didx = 2, 4
          argmaxv_lb(didx) = max(argmaxv_lub_ei(didx, 1), argmaxv_lub_ee(didx, 1))
          argmaxv_ub(didx) = min(argmaxv_lub_ei(didx, 2), argmaxv_lub_ee(didx, 2))
       end do
    else
       !!search upper and lower boundaries on three q axis for the condition(ei:ee, :, :, :)
       argmaxv_lub_ei =  argmaxv4d_lub()
       argmaxv_lb =  argmaxv_lub_ei(:, 1)
       argmaxv_ub =  argmaxv_lub_ei(:, 2)
    end if
    argmaxv_lb(1) = ei
    argmaxv_ub(1) = ee
    print *, "argmaxv_lb:", argmaxv_lb
    print *, "argmaxv_ub:", argmaxv_ub
       

    rectanglar%lb_0 = argmaxv_lb(4) - 1
    rectanglar%lb_1 = argmaxv_lb(3) - 1
    rectanglar%lb_2 = argmaxv_lb(2) - 1
    rectanglar%lb_3 = argmaxv_lb(1) - 1
    rectanglar%ub_0 = argmaxv_ub(4) - 1
    rectanglar%ub_1 = argmaxv_ub(3) - 1
    rectanglar%ub_2 = argmaxv_ub(2) - 1
    rectanglar%ub_3 = argmaxv_ub(1) - 1
 
  end function rectanglar


  function argmaxv3d_lub(eidx)
    integer, intent(in)    :: eidx
    integer :: max_lb(4), min_ub(4)
    integer :: dq, didx
    integer :: argmaxv3d_lub(4, 2)
    print *, "entering function argmaxv3d_lub"

    argmaxv3d_lub = 0 
    call boundaries(max_lb, (/1, 1, 1, 1/), datasize, 1)
    call boundaries(min_ub, datasize, (/1, 1, 1, 1/), -1)
    do dq = 8, 4, -2
      call argmaxvlight(eidx, dq, max_lb, min_ub, argmaxv3d_lub(:, 1), argmaxv3d_lub(:, 2))
      max_lb = argmaxv3d_lub(:, 1) 
      min_ub = argmaxv3d_lub(:, 2) 
      do didx = 2, 4
         if ( max_lb(didx) - int(dq/2) > 0 ) then 
             max_lb(didx) = max_lb(didx) - int(dq/2)
         else
             max_lb(didx) = 1
         end if
         if ( min_ub(didx) + int(dq/2) <= datasize(didx) ) then 
             min_ub(didx) = min_ub(didx) + int(dq/2)
         else
             min_ub(didx) = datasize(didx)
         end if
      end do
    end do
  end function argmaxv3d_lub

  
  function argmaxv4d_lub
    integer :: max_lb(4), min_ub(4)
    integer :: dq, didx
    integer :: argmaxv4d_lub(4, 2)

    print *, "entering function argmaxv4d_lub"

    argmaxv4d_lub = 0 
    call boundaries(max_lb, (/1, 1, 1, 1/), datasize, 1)
    call boundaries(min_ub, datasize, (/1, 1, 1, 1/), -1)
    do dq = 17, 1, -2
      call argmaxv(dq, max_lb, min_ub, argmaxv4d_lub(:, 1), argmaxv4d_lub(:, 2))
      max_lb = argmaxv4d_lub(:, 1) 
      min_ub = argmaxv4d_lub(:, 2) 
      do didx = 2, 4
         if ( max_lb(didx) - dq > 0 ) then 
             max_lb(didx) = max_lb(didx) - dq
         else
             max_lb(didx) = 1
         end if
         if ( min_ub(didx) + dq <= datasize(didx) ) then 
             min_ub(didx) = min_ub(didx) + dq
         else
             min_ub(didx) = datasize(didx)
         end if
      end do
    end do
  end function argmaxv4d_lub


  subroutine boundaries(qb, iniq, endq, dq)
    integer, intent(out) :: qb(4)
    integer, intent(in)  :: iniq(4), endq(4), dq
    integer  qid
    integer cond_order(4, 4)
    qb = iniq

    do qid = iniq(2), endq(2), dq
      if ( sum(cond(ei:ee, qid, :, :)) > 0.1 ) then
           qb(2) = qid
           exit
      end if
    end do
    do qid = iniq(3), endq(3), dq
      if ( sum(cond(ei:ee, :, qid, :)) > 0.1 ) then
           qb(3) = qid
           exit
      end if
    end do
    do qid = iniq(4), endq(4), dq
      if ( sum(cond(ei:ee, :, :, qid)) > 0.1 ) then
           qb(4) = qid
           exit
      end if
    end do
  end subroutine boundaries


  subroutine argmaxv(dq, max_lb, min_ub, argmaxv_lb, argmaxv_ub)
    integer, intent(in)  :: dq
    integer, intent(in)  :: max_lb(4), min_ub(4)
    integer, intent(inout) :: argmaxv_lb(4), argmaxv_ub(4)
    integer  lb2, lb3, lb4, ub1, ub2, ub3, ub4, v, maxlb4(3), diff(4)
    integer  local_maxv(datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_lb(4, datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_ub(4, datasize(2), datasize(3), datasize(4))
    local_maxv = 0

    if ( sum(argmaxv_ub - argmaxv_lb) /=0 ) then
        diff = int((argmaxv_ub - argmaxv_lb) * 0.7)
    else
        diff = 0
    end if
    print *, "diffp:", diff
    print *, "dq:", dq

    do lb2 = max_lb(2), min_ub(2) - diff(2), dq
    do ub2 = lb2 + diff(2), min_ub(2), dq
    !$omp parallel do
    do lb3 = max_lb(3), min_ub(3) - diff(3), dq
    do ub3 = lb3 + diff(3), min_ub(3), dq
    do lb4 = max_lb(4), min_ub(4) - diff(4), dq
    do ub4 = lb4 + diff(4), min_ub(4), dq
      if ( check_area(dq, ei, ee, lb2, ub2, lb3, ub3, lb4, ub4, max_lb, min_ub) ) then
        v = (ee - ei + 1)*(ub2 - lb2 + 1)*(ub3 - lb3 + 1)*(ub4 - lb4 + 1)
        if (real(sum(cond(ei:ee, lb2:ub2, lb3:ub3, lb4:ub4))) >= a_v * v .and. v > local_maxv(lb2, lb3, lb4) ) then
          local_maxv(lb2, lb3, lb4) = v
          local_argmaxv_lb(:, lb2, lb3, lb4) = (/ei, lb2, lb3, lb4/)
          local_argmaxv_ub(:, lb2, lb3, lb4) = (/ee, ub2, ub3, ub4/)
        end if
      end if
    end do
    end do
    end do
    end do

    end do
    end do
    maxlb4 = maxloc(local_maxv)
    if ( maxval(local_maxv) > 0 ) then
      maxlb4 = maxloc(local_maxv)
      argmaxv_lb = local_argmaxv_lb(:, maxlb4(1), maxlb4(2), maxlb4(3))
      argmaxv_ub = local_argmaxv_ub(:, maxlb4(1), maxlb4(2), maxlb4(3))
    end if
    print *, argmaxv_lb
    print *, argmaxv_ub
  end subroutine argmaxv


  function check_area(dq, elb, eub, lb2, ub2, lb3, ub3, lb4, ub4, max_lb, min_ub)
    integer, intent(in) :: dq, elb, eub, lb2, ub2, lb3, ub3, lb4, ub4
    integer, intent(in) :: max_lb(4), min_ub(4)
    logical check_area
    check_area = .false.

    if ( real(sum(cond(elb:eub, max(lb2-1,1):min(lb2+1,datasize(2)), lb3:ub3, lb4:ub4))) >= &
        a_b * (eub-elb+1)*(ub3-lb3+1)*(ub4-lb4+1)*(min(lb2+1,datasize(2))-max(lb2-1,1)+1) ) then  
    if ( real(sum(cond(elb:eub, max(ub2-1,1):min(ub2+1,datasize(2)), lb3:ub3, lb4:ub4))) >= &
        a_b * (eub-elb+1)*(ub3-lb3+1)*(ub4-lb4+1)*(min(ub2+1,datasize(2))-max(ub2-1,1)+1) ) then  
    if ( real(sum(cond(elb:eub, lb2:ub2, max(lb3-1,1):min(lb3+1,datasize(3)), lb4:ub4))) >= &
        a_b * (eub-elb+1)*(ub2-lb2+1)*(ub4-lb4+1)*(min(lb3+1,datasize(3))-max(lb3-1,1)+1) ) then  
    if ( real(sum(cond(elb:eub, lb2:ub2, max(ub3-1,1):min(ub3+1,datasize(3)), lb4:ub4))) >= &
        a_b * (eub-elb+1)*(ub2-lb2+1)*(ub4-lb4+1)*(min(ub3+1,datasize(3))-max(ub3-1,1)+1) ) then  
    if ( real(sum(cond(elb:eub, lb2:ub2, lb3:ub3, max(lb4-1,1):min(lb4+1,datasize(4))))) >= &
        a_b * (eub-elb+1)*(ub2-lb2+1)*(ub3-lb3+1)*(min(lb4+1,datasize(4))-max(lb4-1,1)+1) ) then  
    if ( real(sum(cond(elb:eub, lb2:ub2, lb3:ub3, max(ub4-1,1):min(ub4+1,datasize(4))))) >= &
        a_b * (eub-elb+1)*(ub2-lb2+1)*(ub3-lb3+1)*(min(ub4+1,datasize(4))-max(ub4-1,1)+1) ) then  
    !if ( real(sum(cond(elb:eub, lb2, lb3:ub3, lb4:ub4))) >= a_b * (eub-elb+1)*(ub3-lb3+1)*(ub4-lb4+1) ) then  
    !if ( real(sum(cond(elb:eub, ub2, lb3:ub3, lb4:ub4))) >= a_b * (eub-elb+1)*(ub3-lb3+1)*(ub4-lb4+1) ) then  
    !if ( real(sum(cond(elb:eub, lb2:ub2, lb3, lb4:ub4))) >= a_b * (eub-elb+1)*(ub2-lb2+1)*(ub4-lb4+1) ) then  
    !if ( real(sum(cond(elb:eub, lb2:ub2, ub3, lb4:ub4))) >= a_b * (eub-elb+1)*(ub2-lb2+1)*(ub4-lb4+1) ) then  
    !if ( real(sum(cond(elb:eub, lb2:ub2, lb3:ub3, lb4))) >= a_b * (eub-elb+1)*(ub2-lb2+1)*(ub3-lb3+1) ) then  
    !if ( real(sum(cond(elb:eub, lb2:ub2, lb3:ub3, ub4))) >= a_b * (eub-elb+1)*(ub2-lb2+1)*(ub3-lb3+1) ) then
      check_area = .true.
    end if
    end if
    end if
    end if
    end if
    end if
  end function check_area


  subroutine argmaxvlight(eidx, dq, max_lb, min_ub, argmaxv_lb, argmaxv_ub)
    integer, intent(in)  :: eidx, dq
    integer, intent(in)  :: max_lb(4), min_ub(4)
    integer, intent(inout) :: argmaxv_lb(4), argmaxv_ub(4)
    integer  lb2, lb3, lb4, ub2, ub3, ub4, v, maxlb4(3), diff(4)
    integer  local_maxv(datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_lb(4, datasize(2), datasize(3), datasize(4))
    integer  local_argmaxv_ub(4, datasize(2), datasize(3), datasize(4))
    local_maxv = 0

    if ( sum(argmaxv_ub - argmaxv_lb) /=0 ) then
        diff = int((argmaxv_ub - argmaxv_lb) * 0.7)
    else
        diff = 0
    end if
    print *, "diff:", diff

    do lb2 = max_lb(2), min_ub(2) - diff(2), dq
    do ub2 = lb2 + diff(2), min_ub(2), dq
    !$omp parallel do
    do lb3 = max_lb(3), min_ub(3) - diff(3), dq
    do ub3 = lb3 + diff(3), min_ub(3), dq
    do lb4 = max_lb(4), min_ub(4) - diff(4), dq
    do ub4 = lb4 + diff(4), min_ub(4), dq
      if ( check_area(dq, eidx, eidx, lb2, ub2, lb3, ub3, lb4, ub4, max_lb, min_ub) ) then 
        v = (ub2 - lb2 + 1)*(ub3 - lb3 + 1)*(ub4 - lb4 + 1)
        if (real(sum(cond(eidx, lb2:ub2, lb3:ub3, lb4:ub4))) >= a_v * v .and. v > local_maxv(lb2, lb3, lb4) ) then
          local_maxv(lb2, lb3, lb4) = v
          local_argmaxv_lb(:, lb2, lb3, lb4) = (/eidx, lb2, lb3, lb4/)
          local_argmaxv_ub(:, lb2, lb3, lb4) = (/eidx, ub2, ub3, ub4/)
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


end module costfort4d
