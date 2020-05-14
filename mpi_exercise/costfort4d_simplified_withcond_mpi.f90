!fortran 90 library for 4D bin width optimization
!this library is used by optbinwidth4d_wholefort.py
!Kazuyoshi TATSUMI 2020/01/23
!mpi + openmpi version  2020/05/13

module costfort4d
  use mpi
  implicit none
  integer datasize(4)
  integer maxw(4)
contains

  subroutine cost4d(comm, maxw3, maxw2, maxw1, maxw0, datasize3, datasize2, datasize1, datasize0,& 
      usecond, condparam, data_array, condition, cost, histaves, deltas) bind(C)
    integer, intent(in) :: maxw0
    integer, intent(in) :: maxw1
    integer, intent(in) :: maxw2
    integer, intent(in) :: maxw3
    integer, intent(in) :: datasize0
    integer, intent(in) :: datasize1
    integer, intent(in) :: datasize2
    integer, intent(in) :: datasize3
    double precision, intent(in) :: data_array(datasize0, datasize1, datasize2, datasize3)                         
    integer, intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
    logical(1), intent(in) :: usecond
    double precision, intent(in) :: condparam
    double precision, intent(inout), dimension(maxw0, maxw1, maxw2, maxw3) :: cost, histaves, deltas                         
    double precision, allocatable :: cumdata(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)                    
    integer, allocatable :: hist_cond(:,:,:,:)
    integer comm, psize, rank, ierr, mrank
    integer histsize(4), nw(4)
    integer ax_id1, ax_id2, ax_id3, width_id1, width_id2, width_id3, width_id4
    datasize = (/datasize0, datasize1, datasize2, datasize3/)
    maxw = (/maxw0, maxw1, maxw2, maxw3/)
    cost = 0.0
    histaves = 0.0
    deltas = 0.0

    call MPI_Comm_size(comm, psize, ierr)
    call MPI_Comm_rank(comm, rank, ierr)

    !!You can comment out or remove comment notations  to generate a *.so file which provides a part of the whole processes.
    !!In the latest form, I created following *.so files: 
    !!  1) costfort01d.so for help0d and help1d routines 
    !!  2) each of costfort2d1-6.so for each of the if sentences within double do-loops of help2d
    !!  3) each of costfort3d1-3.so for each of the if sentences within triple do-loops of help3d 
    !!  4) costfort4d-first,second.so for the two divisions of the outermost of the quintet do-loops
    !!These *.so co-operating with optbinwidth**.py enables pseudo-parallel calculations.
    !! 2020 Jan 16 Kazuyoshi TATSUMI
    print *, "entering help0d at ", rank
    call help0d(data_array, condition, usecond, condparam, cost, histaves, deltas)
    print *, "entering help1d at ", rank
    do ax_id1 = 1, 4
      call help1d(data_array, ax_id1, condition, usecond, condparam, cost, histaves, deltas, rank, psize)
    enddo
    print *, "entering help2d at ", rank
    do ax_id1 = 1, 3
    do ax_id2 = ax_id1 + 1, 4
      !if (ax_id1 == 1 .and. ax_id2 == 2) call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
      !if (ax_id1 == 1 .and. ax_id2 == 3) call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
      !if (ax_id1 == 1 .and. ax_id2 == 4) call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
      !if (ax_id1 == 2 .and. ax_id2 == 3) call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
      !if (ax_id1 == 1 .and. ax_id2 == 4) call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
      !if (ax_id1 == 3 .and. ax_id2 == 4) call help2d(data_array, [ax_id1, ax_id2], condition, usecond)
      call help2d(data_array, [ax_id1, ax_id2], condition, usecond, condparam, cost, histaves, deltas, rank, psize)
    enddo
    enddo
    print *, "entering help3d at ", rank
    do ax_id1 = 1, 2
    do ax_id2 = ax_id1 + 1, 3
    do ax_id3 = ax_id2 + 1, 4
      !if (ax_id1 == 1 .and. ax_id2 == 2 .and. ax_id3 == 3) call help3d(data_array, [ax_id1, ax_id2, ax_id3], condition, usecond)
      !if (ax_id1 == 1 .and. ax_id2 == 2 .and. ax_id3 == 4) call help3d(data_array, [ax_id1, ax_id2, ax_id3], condition, usecond)
      !if (ax_id1 == 1 .and. ax_id2 == 3 .and. ax_id3 == 4) call help3d(data_array, [ax_id1, ax_id2, ax_id3], condition, usecond)
      !if (ax_id1 == 2 .and. ax_id2 == 3 .and. ax_id3 == 4) call help3d(data_array, [ax_id1, ax_id2, ax_id3], condition, usecond)
      call help3d(data_array, [ax_id1, ax_id2, ax_id3], condition, usecond, condparam, cost, histaves, deltas, rank, psize)
    enddo
    enddo
    enddo



    print *, "entering hist4d at ", rank
    !cumdata = cumsum4d(cumsum4d(cumsum4d(cumsum4d(data_array, 1), 2), 3), 4)
    cumdata = cumsum4d(data_array, 1)
    cumdata = cumsum4d(cumdata, 2)
    cumdata = cumsum4d(cumdata, 3)
    cumdata = cumsum4d(cumdata, 4)
    if (usecond) then 
        !cum_cond = cumsum4di(cumsum4di(cumsum4di(cumsum4di(condition, 1), 2), 3), 4)
        cum_cond = cumsum4di(condition, 1)
        cum_cond = cumsum4di(cum_cond, 2)
        cum_cond = cumsum4di(cum_cond, 3)
        cum_cond = cumsum4di(cum_cond, 4)
    endif

    !call MPI_Barrier(comm, ierr)
    !if (maxw(1) >= rank+2) then
    mrank = psize - rank - 1
    if (maxw(1) >= mrank+2) then
    !do width_id1 = rank+2, maxw(1) - mod(maxw(1)-rank-2, psize), psize
    do width_id1 = mrank+2, maxw(1) - mod(maxw(1)-mrank-2, psize), psize
    !do width_id1 = 1, maxw(1)
    nw(1) = width_id1
    !print *, "nw(1):", nw(1), "at rank=", rank
    histsize(1) = (datasize(1) - mod(datasize(1), nw(1))) / nw(1)
    do width_id2 = 2, maxw(2)
    nw(2) = width_id2
    histsize(2) = (datasize(2) - mod(datasize(2), nw(2))) / nw(2)
    do width_id3 = 2, maxw(3)
    nw(3) = width_id3
    histsize(3) = (datasize(3) - mod(datasize(3), nw(3))) / nw(3)
    do width_id4 = 2, maxw(4)
    nw(4) = width_id4
    histsize(4) = (datasize(4) - mod(datasize(4), nw(4))) / nw(4)
          hist_array = hist4d(cumdata, histsize, nw)
          if (usecond) then 
              hist_cond = hist4di(cum_cond, histsize, nw)
              !call stat1(pack(hist_array, hist_cond == maxval(hist_cond)), nw)
              call stat1(pack(hist_array, hist_cond >= maxval(hist_cond)*condparam), nw, cost, histaves, deltas)
          else
              call stat(hist_array, nw, cost, histaves, deltas) 
          endif
    end do
    end do
    end do
    end do
    end if

    print *, "minloc Cn:", minloc(cost), "with its value:", minval(cost), "at ", rank

  end subroutine cost4d

  subroutine help0d(data_array, condition, usecond, condparam, cost, histaves, deltas)
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    double precision, intent(inout), dimension(maxw(1), maxw(2), maxw(3), maxw(4)) :: cost, histaves, deltas                         
    double precision, intent(in) :: condparam
    integer nw(4)
    nw = 1
    if (usecond) then
      !call stat1(pack(hist_array, condition == maxval(condition)), nw)
      call stat1(pack(data_array, condition >= maxval(condition)*condparam), nw, cost, histaves, deltas)
    else
      call stat(data_array, nw, cost, histaves, deltas)
    end if
  end subroutine help0d
  subroutine help1d(data_array, ax, condition, usecond, condparam,  cost, histaves, deltas, rank, psize)
    integer, intent(in) :: ax
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    double precision, intent(in) :: condparam
    double precision, intent(inout), dimension(maxw(1), maxw(2), maxw(3), maxw(4)) :: cost, histaves, deltas                         
    logical(1) :: Ishistsizeallocated
    double precision, allocatable :: cumdata1d(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)
    integer, allocatable :: histcond(:,:,:,:)
    integer width_id, histsize, nw(4),  cumdata_size(4), cumdata_order(4)
    integer psize, rank
    nw = 1
    Ishistsizeallocated = .true.
    cumdata1d = cumsum4d(data_array, ax)
    if (usecond) then
        cum_cond = cumsum4di(condition, ax)
    endif
    if (ax == 1) then
      Ishistsizeallocated = .false.
    elseif (ax == 2) then
      cumdata_size = (/datasize(2), datasize(1), datasize(3), datasize(4)/) 
      cumdata_order = (/2, 1, 3, 4/)
    elseif (ax == 3) then
      cumdata_size = (/datasize(3), datasize(2), datasize(1), datasize(4)/) 
      cumdata_order = (/3, 2, 1, 4/)
    elseif (ax == 4) then
      cumdata_size = (/datasize(4), datasize(2), datasize(3), datasize(1)/) 
      cumdata_order = (/4, 2, 3, 1/)
    endif
    if (Ishistsizeallocated) then
      cumdata1d = reshape(cumdata1d, cumdata_size, order = cumdata_order)
      if (usecond) then
        cum_cond = reshape(cum_cond, cumdata_size, order = cumdata_order)
      endif
    endif
    !do width_id = 2, maxw(ax)
    if (maxw(ax) >= rank+2) then
    do width_id = rank+2, maxw(ax) - mod(maxw(ax)-rank-2, psize), psize
      nw(ax) = width_id
      histsize = (datasize(ax) - mod(datasize(ax), nw(ax))) / nw(ax)
      hist_array = hist1d(cumdata1d, histsize, nw(ax))
      if (usecond) then
        histcond = hist1di(cum_cond, histsize, nw(ax))
        !call stat1(pack(hist_array, histcond == maxval(histcond)), nw)
        call stat1(pack(hist_array, histcond >= maxval(histcond)*condparam), nw, cost, histaves, deltas)
      else
        call stat(hist_array, nw, cost, histaves, deltas)
      end if
    enddo
    end if
  end subroutine help1d

  subroutine help2d(data_array, ax, condition, usecond, condparam,  cost, histaves, deltas, rank, psize)
    integer, intent(in) :: ax(:)
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    double precision, intent(in) :: condparam
    double precision, intent(inout), dimension(maxw(1), maxw(2), maxw(3), maxw(4)) :: cost, histaves, deltas                         
    logical(1) :: Ishistsizeallocated
    double precision, allocatable :: cumdata2d(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)
    integer, allocatable :: histcond(:,:,:,:)
    integer nw(4), histsize(4), cumdata_size_id, width_id1, width_id2
    integer, allocatable ::  cumdata_size(:,:), cumdata_order(:,:)
    integer psize, rank, mrank
    nw = 1
    !cumdata2d = cumsum4d(cumsum4d(data_array, ax(1)), ax(2))
    cumdata2d = cumsum4d(data_array, ax(1))
    cumdata2d = cumsum4d(cumdata2d, ax(2))
    if (usecond) then 
        !cum_cond = cumsum4di(cumsum4di(condition, ax(1)), ax(2))
        cum_cond = cumsum4di(condition, ax(1))
        cum_cond = cumsum4di(cum_cond, ax(2))
    endif
    if (ax(1) == 1 .and. ax(2) == 2) then
      Ishistsizeallocated = .false.
    else
      Ishistsizeallocated = .true.
      if (ax(1) == 1 .and. ax(2) == 3) then
        allocate(cumdata_size(1,4), cumdata_order(1,4))
        cumdata_size(1,:) = (/datasize(1), datasize(3), datasize(2), datasize(4)/)
        cumdata_order(1,:) = (/1, 3, 2, 4/)
      elseif (ax(1) == 1 .and. ax(2) == 4) then
        allocate(cumdata_size(1,4), cumdata_order(1,4))
        cumdata_size(1,:) = (/datasize(1), datasize(4), datasize(3), datasize(2)/)
        cumdata_order(1,:) = (/1, 4, 3, 2/)
      elseif (ax(1) == 2 .and. ax(2) == 3) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(2), datasize(1), datasize(3), datasize(4)/)
        cumdata_order(1,:) = (/2, 1, 3, 4/)
        cumdata_size(2,:) = (/datasize(2), datasize(3), datasize(1), datasize(4)/)
        cumdata_order(2,:) = (/1, 3, 2, 4/)
      elseif (ax(1) == 2 .and. ax(2) == 4) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(2), datasize(1), datasize(3), datasize(4)/)
        cumdata_order(1,:) = (/2, 1, 3, 4/)
        cumdata_size(2,:) = (/datasize(2), datasize(4), datasize(3), datasize(1)/)
        cumdata_order(2,:) = (/1, 4, 3, 2/)
      elseif (ax(1) == 3 .and. ax(2) == 4) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(3), datasize(2), datasize(1), datasize(4)/)
        cumdata_order(1,:) = (/3, 2, 1, 4/)
        cumdata_size(2,:) = (/datasize(3), datasize(4), datasize(1), datasize(2)/)
        cumdata_order(2,:) = (/1, 4, 3, 2/)
      endif 
    endif
    if (Ishistsizeallocated) then
      do cumdata_size_id = 1, size(cumdata_size,1)
         cumdata2d = reshape(cumdata2d, cumdata_size(cumdata_size_id,1:4), order = cumdata_order(cumdata_size_id,1:4))
         if (usecond) cum_cond = reshape(cum_cond, cumdata_size(cumdata_size_id,1:4), &
                                         order = cumdata_order(cumdata_size_id, 1:4))
      enddo
    endif
    !do width_id1 = 2, maxw(ax(1))
    !mrank = psize - rank -1
    if ( maxw(ax(1)) >=  rank  + 2 ) then
    !if ( maxw(ax(1)) >=  mrank  + 2 ) then
    do width_id1 = rank+2, maxw(ax(1)) - mod(maxw(ax(1))-rank-2, psize), psize
    !do width_id1 = mrank+2, maxw(ax(1)) - mod(maxw(ax(1))-mrank-2, psize), psize
      if (maxw(ax(1)) < width_id1 ) print *, "error width_id1", width_id1," exceeds maxw(ax(1))", maxw(ax(1)), ax(1)
      nw(ax(1)) = width_id1
    do width_id2 = 2, maxw(ax(2))
      nw(ax(2)) = width_id2
      histsize = (datasize - mod(datasize, nw)) / nw
      hist_array = hist2d(cumdata2d, [histsize(ax(1)), histsize(ax(2))], [nw(ax(1)), nw(ax(2))])
      if (usecond) then
        histcond = hist2di(cum_cond, [histsize(ax(1)), histsize(ax(2))], [nw(ax(1)), nw(ax(2))])
        !call stat1(pack(hist_array, histcond == maxval(histcond)), nw)
        call stat1(pack(hist_array, histcond >= maxval(histcond)*condparam), nw, cost, histaves, deltas)
      else
        call stat(hist_array, nw, cost, histaves, deltas)
      endif
    enddo
    enddo
    endif
  end subroutine help2d

  subroutine help3d(data_array, ax, condition, usecond, condparam, cost, histaves, deltas, rank, psize)
    integer, intent(in) :: ax(:)
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    integer, intent(in) :: condition(datasize(1), datasize(2), datasize(3), datasize(4))
    logical(1), intent(in) :: usecond
    double precision, intent(in) :: condparam
    double precision, intent(inout), dimension(maxw(1), maxw(2), maxw(3), maxw(4)) :: cost, histaves, deltas                         
    logical(1) :: Ishistsizeallocated
    double precision, allocatable :: cumdata3d(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)
    integer, allocatable :: histcond(:,:,:,:)
    integer nw(4), histsize(4), cumdata_size_id, width_id1, width_id2, width_id3
    integer, allocatable ::  cumdata_size(:,:),cumdata_order(:,:)
    integer psize, rank
    nw = 1
    !cumdata3d = cumsum4d(cumsum4d(cumsum4d(data_array, ax(1)), ax(2)), ax(3))
    cumdata3d = cumsum4d(data_array, ax(1))
    cumdata3d = cumsum4d(cumdata3d, ax(2))
    cumdata3d = cumsum4d(cumdata3d, ax(3))
    if (usecond) then 
        !cum_cond = cumsum4di(cumsum4di(cumsum4di(condition, ax(1)), ax(2)), ax(3))
        cum_cond = cumsum4di(condition, ax(1))
        cum_cond = cumsum4di(cum_cond, ax(2))
        cum_cond = cumsum4di(cum_cond, ax(3))
    endif
    if (ax(1) == 1 .and. ax(2) == 2 .and. ax(3) == 3) then
      Ishistsizeallocated = .false.
    else
      Ishistsizeallocated = .true.
      if (ax(1) == 1 .and. ax(2) == 2 .and. ax(3) == 4) then
        allocate(cumdata_size(1,4), cumdata_order(1,4))
        cumdata_size(1,:) = (/datasize(1), datasize(2), datasize(4), datasize(3)/)
        cumdata_order(1,:) = (/1, 2, 4, 3/)
      elseif (ax(1) == 1 .and. ax(2) == 3 .and. ax(3) == 4) then
        allocate(cumdata_size(2,4), cumdata_order(2,4))
        cumdata_size(1,:) = (/datasize(1), datasize(3), datasize(2), datasize(4)/)
        cumdata_order(1,:) = (/1, 3, 2, 4/)
        cumdata_size(2,:) = (/datasize(1), datasize(3), datasize(4), datasize(2)/)
        cumdata_order(2,:) = (/1, 2, 4, 3/)
      elseif (ax(1) == 2 .and. ax(2) == 3 .and. ax(3) == 4) then
        allocate(cumdata_size(3,4), cumdata_order(3,4))
        cumdata_size(1,:) = (/datasize(2), datasize(1), datasize(3), datasize(4)/)
        cumdata_order(1,:) = (/2, 1, 3, 4/)
        cumdata_size(2,:) = (/datasize(2), datasize(3), datasize(1), datasize(4)/)
        cumdata_order(2,:) = (/1, 3, 2, 4/)
        cumdata_size(3,:) = (/datasize(2), datasize(3), datasize(4), datasize(1)/)
        cumdata_order(3,:) = (/1, 2, 4, 3/)
      end if
    endif
    if (Ishistsizeallocated) then
      do cumdata_size_id = 1, size(cumdata_size,1)
        cumdata3d = reshape(cumdata3d, cumdata_size(cumdata_size_id,1:4), order = cumdata_order(cumdata_size_id,1:4))
        if (usecond) cum_cond = reshape(cum_cond, cumdata_size(cumdata_size_id,1:4), order = cumdata_order(cumdata_size_id,1:4))
      enddo
    endif
    !do width_id1 = 2, maxw(ax(1))
    if (maxw(ax(1)) >= rank+2) then
    do width_id1 = rank+2, maxw(ax(1)) - mod(maxw(ax(1))-rank-2, psize), psize
      nw(ax(1)) = width_id1
    do width_id2 = 2, maxw(ax(2))
      nw(ax(2)) = width_id2
    do width_id3 = 2, maxw(ax(3))
      nw(ax(3)) = width_id3
      histsize = (datasize - mod(datasize, nw)) / nw
      hist_array = hist3d(cumdata3d, [histsize(ax(1)), histsize(ax(2)), histsize(ax(3))], [nw(ax(1)), nw(ax(2)), nw(ax(3))])
      if (usecond) then
        histcond = hist3di(cum_cond, [histsize(ax(1)), histsize(ax(2)), histsize(ax(3))], [nw(ax(1)), nw(ax(2)), nw(ax(3))])
        !call stat1(pack(hist_array, histcond == maxval(histcond)), nw)
        call stat1(pack(hist_array, histcond >= maxval(histcond)*condparam), nw, cost, histaves, deltas)
      else
        call stat(hist_array, nw, cost, histaves, deltas) 
      endif
    enddo
    enddo
    enddo
    endif
  end subroutine help3d

  subroutine stat(k, nw, cost, histaves, deltas)
    double precision, intent(in) :: k(:,:,:,:)
    integer, intent(in) :: nw(4)
    double precision, intent(inout), dimension(maxw(1), maxw(2), maxw(3), maxw(4)) :: cost, histaves, deltas                         
    double precision kave, v
    kave = sum(k) / dble(size(k)); v = sum((k - kave)**2) / dble(size(k))
    !print *, "cost with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / (dble(product(nw))**2)
    histaves(nw(1), nw(2), nw(3), nw(4)) = kave
    deltas(nw(1), nw(2), nw(3), nw(4)) = product(nw)
    cost(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / (dble(product(nw))**2)
  end subroutine stat
  subroutine stat1(knonzero, nw, cost, histaves, deltas)
    double precision, intent(in) :: knonzero(:)
    integer, intent(in) :: nw(4)
    double precision, intent(inout), dimension(maxw(1), maxw(2), maxw(3), maxw(4)) :: cost, histaves, deltas                         
    double precision kave, v
    kave = sum(knonzero) / dble(size(knonzero)); v = sum((knonzero - kave)**2) / dble(size(knonzero))
    !print *, "cost with ", nw(1), nw(2), nw(3), nw(4), ":", (2.0 * kave - v) / (dble(product(nw))**2)
    histaves(nw(1), nw(2), nw(3), nw(4)) = kave
    deltas(nw(1), nw(2), nw(3), nw(4)) = product(nw)
    cost(nw(1), nw(2), nw(3), nw(4)) = (2.0 * kave - v) / (dble(product(nw))**2)
  end subroutine stat1
          

  function cumsum4d(data_array, ax)
    integer, intent(in) :: ax
    double precision, intent(in) :: data_array(datasize(1), datasize(2), datasize(3), datasize(4))
    double precision  :: cumsum4d(datasize(1), datasize(2), datasize(3),datasize(4))
    integer :: cum_id
    cumsum4d = data_array
    if (ax == 4) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(:, :, :, cum_id + 1) = cumsum4d(:, :, :, cum_id) + data_array(:, :, :, cum_id + 1)
       end do
    else if (ax == 3) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(:, :, cum_id + 1, :) = cumsum4d(:, :, cum_id, :) + data_array(:, :, cum_id + 1, :)
       end do
    else if (ax == 2) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(:, cum_id + 1, :, :) = cumsum4d(:, cum_id, :, :) + data_array(:, cum_id + 1, :, :)
       end do
    else if (ax == 1) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4d(cum_id + 1, :, :, :) = cumsum4d(cum_id, :, :, :) + data_array(cum_id + 1, :, :, :)
       end do
    end if
  end function cumsum4d       
  function cumsum4di(d, ax)
    integer, intent(in) :: ax
    integer, intent(in) :: d(datasize(1), datasize(2), datasize(3), datasize(4))
    integer :: cumsum4di(datasize(1), datasize(2), datasize(3), datasize(4))
    integer :: cum_id
    cumsum4di = d
    if (ax == 4) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(:, :, :, cum_id + 1) = cumsum4di(:, :, :, cum_id) + d(:, :, :, cum_id + 1)
       end do
    else if (ax == 3) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(:, :, cum_id + 1, :) = cumsum4di(:, :, cum_id, :) + d(:, :, cum_id + 1, :)
       end do
    else if (ax == 2) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(:, cum_id + 1, :, :) = cumsum4di(:, cum_id, :, :) + d(:, cum_id + 1, :, :)
       end do
    else if (ax == 1) then
       do cum_id = 1, datasize(ax) - 1
          cumsum4di(cum_id + 1, :, :, :) = cumsum4di(cum_id, :, :, :) + d(cum_id + 1, :, :, :)
       end do
    end if
  end function cumsum4di

  function hist4d(cumdata, N, nw)
    double precision, intent(in) :: cumdata(:,:,:,:) 
    integer, intent(in) :: N(4), nw(4)
    integer :: i, j, h, l, ihead, jhead, hhead, lhead
    double precision :: hist4d(N(1), N(2), N(3),N(4))
    !$omp parallel do
    do i = 1, N(1)
    ihead = i*nw(1) 
    do j = 1, N(2)
    jhead = j*nw(2) 
    do h = 1, N(3)
    hhead = h*nw(3) 
    do l = 1, N(4)
    lhead = l*nw(4) 
      if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead)
      else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead - nw(1), jhead, hhead, lhead)
      else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead)
      else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead)
      else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4))
      else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead)
      else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead)
      else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4))
      else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead)
      else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4))
      else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4))
      else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead)
      else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4))
      else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4))
      else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                       - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                       + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4))
      else
         hist4d(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                       - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                       - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                       + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4)) &
                       - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead) &
                       + cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead - nw(4))
      end if
    end do
    end do
    end do
    end do
  end function hist4d

  function hist4di(cumdata, N, nw)
    integer, intent(in) :: cumdata(:,:,:,:) 
    integer, intent(in) :: N(4), nw(4)
    integer :: i, j, h, l, ihead, jhead, hhead, lhead
    integer :: hist4di(N(1),N(2),N(3),N(4))
    !$omp parallel do
    do i = 1, N(1)
    ihead = i*nw(1) 
    do j = 1, N(2)
    jhead = j*nw(2) 
    do h = 1, N(3)
    hhead = h*nw(3) 
    do l = 1, N(4)
    lhead = l*nw(4) 
       if ( i == 1 .and. j == 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead)
       else if ( j == 1 .and. i /= 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead - nw(1), jhead, hhead, lhead)
       else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead)
       else if ( j == 1 .and. h /= 1 .and. i == 1 .and. l == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead)
       else if ( j == 1 .and. l /= 1 .and. h == 1 .and. i == 1 ) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4))
       else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead)
       else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead)
       else if ( i /= 1 .and. j == 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4))
       else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead)
       else if ( i == 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4))
       else if ( i == 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4))
       else if ( i /= 1 .and. j /= 1 .and. h /= 1 .and. l == 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead)
       else if ( i /= 1 .and. j /= 1 .and. h == 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4))
       else if ( i /= 1 .and. j == 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4))
       else if ( i == 1 .and. j /= 1 .and. h /= 1 .and. l /= 1) then
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead, jhead - nw(2), hhead, lhead) - cumdata(ihead, jhead, hhead - nw(3), lhead) &
                        - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4))
       else
          hist4di(i, j, h, l) = cumdata(ihead, jhead, hhead, lhead) &
                        - cumdata(ihead - nw(1), jhead, hhead, lhead) - cumdata(ihead, jhead - nw(2), hhead, lhead) &
                        - cumdata(ihead, jhead, hhead - nw(3), lhead) - cumdata(ihead, jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead, hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead, jhead - nw(2), hhead, lhead - nw(4)) &
                        + cumdata(ihead, jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead, jhead - nw(2), hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead, hhead - nw(3), lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead, lhead - nw(4)) &
                        - cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead) &
                        + cumdata(ihead - nw(1), jhead - nw(2), hhead - nw(3), lhead - nw(4))
       end if
    end do
    end do
    end do
    end do
  end function hist4di

  function hist1d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N, nw
    integer :: i, ihead
    double precision :: hist1d(N, size(B,2), size(B,3), size(B,4))
    hist1d(1,:,:,:) = B(nw,:,:,:)
    !$omp parallel do
    do i = 2, N
       ihead = i*nw
       hist1d(i,:,:,:) = B(ihead,:,:,:) - B(ihead - nw,:,:,:)
    end do
  end function hist1d
  function hist1di(C, N, nw)
    integer, intent(in) :: C(:,:,:,:) 
    integer, intent(in) :: N, nw
    integer :: i, ihead
    integer :: hist1di(N, size(C,2), size(C,3), size(C,4))
    hist1di(1,:,:,:) = C(nw,:,:,:)
    !$omp parallel do
    do i = 2, N
       ihead = i*nw
       hist1di(i,:,:,:) = C(ihead,:,:,:) - C(ihead - nw,:,:,:)
    end do
  end function hist1di

  function hist2d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(2), nw(2)
    integer :: i, ihead, j, jhead
    double precision :: hist2d(N(1), N(2), size(B,3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2) 
          if ( i == 1 .and. j == 1) then
              hist2d(i, j, :, :) = B(ihead, jhead,:, :)
          else if ( i /= 1 .and. j == 1) then
              hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :)
          else if ( i == 1 .and. j /= 1) then
              hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead, jhead - nw(2), :, :)
          else 
              hist2d(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :) - B(ihead, jhead - nw(2), :, :) &
                                   + B(ihead - nw(1), jhead - nw(2), :, :)
          end if
       end do
    end do
  end function hist2d
  function hist2di(B, N, nw)
    integer, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(2), nw(2)
    integer :: i, ihead, j, jhead
    integer :: hist2di(N(1), N(2), size(B,3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2) 
          if ( i == 1 .and. j == 1) then
              hist2di(i, j, :, :) = B(ihead, jhead,:, :)
          else if ( i /= 1 .and. j == 1) then
              hist2di(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :)
          else if ( i == 1 .and. j /= 1) then
              hist2di(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead, jhead - nw(2), :, :)
          else 
              hist2di(i, j, :, :) = B(ihead, jhead, :, :) - B(ihead - nw(1), jhead, :, :) - B(ihead, jhead - nw(2), :, :) &
                                   + B(ihead - nw(1), jhead - nw(2), :, :)
          end if
       end do
    end do
  end function hist2di

  function hist3d(B, N, nw)
    double precision, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(3), nw(3)
    integer :: i, ihead, j, jhead, h, hhead
    double precision :: hist3d(N(1), N(2), N(3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2)
          do h = 1, N(3)
             hhead = h*nw(3)
             if ( i == 1 .and. j == 1 .and. h == 1) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :)
             else if ( j == 1 .and. i /= 1 .and. h == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :)
             else if ( i == 1 .and. j /= 1 .and. h == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :)
             else if ( i == 1 .and. h /= 1 .and. j == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead, hhead - nw(3), :)
             else if ( i /= 1 .and. j /= 1 .and. h == 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                     - B(ihead, jhead - nw(2), hhead, :) + B(ihead - nw(1), jhead - nw(2), hhead, :)
             else if ( i /= 1 .and. j == 1 .and. h /= 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :)
             else if ( i == 1 .and. j /= 1 .and. h /= 1 ) then
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead, jhead - nw(2), hhead - nw(3), :)
             else
                hist3d(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead - nw(2), hhead, :) - B(ihead, jhead, hhead - nw(3), :) &
                                    + B(ihead, jhead - nw(2), hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :) &
                                    + B(ihead - nw(1), jhead - nw(2), hhead, :) - B(ihead - nw(1), jhead - nw(2), hhead - nw(3), :)
             end if
          end do
       end do
    end do
  end function hist3d
  function hist3di(B, N, nw)
    integer, intent(in) :: B(:,:,:,:) 
    integer, intent(in) :: N(3), nw(3)
    integer :: i, ihead, j, jhead, h, hhead
    integer :: hist3di(N(1), N(2), N(3), size(B,4))
    !$omp parallel do
    do i = 1, N(1)
       ihead = i*nw(1) 
       do j = 1, N(2)
          jhead = j*nw(2)
          do h = 1, N(3)
             hhead = h*nw(3)
             if ( i == 1 .and. j == 1 .and. h == 1) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :)
             else if ( j == 1 .and. i /= 1 .and. h == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :)
             else if ( i == 1 .and. j /= 1 .and. h == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :)
             else if ( i == 1 .and. h /= 1 .and. j == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead, hhead - nw(3), :)
             else if ( i /= 1 .and. j /= 1 .and. h == 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                     - B(ihead, jhead - nw(2), hhead, :) + B(ihead - nw(1), jhead - nw(2), hhead, :)
             else if ( i /= 1 .and. j == 1 .and. h /= 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :)
             else if ( i == 1 .and. j /= 1 .and. h /= 1 ) then
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead, jhead - nw(2), hhead, :) &
                                    - B(ihead, jhead, hhead - nw(3), :) + B(ihead, jhead - nw(2), hhead - nw(3), :)
             else
                hist3di(i, j, h, :) = B(ihead, jhead, hhead, :) - B(ihead - nw(1), jhead, hhead, :) &
                                    - B(ihead, jhead - nw(2), hhead, :) - B(ihead, jhead, hhead - nw(3), :) &
                                    + B(ihead, jhead - nw(2), hhead - nw(3), :) + B(ihead - nw(1), jhead, hhead - nw(3), :) &
                                    + B(ihead - nw(1), jhead - nw(2), hhead, :) - B(ihead - nw(1), jhead - nw(2), hhead - nw(3), :)
             end if
          end do
       end do
    end do
  end function hist3di

end module costfort4d
