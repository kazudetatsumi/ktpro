!fortran 90 library for 4D bin width optimization
!this library is used by optbinwidth4d_wholefort.py
!Kazuyoshi TATSUMI 2020/01/23

module costfort4d
 use mpi
 implicit none
 integer datasize(4)
 integer nw(4)
contains

  subroutine cost4d(comm, maxnw3, maxnw2, maxnw1, maxnw0, datasize3, datasize2, datasize1, datasize0,& 
    sumpdf, data_array, pdfs, condition, ise) bind(C)
    integer, intent(in) :: maxnw0
    integer, intent(in) :: maxnw1
    integer, intent(in) :: maxnw2
    integer, intent(in) :: maxnw3
    integer, intent(in) :: datasize0
    integer, intent(in) :: datasize1
    integer, intent(in) :: datasize2
    integer, intent(in) :: datasize3
    double precision, intent(in) :: sumpdf
    double precision, intent(in) :: data_array(datasize0, datasize1, datasize2, datasize3)                         
    double precision, intent(in) :: pdfs(datasize0, datasize1, datasize2, datasize3)                         
    integer, intent(in) :: condition(datasize0, datasize1, datasize2, datasize3)                         
	double precision, intent(inout), dimension(maxnw0, maxnw1, maxnw2, maxnw3) :: ise
    double precision, allocatable :: cumdata(:,:,:,:)
    integer, allocatable :: cum_cond(:,:,:,:)
    double precision, allocatable :: hist_array(:,:,:,:)                    
	double precision sumdata, sumpdfs
	integer :: nws(maxnw0*maxnw1*maxnw2*maxnw3, 4)
    integer, allocatable :: hist_cond(:,:,:,:)
    integer comm, psize, rank, ierr, mrank
    integer histsize(4), sumhist, maxcond
    integer hx, hy, hz, hw, px, py, pz, pw, nw0, nw1, nw2, nw3, idx1d, maxw
    call MPI_Comm_size(comm, psize, ierr)
    call MPI_Comm_rank(comm, rank, ierr)
	ise(:,:,:,:) = 0.0
	idx1d = 0 
	do nw0 = 1, maxnw0
	do nw1 = 1, maxnw1
	do nw2 = 1, maxnw2
	do nw3 = 1, maxnw3
	idx1d = idx1d + 1
	nws(idx1d, 1) = nw0
	nws(idx1d, 2) = nw1
	nws(idx1d, 3) = nw2
	nws(idx1d, 4) = nw3
	enddo
	enddo
	enddo
	enddo



	!do idx1d = 1, maxnw0*maxnw1*maxnw2*maxnw3

	mrank = psize - rank - 1
	maxw = maxnw0*maxnw1*maxnw2*maxnw3
	if (maxw >= mrank+1) then
	do idx1d = mrank+1, maxw - mod(maxw-mrank-1, psize), psize

	nw = (/nws(idx1d, 1), nws(idx1d, 2), nws(idx1d, 3), nws(idx1d, 4)/)
	!nw = (/nw0, nw1, nw2, nw3/)
    datasize = (/datasize0, datasize1, datasize2, datasize3/)
	sumpdfs = 0
    cumdata = cumsum4d(cumsum4d(cumsum4d(cumsum4d(data_array, 1), 2), 3), 4)
	cum_cond = cumsum4di(cumsum4di(cumsum4di(cumsum4di(condition, 1), 2), 3), 4)
    histsize(1) = (datasize(1) - mod(datasize(1), nw(1))) / nw(1)
    histsize(2) = (datasize(2) - mod(datasize(2), nw(2))) / nw(2)
    histsize(3) = (datasize(3) - mod(datasize(3), nw(3))) / nw(3)
    histsize(4) = (datasize(4) - mod(datasize(4), nw(4))) / nw(4)
    hist_array = hist4d(cumdata, histsize)
    hist_cond = hist4di(cum_cond, histsize)
	maxcond = maxval(hist_cond)
	sumdata = sum(pack(hist_array, hist_cond == maxcond))
	sumhist = sum(pack(hist_cond, hist_cond == maxcond)/maxcond)
	do hx = 1, histsize(1)
	do hy = 1, histsize(2)
	do hz = 1, histsize(3)
	do hw = 1, histsize(4)
	  if (hist_cond(hx, hy, hz, hw) == maxcond) then
	    sumpdfs = sumpdfs + sum(pdfs((hx-1)*nw(1)+1:hx*nw(1),(hy-1)*nw(2)+1:hy*nw(2),(hz-1)*nw(3)+1:hz*nw(3),(hw-1)*nw(4)+1:hw*nw(4)))
	  endif
	enddo
	enddo
	enddo
	enddo
	do hx = 1, histsize(1)
	do hy = 1, histsize(2)
	do hz = 1, histsize(3)
	do hw = 1, histsize(4)
	  if (hist_cond(hx, hy, hz, hw) == maxcond) then
		do px = (hx-1)*nw(1)+1, hx*nw(1) 
		do py = (hy-1)*nw(2)+1, hy*nw(2)
		do pz = (hz-1)*nw(3)+1, hz*nw(3)
		do pw = (hw-1)*nw(4)+1, hw*nw(4)
		  if (sumdata > 0.0) then
		  ise(nw(1), nw(2), nw(3), nw(4)) = ise(nw(1), nw(2), nw(3), nw(4)) + (pdfs(px, py, pz, pw) -  &
			                        hist_array(hx, hy, hz, hw)/(sumdata/sumpdfs*sumpdf*product(nw)))**2
		  else
		  ise(nw(1), nw(2), nw(3), nw(4)) = 1.0
		  endif
	    enddo
	    enddo
	    enddo
	    enddo
      endif
	enddo
	enddo
	enddo
	enddo
	ise(nw(1), nw(2), nw(3), nw(4)) = ise(nw(1), nw(2), nw(3), nw(4)) / (product(nw)*sumhist) 
	print *, idx1d, rank, ise(nw(1), nw(2), nw(3), nw(4)), nw, sumhist, sumdata, sumpdfs, sumpdf

	enddo
    endif
	!enddo
	!enddo
	!enddo


	end subroutine cost4d


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

  function hist4d(cumdata, N)
    double precision, intent(in) :: cumdata(:,:,:,:) 
    integer, intent(in) :: N(4)
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

  function hist4di(cumdata, N)
    integer, intent(in) :: cumdata(:,:,:,:) 
    integer, intent(in) :: N(4)
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


end module costfort4d
