subroutine cost1d(maxw, lenA, A, Cn, kaves, deltas) bind(C)
    integer, intent(in) :: maxw
    integer, intent(in) :: lenA
    double precision, intent(in), dimension(lenA) :: A                        
    double precision, intent(inout), dimension(maxw) :: Cn, kaves, deltas                        
    double precision, dimension(lenA) :: k                        
    integer nw0, N0, i, ihead
    real kave, v
    k = 0.0
    Cn = 0.0
    kaves = 0.0
    deltas = 0.0
    do nw0 = 1, maxw
       N0 = (lenA - mod(lenA, nw0)) / nw0
       do i = 1, N0
          ihead = i*nw0 
          if ( i == 1) then
              k(i) = A(ihead)
          else
              k(i) = A(ihead) - A(ihead - nw0)
          end if
       end do
       kave = sum(k(1:N0)) / real(N0)
       v = sum((k(1:N0) - kave)**2) / real(N0)
       Cn(nw0) = (2.0 * kave - v) / (real(nw0)**2)
       deltas(nw0) = real(nw0)
       kaves(nw0) = kave
    end do
    print *, 'size of Cn',size(Cn)
    print *, 'len0', nw0
    print *, 'maxw', maxw
end subroutine cost1d
