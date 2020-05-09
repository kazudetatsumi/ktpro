integer, parameter :: maxw = 10
integer nw
double precision, dimension(maxw, maxw) ::  A
  A = 0.0
  nw = 11
  A(nw, 129) = 0.5 
  print *, A


stop
end
