  subroutine ramanstrength(mcv,nphonatc,nphonatptr,ncfoc,iocptr,eigr,maxd2,ramstrength)
!
!  Calculates the Raman strengths for the modes. Gamma point version.
!
!  10/13 Created by separating from oscillatorstrengthg
!
!  On entry :
!
!  mcv         = no. of modes
!  nphonatc    = total number of cores in relevant regions
!  nphonatptr  = pointer from reduce to full atom set
!  ncfoc       = number of condensed core sites
!  iocptr      = pointer that connects sites to condensed sites
!  eigr        = eigenvectors of dynamical matrix
!  maxd2       = left-hand dimension of eigr
!
!  On exit : 
!
!  ramstrength = Raman susceptibility tensors, if lraman = .true.
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2013
!
!  Julian Gale, NRI, Curtin University, October 2013
!
  use control,    only : lraman
  use current
  use element
  use iochannels
  use properties, only : ramanasus
  use species,    only : massspec
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: iocptr(*)
  integer(i4), intent(in)    :: maxd2
  integer(i4), intent(in)    :: mcv
  integer(i4), intent(in)    :: nphonatc
  integer(i4), intent(in)    :: nphonatptr(*)
  integer(i4), intent(in)    :: ncfoc
  real(dp),    intent(in)    :: eigr(maxd2,mcv)
  real(dp),    intent(out)   :: ramstrength(3,3,mcv)
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ix
  integer(i4)                :: iy
  integer(i4)                :: iz
  integer(i4)                :: j
  integer(i4)                :: m
!
  if (lraman) then
!**********************************
!  Oscillator strengths per mode  *
!**********************************
!
!  Loop over modes
!
    do m = 1,mcv
      ramstrength(1:3,1:3,m) = 0.0_dp
!
!  Loop over full sites
!
      do i = 1,ncfoc
        ix = 3*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
!
!  Find all cores associated with full site
!
        do j = 1,nphonatc
          if (iocptr(j).eq.i) then
!
!  Multiple inverse mass weighted eigenvectors by Born charges
!
            ramstrength(1:3,1:3,m) = ramstrength(1:3,1:3,m) + ramanasus(1:3,1:3,1,j)*eigr(ix,m) &
                                                            + ramanasus(1:3,1:3,2,j)*eigr(iy,m) &
                                                            + ramanasus(1:3,1:3,3,j)*eigr(iz,m) 
          endif
        enddo
      enddo
    enddo
  else
!
!  Loop over modes
!
    do m = 1,mcv
      ramstrength(1:3,1:3,m) = 0.0_dp
    enddo
  endif
!
  return
  end
