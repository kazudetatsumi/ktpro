  subroutine outderv
!
!  Outputs second derivatives for debugging
!
!   5/02 Created from strfin
!   5/09 Output banners modified and space added after strain second derivs
!   7/13 Modified to allow for the fact that the second derivative matrix
!        may not have been symmetrised prior to this call
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
!  Julian Gale, NRI, Curtin University, July 2013
!
  use control
  use current
  use derivatives
  use iochannels
  use symmetry
  implicit none
!
!  Local variables
!
  integer(i4) :: i
  integer(i4) :: j
  integer(i4) :: k
  integer(i4) :: maxlim
  integer(i4) :: matom
!
  if (lsymderv2) then
    matom = nasym
  else
    matom = numat
  endif
  maxlim = 3*matom
  if (nbsm.gt.0) maxlim = maxlim + matom
!
!  Write out a space of any of the options are to be output
!
  if (index(keyword,'derv2').ne.0.or.index(keyword,'derv3').ne.0) then
    write(ioout,'(/)')
  endif
!********************
!  Strain - strain  *
!********************
  if (index(keyword,'derv2').ne.0.and.lstr) then
    write(ioout,'(''  Strain-strain second derivative matrix : (eV)'',/)')
    do i = 1,nstrains
      write(ioout,'(6f12.5)')(sderv2(j,i),j=1,nstrains)
    enddo
    write(ioout,'(/)')
  endif
!**********************
!  Internal - strain  *
!**********************
  if (index(keyword,'derv3').ne.0.and.lstr) then
    write(ioout,'(''  Mixed strain-internal second derivative matrix : '',''(eV/Angstrom)'',/)')
    do i = 1,maxlim
      write(ioout,'(9f12.5)')(derv3(i,j),j=1,nstrains)
    enddo
    write(ioout,'(/)')
  endif
!************************
!  Internal - internal  *
!************************
  if (index(keyword,'derv2').ne.0) then
    write(ioout,'(''  Internal-internal second derivative matrix : '',''(eV/Angstrom**2)'',/)')
    do i = 1,maxlim
      write(ioout,'(9f12.5)')(derv2(j,i),j=1,i),(derv2(i,k),k=i+1,maxlim)
    enddo
    write(ioout,'(/)')
  endif
!
  return
  end
