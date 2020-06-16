  subroutine outshengBTEcontrol(iout)
!
!  Outputs file CONTROL for ShengBTE.
!
!   3/14 Created
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
!  Copyright Curtin University 2014
!
!  Julian Gale, NRI, Curtin University, March 2014
!
  use constants
  use control
  use current
  use files
  use iochannels
  use ksample
  use m_pdfneutron, only : lshrinkset
  use properties,   only : dicons
  use species,      only : nspec, symspec
  implicit none
!
!  Passed variables
!
  integer(i4),                    intent(in)   :: iout
!
!  Local variables
!
  character(len=5)                             :: lab
  integer(i4)                                  :: i
  integer(i4), dimension(:), allocatable       :: itmp
  integer(i4)                                  :: itype
  integer(i4)                                  :: j
  integer(i4)                                  :: status
  logical                                      :: lnboxzero
  real(dp)                                     :: fmaxdos
!
!  Open file
!
  open(iout,file='CONTROL',status='unknown',form='formatted')
!-------------------------------------------------------------------------------
!  Allocations block
!-------------------------------------------------------------------------------
  write(iout,'(''&allocations'')')
  write(iout,'(8x,''nelements='',i6)',advance='no') nspec
  write(iout,'('','')')
  write(iout,'(8x,''natoms='',i6)',advance='no') numat
  write(iout,'('','')')
  if (lshrinkset(ncf)) then
    write(iout,'(8x,''ngrid(:)='',3i4)') nxks(ncf),nyks(ncf),nzks(ncf)
  else
    write(iout,'(8x,''ngrid(:)=1 1 1'')') 
  endif
  write(iout,'(''&end'')')
!-------------------------------------------------------------------------------
!  Crystal block
!-------------------------------------------------------------------------------
  write(iout,'(''&crystal'')')
!
!  Output lattice factor as conversion from Angstroms to nm
!
  write(iout,'(8x,''lfactor=0.1'')',advance='no') 
  write(iout,'('','')')
!
!  Lattice vectors
!
  write(iout,'(8x,''lattvec(:,1)='',3(f12.6,1x))',advance='no') (rv(j,1),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''lattvec(:,2)='',3(f12.6,1x))',advance='no') (rv(j,2),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''lattvec(:,3)='',3(f12.6,1x))',advance='no') (rv(j,3),j=1,3)
  write(iout,'('','')')
!
!  Elements
!
  write(iout,'(8x,''elements='')',advance='no') 
  do i = 1,numat
    write(iout,'(''"'')',advance='no') 
    call label(nat(i),0_i4,lab)
    write(iout,'(a)',advance='no') trim(lab)
    if (i.eq.numat) then
      write(iout,'(''"'')') 
    else
      write(iout,'(''" '')',advance='no') 
    endif
  enddo
  write(iout,'(8x,''types='')',advance='no')
  do i = 1,numat
    write(iout,'(i3,1x))',advance='no') nspecptr(nrelat(i))
  enddo
  write(iout,'('','')')
!
!  Atomic fractional coordinates
!
  do i = 1,numat
    write(iout,'(8x,''positions(:,'',i3,'')='',3(f12.6,1x))',advance='no') i,xfrac(i),yfrac(i),zfrac(i)
    write(iout,'('','')')
  enddo
!
!  Static dielectric constant tensor
!
  write(iout,'(8x,''epsilon(:,1)='',3(f12.6,1x))',advance='no') (dicons(j,1),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''epsilon(:,2)='',3(f12.6,1x))',advance='no') (dicons(j,2),j=1,3)
  write(iout,'('','')')
  write(iout,'(8x,''epsilon(:,3)='',3(f12.6,1x))',advance='no') (dicons(j,3),j=1,3)
  write(iout,'('','')')
!
!  Born effective charges
!
  do i = 1,numat
    write(iout,'(8x,''born(:,1,'',i3,'')='',3(f12.6,1x))',advance='no') i,(bornq(j,1,i),j=1,3)
    write(iout,'('','')')
    write(iout,'(8x,''born(:,2,'',i3,'')='',3(f12.6,1x))',advance='no') i,(bornq(j,2,i),j=1,3)
    write(iout,'('','')')
    write(iout,'(8x,''born(:,3,'',i3,'')='',3(f12.6,1x))',advance='no') i,(bornq(j,3,i),j=1,3)
    write(iout,'('','')')
  enddo
  write(iout,'(8x,''scell(:)=1 1 1'')') 

  write(iout,'(''&end'')')
!-------------------------------------------------------------------------------
!  Parameters block
!-------------------------------------------------------------------------------
  write(iout,'(''&parameters'')')

  write(iout,'(''&end'')')
!-------------------------------------------------------------------------------
!  Flags block
!-------------------------------------------------------------------------------
  write(iout,'(''&flags'')')

  write(iout,'(''&end'')')
!
!  Close file
!
  close(iout)

  return
  end
