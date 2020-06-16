  subroutine harmonicrelax(nvari,nmin,xc,fc,gc,hessian,maxhess,imode)
!
!  Estimate the energy change on relaxation using the harmonic approximation
!
!  imode   = 1 => bulk calculation
!  imode   = 2 => defect calculation
!
!   3/14 Created from minimise
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
!  Julian Gale, NRI, Curtin University, March 2014
!
  use control
  use current
  use defects
  use derivatives, only : derv2
  use energies,    only : erelax, efreeze, fcsave
  use iochannels
  use optimisation
  use parallel

  implicit none
!
!  Passed variables
!
  integer(i4),          intent(in)            :: imode
  integer(i4),          intent(in)            :: nmin
  integer(i4),          intent(in)            :: nvari
  integer(i4),          intent(in)            :: maxhess
  real(dp),             intent(inout)         :: fc
  real(dp),             intent(inout)         :: gc(nvari)
  real(dp),             intent(inout)         :: hessian(maxhess)
  real(dp),             intent(inout)         :: xc(nvari)
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: ind
  integer(i4)                                 :: info
  integer(i4)                                 :: iflag
  integer(i4)                                 :: j
  integer(i4)                                 :: nlvar
  integer(i4)                                 :: status
  integer(i4),              allocatable, save :: kpvt(:)
  real(dp)                                    :: ddot
  real(dp)                                    :: fcin
  real(dp)                                    :: funct1
  real(dp),                 allocatable, save :: pvect(:)
  real(dp),                 allocatable, save :: temp(:)
!
!  Set local number of variables
!
  nlvar = nvari - nmin + 1
!
!  If number of variables is zero then relaxation energy is zero
!
  if (nlvar.eq.0) then
    erelax = 0.0_dp
    return
  endif
!
!  Allocate local memory
!
  allocate(pvect(nvari),stat=status)
  if (status/=0) call outofmemory('minimise','pvect')
  allocate(kpvt(nlvar),stat=status)
  if (status/=0) call outofmemory('nrhess','kpvt')
  allocate(temp(3*nlvar),stat=status)
  if (status/=0) call outofmemory('nrhess','temp')
!
!  Calculate initial function, gradients and necessary second derivatives
!
  iflag = 2
  fcin = fc
  if (imode.eq.1) then
    call funct(iflag,nvari,xc,fc,gc)
  else
    call deffun(iflag,nvari,xc,fc,gc)
  endif
  fcsave = fc
  if (lfreeze) then
    efreeze = fcin - fc
    fc = fc + efreeze
  else
    efreeze = 0.0_dp
  endif
!***************************
!  Generate exact hessian  *
!***************************
  iflag = 2
  if (imode.eq.1) then
    call funct(iflag,nvari,xc,funct1,gc)
  else
    call deffun(iflag,nvari,xc,funct1,gc)
  endif
!
!  Set Hessian
!
  if (imode.eq.1) then
    if (ndim.ge.1) then
      if (lfreeze) then
        call sec3f
      else
        call sec3
      endif
    else
      if (lfreeze) then
        call sec0f
      else
        call sec0
      endif
    endif
  else
    call defsec
  endif

!  Pack hessian and save diagonal elements
!
  ind = 0
  do i = nmin,nvar
    do j = nmin,i
      ind = ind + 1
      hessian(ind) = derv2(j,i)
    enddo
  enddo
!
!  Factorise matrix
!
  call dsptrf('U',nlvar,hessian,kpvt,info)
!
!  Check for singularities
!
  if (info.gt.0) then
    call outerror('hessian inversion failed during estimation of relaxation energy',0_i4)
    call stopnow('harmonicrelax')
  else
!
!  Complete inversion
!
    call dsptri('U',nlvar,hessian,kpvt,temp,info)
  endif
!
!  Multiply inverse Hessian by gradient to get displacement vector
!
  call dspmv('U',nlvar,1.0_dp,hessian,gc(nmin),1_i4,0.0_dp,pvect(nmin),1_i4)
!
!  Compute predicted energy charge within the harmonic approximation
!
  erelax = - ddot(nlvar,pvect(nmin),1_i4,gc(nmin),1_i4)
  erelax = 0.5_dp*erelax
!
!  Free local memory
!
  deallocate(temp,stat=status)
  if (status/=0) call deallocate_error('nrhess','temp')
  deallocate(kpvt,stat=status)
  if (status/=0) call deallocate_error('nrhess','kpvt')
  deallocate(pvect,stat=status)
  if (status/=0) call deallocate_error('minimise','pvect')
!
  return
  end
