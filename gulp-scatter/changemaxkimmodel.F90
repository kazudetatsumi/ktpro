  subroutine changemaxkimmodel
!
!  Alters the size of the arrays associated with maxkimmodel
!
!  10/12 Created
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
!  Copyright Curtin University 2012
!
!  Julian Gale, NRI, Curtin University, October 2012
!
  use configurations, only : maxcfg
  use current,        only : maxat
  use kim_models
  use reallocate
#ifdef KIM
  use reallocate_kim
#endif
  implicit none
!
  integer(i4)       :: ierror, i, j
  integer(i4), save :: oldmaxkimmodel = 0
!
#ifdef KIM
  call realloc_kim(pkim_model,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','pkim_model')
#endif
  call realloc_ch80(kim_model,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_model')
  call realloc(lkim_model_cfg_OK,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','lkim_model_cfg_OK')
  call realloc(kim_cutoff,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_cutoff')
  call realloc(kim_nbc,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_nbc')
#ifdef KIM
  call realloc_kim_len(kim_test_descriptor,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_test_descriptor')
#endif
!
!  Initialise new part of data array
!
  if (maxkimmodel.gt.oldmaxkimmodel) then
    do i = 1,maxcfg
      do j = oldmaxkimmodel+1,maxkimmodel
        lkim_model_cfg_OK(j,i) = .true.
        kim_nbc(j,i) = 0
#ifdef KIM
        pkim_model(j,i) = 0
#endif
      enddo
    enddo
    do i = oldmaxkimmodel+1,maxkimmodel
      kim_model(i) = ' '
      kim_cutoff(i) = 0.0_dp
    enddo
  endif
!
!  Save current value of maxkimmodel for next call
!
  oldmaxkimmodel = maxkimmodel
!
  return
  end
