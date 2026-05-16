  subroutine potwordk(iin,word,lwordok,iline,line,l55,l1000,linr,lint,llibrary,lfflags,ncurr)
!
!  Processes potential input for OpenKIM potential related information.
!
!  iin = input fortran channel
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
  use element,    only : maxele
  use gulpinput
  use kim_models, only : kim_model, lkim_model, nkimmodel, maxkimmodel
  use parallel
  implicit none
!
!  Passed variables
!
  character(len=20)            :: word
  character(len=maxlinelength) :: line
  integer(i4)                  :: iin
  integer(i4)                  :: iline
  integer(i4)                  :: ncurr
  logical                      :: l55
  logical                      :: l1000
  logical                      :: lfflags
  logical                      :: linr
  logical                      :: lint
  logical                      :: llibrary
  logical                      :: lwordok
!
!  Local variables
!

!
!  Initialise local variables
!
  if (index(word,'kim_m').eq.1) goto 100
  return
!******************
!  OpenKIM model  *
!******************
100 continue
!
!  Read next line
!
  line = ' '
  read(iin,'(a)',end=108) line
  iline = iline + 1
!
!  Increment number of KIM models
!
  nkimmodel = nkimmodel + 1
!
!  Check against value of maxkimmodel
!
  if (nkimmodel.gt.maxkimmodel) then
    maxkimmodel = nkimmodel
    call changemaxkimmodel
  endif
!
!  Assign string to KIM model name
!
  kim_model(nkimmodel) = line
!
!  Set flag to true that indicates we have an OpenKIM model
!
  lkim_model = .true.
!
  lwordok = .true.
  return
!
!  Escape point for end of file
!
108 call outerror('End of input when expecting OpenKIM model name',iline)
  call stopnow('potwordk')
!
  end
