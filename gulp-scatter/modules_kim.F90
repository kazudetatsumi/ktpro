!************************************
!  Module for GULP-KIM interaction  *
!************************************
!
!  10/12 Created from modules.f90
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
!  NB: In this module the precision of some quantities is not set using the GULP defined types
!      since they must remain compatable with OpenKIM
!

!
!  OpenKIM models
!
  module kim_models
    use datatypes
#ifdef KIM
    use KIM_API
#endif
    integer(i4),                                parameter     :: kim_len = 10000                ! Length of KIM descriptor string
    integer(i4),                                         save :: maxkimmodel = 1                ! Maximum number of KIM models
    character(len=80),        dimension(:),     pointer, save :: kim_model => null()            ! Names of KIM models
#ifdef KIM
    integer(kind=kim_intptr), dimension(:,:),   pointer, save :: pkim_model => null()           ! Integer references for KIM pointers for each model/configuration
#endif
    integer,                  dimension(:,:),   pointer, save :: kim_nbc => null()              ! Array containing integer code NBC setting : 0 => cluster, 1=> RVEC_F
    integer,                  dimension(:),     pointer, save :: kim_ntype => null()            ! Array to hold number of particle types for KIM for each configuration
    integer,                  dimension(:,:),   pointer, save :: kim_types => null()            ! Array to hold particle types for KIM of each atom in each configuration
    logical,                  dimension(:,:),   pointer, save :: lkim_model_cfg_OK => null()    ! Stores logical check that configuration is OK for KIM model
    real*8,                   dimension(:,:),   pointer, save :: kim_coord => null()            ! Array to hold coordinates for KIM
    real*8,                   dimension(:),     pointer, save :: kim_cutoff => null()           ! Array to hold cutoffs for KIM models
    real*8,                   dimension(:,:),   pointer, save :: kim_forces => null()           ! Array to hold forces for KIM
!
    character(len=kim_len),   dimension(:,:),   pointer, save :: kim_test_descriptor => null()  ! KIM descriptor for test - used instead of a .kim file
    integer(i4),                                         save :: nkimmodel = 0                  ! Number of KIM models
    integer(i4),                                         save :: nlibnkimmodel = 0              ! Number of KIM models prior to libraries
    integer(i4),                                         save :: kim_types_dummy(1)             ! Dummy variable for KIM particle types
    logical,                                             save :: lkim_model = .false.           ! If true then OpenKIM model to be used
    real*8,                                              save :: kim_energy                     ! Variable to hold KIM energy
    real*8,                                              save :: kim_virial(6)                  ! Variable to hold KIM virial
  end module kim_models
