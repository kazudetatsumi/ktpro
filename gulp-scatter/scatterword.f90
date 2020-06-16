  subroutine scatterword(nru,word,lwordok,iline,line,l55,l1000,ncurr)
!
!  Processes input for scatter related words
!
!  nru = fortran channel for reading input
!
!  2/10 Created from 'phonword.f90'
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
!  Copyright Curtin University 2010
!
!  Julian Gale, NRI, Curtin University, September 2010
!
  use configurations
  use gulpinput
  use iochannels
  use parallel
  use scatterdata
  implicit none
!
!  Passed variables
!
  character(len=20)            :: word
  character(len=maxlinelength) :: line
  integer(i4)                  :: iline
  integer(i4)                  :: ncurr
  integer(i4)                  :: nru
  logical                      :: lwordok
  logical                      :: l55
  logical                      :: l1000
!
!  Local variables
!
  integer(i4)                  :: i
  integer(i4)                  :: j
  real(dp)                     :: ords(6)
!
!  Input for scatter
!

  if (index(word,'k_ty').eq.1) goto 740
  if (index(word,'s_op').eq.1) goto 750
  if (index(word,'s_ty').eq.1) goto 760
  if (index(word,'s_an').eq.1) goto 770
  if (index(word,'m_tr').eq.1) goto 780
  if (index(word,'e_tr').eq.1) goto 790
  if (index(word,'sc_f').eq.1) goto 800
  if (index(word,'s_ou').eq.1) goto 820
  return

!*****************************************************************************************
!  k_type:  Kernel Selection                                                             *
!  Select the kind of S(Q,w) calculation to be included with the data.                   *
!                                                                                        *
!  coherent:    Calculates the coherent inelastic contribution to S(Q,w) only (default). *
!  incoherent:  Calculates the incoherent inelastic contribution to S(Q,w) only.         *
!  t_scatter:   Calculates the coherent & incoherent inelastic contributions to S(Q,w).  *
!  com_scatter: Calculates all contributions to S(Q,w) - complete scattering (sets flags)*
!  k_empty:     Runs the sampling process, but doesn't calculate S(Q,w).                 *
!*****************************************************************************************
740 continue
  do i = 1,nword
  if (index(words(i),'coh').ne.0) then
    q_caltype = 1
    elseif (index(words(i),'inc').ne.0) then
      q_caltype = 2
    elseif (index(words(i),'t_sc').ne.0) then
      q_caltype = 3
    elseif (index(words(i),'com_s').ne.0) then
      q_caltype = 4
    elseif (index(words(i),'k_em').ne.0) then
      q_caltype = 5
    endif
  enddo
  lwordok = .true.
  return

!*****************************************************************************************
!  s_options:  Scattering options                                                        *
!  Select from the various additional options available.                                 *
!                                                                                        *
!  d_dwaller:   Dynamic debye-waller factor calculation (default).                       *
!  l_dwaller:   Library debye-waller factor, taken from a prior calculation (file).      *
!  u_dwaller:   User-defined debye-waller factor, input in option file. *                *
!  m_phonon:    Multi-phonon scattering correction to S(Q,w).                            *
!  m_scatter:   Multiple scattering correction to S(Q,w).                                *
!  f_lsquare:   Fitting force constants using least squares.                             *
!  f_geofit:    Fitting force constants using the GEO algorithm.                         *
!  f_gimethod:  Fitting force constants using GSM method.                                *
!                                                                                        *
!              * denotes the requirement for a sc_filename.sopt option file              *
!                to be parsed along with the .gin file.                                  *
!*****************************************************************************************
750 continue
  do i = 1,nword
    if (index(words(i),'d_dw').ne.0) then
      lscatopt(1)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'l_dw').ne.0) then
      lscatopt(2)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'u_dw').ne.0) then
      lscatopt(3)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'m_ph').ne.0) then
      lscatopt(4)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'m_sc').ne.0) then
      lscatopt(5)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'f_ls').ne.0) then
      lscatopt(6)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'f_ge').ne.0) then
      lscatopt(7)=.true.
      lscatoptflag = .true.
    endif
    if (index(words(i),'f_gi').ne.0) then
      lscatopt(8)=.true.
      lscatoptflag = .true.
    endif
  enddo
  lwordok = .true.
  return

!*****************************************************************************************
!  s_type -    Sampling Type                                                             *
!  Select the kind of reciprocal space sampling algorithm to use.                        *
!                                                                                        *
!  rs_onion:    Reciprocal space onion method for sampling reciprocal space (default).   *
!  ls_set:      Linear samplingb set method for sampling reciprocal space.               *
!  rs_cylinder: Cylindrical sampling method for sampling reciprocal space.               *
!  mp_sample:   Monkhorse-pack method for sampling reciprocal space.                     *
!  mc_sample:   Monte Carlo method for sampling reciprocal space.                        *
!  sc_sample:   Single crystal sampling, using gamma point to <a*b*c*> as a single vector*
!  pe_file:     Parse existing file.  Code reads in (Q,w) space sample from file         *
!*****************************************************************************************
760 continue
  do i = 1,nword
    if (index(words(i),'rs_o').ne.0) then
      q_runtype = 1
      out_eig_flag = .true.
    elseif (index(words(i),'ls_s').ne.0) then
      q_runtype = 2
      out_eig_flag = .true.
    elseif (index(words(i),'rs_c').ne.0) then
      q_runtype = 3
      out_eig_flag = .true.
    elseif (index(words(i),'mp_s').ne.0) then
      q_runtype = 4
      out_eig_flag = .true.
    elseif (index(words(i),'mc_s').ne.0) then
      q_runtype = 5
      out_eig_flag = .true.
    elseif (index(words(i),'sc_s').ne.0) then
      q_runtype = 6
      out_eig_flag = .true.
    elseif (index(words(i),'pe_f').ne.0) then
      q_runtype = 7
    endif
  enddo
  lwordok = .true.
  return

!*****************************************************************************************
!  s_analysis - Pre-analysis options for output                                          *
!  Select the type of pre-analysis options to use, all output as additional files.       *
!                                                                                        *
!  i_slice:  Takes a set of slices through the data and saves as qfilename(N).slc        *
!  i_mlocus: Takes a set of fixed-angle loci through the data and saves as q~(N).lci     *
!  o_tc:                                                                                 *
!  o_mp:     Outputs contribution from multi-phonon terms to qfilename.mpc               *
!  o_ms.     Outputs contribution from multiple scattering terms to qfilename.msc        *
!  o_dw:     Outputs contribution from debye-waller convolution terms to qfilename.dwf   *
!  o_ge:                                                                                 *
!  o_o1:     Option 1 (unused - for later expansion)                                     *
!  o_o2:     Option 1 (unused - for later expansion)                                     *
!                                                                                        *
!*****************************************************************************************
770 continue
  do i = 1,nword
    if (index(words(i),'i_sl').ne.0) then
      lscatanopt(1) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'i_ml').ne.0) then
      lscatanopt(2) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'o_tc').ne.0) then
      lscatanopt(3) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'0_mp').ne.0) then
      lscatanopt(4) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'o_ms').ne.0) then
      lscatanopt(5) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'o_dw').ne.0) then
      lscatanopt(6) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'o_ge').ne.0) then
      lscatanopt(7) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'o_o1').ne.0) then
      lscatanopt(8) = .true.
      lscatanoptflag = .true.
    endif
    if (index(words(i),'o_o2').ne.0) then
      lscatanopt(9) = .true.
      lscatanoptflag = .true.
    endif
  enddo
  lwordok = .true.
  return

!*****************************************************************************************
!  m_transfer -Momentum transfer input                                                   *
!  Options for the range of Q, increment size for various sampling options               *
!                                                                                        *
!  rs_onion  [  Minimum Q, maximum Q, Number of layers, number of angular steps   ]      *
!  ls_set    [  Number of sampling vectors,                                              *
!               vector_origin1<a* b* c*>, vector_ordinate1<a* b* c*>,number of steps     *
!               ...                                                                      *
!               vector_originN<a* b* c*>, vector_ordinateN<a* b* c*>,number of steps ]   *
!  rs_cylinder [Minimum Q (axis), maximum Q(axis), basal radius, number of layers,       *
!               number of angular steps ]                                                *
!  mp_sample  [ Minimum Q, maximum Q, sampling shrinking factor in reciprocal lattice    *
!               na*, nb*, nc*, weighting factors for each element wa*, wb*, wc*  ]       *
!  mc_sample  [ Minimum Q, maximum Q, isotropic sampling factor  ]                       *
!  sc_sample  [ Minimum Q, maximum Q, vector_ordinal<a* b* c*>   ]                       *
!*****************************************************************************************
780 if (q_runtype.eq.1) then
!
!  Input for rs_onion sampling method
!
    out_eig_flag = .true.
    if (nfloat.ge.4) then
      q_qmin     = floats(1)
      q_qmax     = floats(2)
      nq_step    = nint(floats(3))
      nq_intstep = nint(floats(4))
    elseif (nfloat.eq.3) then
      q_qmin     = 0.0_dp
      q_qmax     = floats(1)
      nq_step    = nint(floats(2))
      nq_intstep = nint(floats(3))
    elseif (nfloat.eq.2) then
      q_qmin     = 0.0_dp
      q_qmax     = floats(1)
      nq_step    = 100
      nq_intstep = nint(floats(2))
    elseif (nfloat.eq.1) then
      q_qmin     = 0.0_dp
      q_qmax     = floats(1)
      nq_step    = 100
      nq_intstep = 100
    else
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nfloat.lt.4) then
        call outerror('Insufficient input',iline)
        call stopnow('scatterword')
      endif
      q_qmin = floats(1)
      q_qmax = floats(2)
      nq_step = nint(floats(3))
      nq_intstep = nint(floats(4))
    endif
  elseif (q_runtype.eq.2) then
!
!  Input Q for linear sampling method 
!
    out_eig_flag = .true.
    if (nfloat.gt.0) then
      n_qvector = nint(floats(1))
    else
      n_qvector = 1
    endif
    if (n_qvector.gt.maxqvector) then
      maxqvector = n_qvector
      call changemaxqvector
    endif
    do i = 1,n_qvector
      line = '  '
      read(nru,'(a)') line
      iline = iline + 1
      call linepro(nru,line,iline)
      if (nfloat.lt.7) then
        call outerror('Each lss input requires 6 vector co-ords + number of qpoints per vector.',iline)
        call stopnow('scatterword')
      endif
      do j = 1,6
        q_ordinate(j,i) = floats(j)
      enddo
      n_qpoint(i) = nint(floats(7))
      if (n_qpoint(i).gt.n_qpointmax) n_qpointmax = n_qpoint(i)
    enddo 
  elseif (q_runtype.eq.3) then
!
!  Rs_cylinder sampling method
!
  elseif (q_runtype.eq.4) then
!
  elseif (q_runtype.eq.5) then
!
  elseif (q_runtype.eq.6) then
!
!  Single vector from the gamma point sampling method
!
    out_eig_flag = .true.
    n_qvector = 1
    if (n_qvector.gt.maxqvector) then
      maxqvector = n_qvector
      call changemaxqvector
    endif
    if (nfloat.ge.4) then
      ords(1) = 0.0_dp
      ords(2) = 0.0_dp
      ords(3) = 0.0_dp
      ords(4) = floats(1)
      ords(5) = floats(2)
      ords(6) = floats(3)
      n_qpoint(1) = nint(floats(4))
      do j = 1,6
        q_ordinate(j,1) = ords(j)
      enddo
    elseif (nfloat.eq.3) then
      ords(1) = 0.0_dp
      ords(2) = 0.0_dp
      ords(3) = 0.0_dp
      ords(4) = floats(1)
      ords(5) = floats(2)
      ords(6) = floats(3)
      n_qpoint(1) = 100_i4
    else
      call outerror('Must specify three co-ordinates for Q direction plus (optional) increments',iline)
      call stopnow('scatterword')
    endif
    n_qpointmax = n_qpoint(1)
  endif
!    
  if (q_runtype.eq.7) then
!
!  pe_file opens the file q_pe_filename.extension and reads contents into memory, depending on format.
!
  endif
  lwordok = .true.
  return
    
!*****************************************************************************************
!  e_transfer - Energy transfer input                                                    *
!  generic across the set of Q sampling                                                  *
!                                                                                        *  
!  e_transfer [ minimum energy, maximum energy, number of steps                          *
!               units {cm-1, THz, meV (default)} ]                                       *
!                                                                                        *
!*****************************************************************************************
790 if (nfloat.ge.3) then
    q_wmin  = floats(1)
    q_wmax  = floats(2)
    nw_step = nint(floats(3))
  elseif (nfloat.eq.2) then
    q_wmin  = 0.0_dp
    q_wmax  = floats(1)
    nw_step = nint(floats(2))
  elseif (nfloat.eq.1) then
    q_wmin  = 0.0_dp
    q_wmax  = nint(floats(1))
    nw_step = 100
  elseif (nfloat.lt.1) then
    call outerror('Must specify min & max freqs and no. of steps',iline)
    call stopnow('scatterword')
  endif
  lwordok = .true.
  return
     
!*****************************************************************************************
!  sc_filename - Name of the scattering file to be used for all output headers.          *
!*****************************************************************************************
800 if (nword.ge.2) then
    q_filename = words(2)
  else
    line = ' '
    read(nru,'(a)') line
    iline = iline + 1
    call linepro(nru,line,iline)
    if (nword.ge.1) then
      q_filename = words(1)
    endif
  endif
  lwordok = .true.
  return

!*****************************************************************************************
!  s_output -     Suppression/engagement of output options                               *
!*****************************************************************************************
820 continue
  do i = 1,nword
    if (index(words(i),'n_ei').ne.0) then
      out_eig_flag = .false.
    endif
  enddo
  lwordok = .true.
  return
!
  end
