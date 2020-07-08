  subroutine scatter(lprint,lfitrun,fsumsq)
!
!  Calculates the phonons at a given set of k points determined by
!  the single phonon scattering condition and then calculates the
!  coherent and incoherent inelastic cross-sections along fixed Q 
!  or for a range of Q values determined by sampling method.
!
!  Additional functionality has been added to add corrections for multi-phonon and
!  multiple scattering, a range of sampling methods and tools for S(Q,w) library
!  building for models, pre-analysis tools from the 'Scattools' program and I/O
!  options for further analysis.
!  
!  The keyword would be: scatter
!  and the <options>, with [switches] are:
!
!  <k_type>:  Kernel Selection                                                             
!  Select the kind of S(Q,w) calculation to be included with the data.                   
!                                                                                        
!  [coherent]:    Calculates the coherent inelastic contribution to S(Q,w) only (default). 
!  [incoherent]:  Calculates the incoherent inelastic contribution to S(Q,w) only.         
!  [t_scatter]:   Calculates the coherent & incoherent inelastic contributions to S(Q,w).  
!  [com_scatter]: Calculates all contributions to S(Q,w) - complete scattering (sets flags)
!  [k_empty]:     Runs the sampling process, but doesn't calculate S(Q,w).        
!
!  <s_options>:  Scattering options                                                        
!  Select from the various additional options available.                                 
!                                                                                        
!  [d_dwaller]:   Dynamic debye-waller factor calculation (default).                       
!  [l_dwaller]:   Library debye-waller factor, taken from a prior calculation (file).      
!  [u_dwaller]:   User-defined debye-waller factor, input in option file.                  
!  [m_phonon]:    Multi-phonon scattering correction to S(Q,w).                            
!  [m_scatter]:   Multiple scattering correction to S(Q,w).                                
!  [f_lsquare]:   Fitting force constants using least squares.                             
!  [f_geofit]:    Fitting force constants using the GEO algorithm.                         
!  [f_gimethod]:  Fitting force constants using GSM method.                                
!                                                                                        
!              * denotes the requirement for a sc_filename.sopt option file              
!                to be parsed along with the .gin file.                                  
!
!  <s_type> -    Sampling Type                                                             
!  Select the kind of reciprocal space sampling algorithm to use.                                                                           *
!                                                                                        
!  [rs_onion]:    Reciprocal space onion method for sampling reciprocal space (default).   
!  [ls_set]:      Linear samplingb set method for sampling reciprocal space.               
!  [rs_cylinder]: Cylindrical sampling method for sampling reciprocal space.               
!  [mp_sample]:   Monkhorse-pack method for sampling reciprocal space.                     
!  [mc_sample]:   Monte Carlo method for sampling reciprocal space.                        
!  [sc_sample]:   Single crystal sampling, using gamma point to <a*b*c*> as a single vector
!  [pe_file]:     Parse existing file.  Code reads in (Q,w) space sample from file         
!
!  <s_analysis> - Pre-analysis options for output                                          
!                                                                                        
!  [i_slice]:  Takes a set of slices through the data and saves as qfilename(N).slc        
!  [i_mlocus]: Takes a set of fixed-angle loci through the data and saves as q~(N).lci     
!  [o_tc]:     Threshold comparison function                                                                            
!  [o_mp]:     Subtract multi-phonon contribution from main scan.  Two files are produced: a 
!              qfilename_mpc.sqw and a qfilename.mpc file containing the contribution itself.
!  [o_ms].     Subtract multiple scattering contribution from main scan.  Two files are produced: a 
!              qfilename_msc.sqw and a qfilename.msc file containing the contribution itself.                                                                            
!  [o_dw]:     Subtract debye-waller contribution from main scan.  Two files are produced: a 
!              qfilename_dwc.sqw and a qfilename.dwc file containing the contribution itself.                                   
!  [o_ge]:     Use GEO matching, according to the GEO grid used.  One file is produced with the grid
!              info in a table format, and then all the grids afterward.
! 
!  m_transfer -Momentum transfer input                                                   
!  Options for the range of Q, increment size for various sampling option switches, given as:              
!                                                                                        
!  rs_onion  [  Minimum Q, maximum Q, Number of layers, number of angular steps   ]      
!  ls_set    [  Number of sampling vectors, number of steps along vector                 
!               vector_origin1<a* b* c*>, vector_ordinate1<a* b* c*>,                    
!               ...                                                                      
!               vector_originN<a* b* c*>, vector_ordinateN<a* b* c*> ]                              
!  rs_cylinder [Minimum Q (axis), maximum Q(axis), basal radius, number of layers,       
!               number of angular steps ]                                                
!  mp_sample  [ Minimum Q, maximum Q, sampling shrinking factor in reciprocal lattice    
!               na*, nb*, nc*, weighting factors for each element wa*, wb*, wc*  ]       
!  mc_sample  [ Minimum Q, maximum Q, isotropic sampling factor  ]                       
!  sc_sample  [ Minimum Q, maximum Q, vector_ordinal<a* b* c*>   ]                       
!
!  e_transfer - Energy transfer input                                                    
!  generic across the set of Q sampling                                                  
!                                                                                          
!  e_transfer [ minimum energy, maximum energy, number of steps                          
!               units {cm-1, THz, meV (default)} ]                                       
!                                                                                        
!  sc_filename - Name of the scattering file to be used for all output headers.          
!
!  Input options:
!
!  lfitrun  - if true then this is a fitting call and fsumsq should be computed
!  fsumsq   - if lfitrun is true then this is the sum of the squares of the 
!             differences between endarray calculated and input.
!
!   9/10 Standardisation for GULP and merge with main trunk performed
!   9/10 Transpose of kv arrays used to correct results
!   9/10 Fitting of S(Q,omega) added
!   9/10 K point parallelisation added for calculation of S(Q,omega) added
!   9/10 Writing of files suppressed if scatter is called with lprint = .false.
!        This means it is a fitting run and so writing files costs CPU time and
!        is unnecessary.
!   1/12 Phonon calls modified with dummy parameters
!   6/20 cputime -> g_cpu_time renamed       KT@JPARC
!   6/20 matinv  -> matrix_inversion renamed KT@JPARC
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. !opies should be
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
!  Daniel L. Roach, University of Salford, 2010
!  Julian Gale,  Curtin University, 2012
!
  use constants
  use control
  use current
  use element
  use files
  use frequencies,     only : freq
  use general
  use iochannels
  use ksample_scatter, only : xskpt, yskpt, zskpt, wskpt, maxskpt, nskpt
  use parallel
  use partial
  use properties
  use scatterdata
  use shell
  use symmetry,        only : lsymderv2
  use times
  use xcgc
  implicit none
! 
!  Passed variables
!
  logical,      intent(in)                          :: lprint
  logical,      intent(in)                          :: lfitrun
  real(dp),     intent(out)                         :: fsumsq
!
!  Local arrays
!
  character(len=20)                                 :: flabel
  character(len=80)                                 :: outdump_name
  character(len=80)                                 :: outeig_name
  character(len=90)                                 :: title_str1
  integer(i4)                                       :: i
  integer(i4)                                       :: ifail
  integer(i4)                                       :: icount
  integer(i4)                                       :: j
  integer(i4)                                       :: ni
  integer(i4)                                       :: qcount
  integer(i4)                                       :: status
  integer(i4)                                       :: shellcount
  integer(i4)                                       :: thetacount
  integer(i4)                                       :: phicount
  integer(i4)                                       :: qpt_numin
!
  logical                                           :: lsavefreqout
  logical                                           :: lprintloc
!
  real(dp)                                          :: g_cpu_time
  real(dp)                                          :: fc
  real(dp)                                          :: phi
  real(dp)                                          :: Qmodu
  real(dp)                                          :: qv(3,3)
  real(dp)                                          :: init_par(4)
  real(dp)                                          :: t1s
  real(dp)                                          :: t2s
  real(dp)                                          :: t1t
  real(dp)                                          :: t2t
  real(dp)                                          :: tnonscatter
  real(dp)                                          :: theta
  real(dp)                                          :: wrk
  real(dp),     dimension(:),     allocatable, save :: dwaller
  real(dp),     dimension(:,:),   allocatable, save :: endarray
  real(dp),     dimension(:),     allocatable, save :: w1
!
  t1t = g_cpu_time()
  tnonscatter = 0.0_dp
  lfirst = .true.
  ldebug = (index(keyword,'debu').ne.0)
  lscatsurfaceweight = .false.
  lprintloc = (lprint.and.ioproc)
!
!  Initialise fsumsq
!
  fsumsq = 0.0_dp
!
!  Save original value of this logical to restore at the end
!
  lsavefreqout = lfreqout
  lfreqout = .false.
!  
  title_str1 = '# Q Step   Q Modulus     Theta       Phi          Mode    Freq (cm-1)             S(Q,w)'
  do i = 1,4
    init_par(i) = 0.0_dp
  enddo

!***********************************************************************************************
!  Local memory allocation section
!***********************************************************************************************

  allocate(dwaller(numat),stat = status)
  if (status/=0) call outofmemory('scatter','dwaller')
!
!  Allocate memory associate with sampling arrays
!
  if (q_runtype .eq. 1) then
!
!  Memory allocation for rs_onion:
!  Array is ARRAYNAME( Individual values of Q, q and tau vector - Qx, Qy and Qz , 
!                     Number of phi increments, Number of theta increments, number of Q increments)
!  TESTING -   set holding variables for rs_onion
    maxHold1 = 3_i4        ! x,y and z components of the vector in the hold array
    maxHold2 = nq_intstep  ! number of phi increments
    maxHold3 = nq_intstep  ! number of theta increments
    maxHold4 = nq_step     ! number of Q shells
    call changemaxhold
!
  elseif (q_runtype .eq. 2) then
!  Memory allocation for ls_set:
!  Array is ARRAYNAME(Individual values of Q, q and tau vector - Qx, Qy and Qz, 
!                     Number of Q increments, Number of vectors in set, 1... as this is redundant)
!  TESTING
    maxHold1 = 3_i4
    maxHold2 = n_qpointmax
    maxHold3 = n_qvector
    maxHold4 = 1_i4
    call changemaxhold
!   
  elseif (q_runtype .eq. 3) then
!  NONFUNCTIONAL
!
  elseif (q_runtype .eq. 4) then
!  NONFUNCTIONAL
!
  elseif (q_runtype .eq. 5) then
!  NONFUNCTIONAL
!
  elseif (q_runtype .eq. 6) then
!  Memory allocation for sc_sample:
!  Array is ARRAYNAME(Individual values of Q, q and tau vector - Qx, Qy and Qz 
!                      Number of Q increments,1... as this is redundant,1... as this is redundant)
!  TESTING
    maxHold1 = 3_i4  ! x,y and z components of the vector in the hold array
    maxHold2 = n_qpointmax ! Number of q points sampled (in this case, for the entire Q vector)
    maxHold3 = 1_i4  ! effectively not used for this config
    maxHold4 = 1_i4  ! effectively not used for this config
    call changemaxhold
!
  elseif (q_runtype .eq. 7) then
!  Memory allocation for pe_file:
!  Array is defined by the type of input file; with reference to the type of switch present.
!  NONFUNCTIONAL
  endif
!
!  Gather scattering lengths into local array based on isotopes  
!
  do i = 1,ncore
    ni = nat(i)
    if (q_caltype.eq.1) then
      scatlencoh(i) = b_coh(1,ni)
      do j = 1,maxisotopes
        if (abs(atmass(ni)-b_mass(j,ni)).lt.0.2_dp) then
          scatlencoh(i) = b_coh(j,ni)
        endif
      enddo
    elseif (q_caltype.eq.2) then
      scatleninc(i) = b_inc(1,ni)
      do j = 1,maxisotopes
        if (abs(atmass(ni)-b_mass(j,ni)).lt.0.2_dp) then
          scatleninc(i) = b_inc(j,ni)
        endif
      enddo
    elseif (q_caltype.eq.3) then
      scatlencoh(i) = b_coh(1,ni)
      scatleninc(i) = b_inc(1,ni)
      do j = 1,maxisotopes
        if (abs(atmass(ni)-b_mass(j,ni)).lt.0.2_dp) then
          scatlencoh(i) = b_coh(j,ni)
          scatleninc(i) = b_inc(j,ni)
        endif
      enddo
    elseif (q_caltype.eq.4) then
      scatlencoh(i) = b_coh(1,ni)
      scatleninc(i) = b_inc(1,ni)
      do j = 1,maxisotopes
        if (abs(atmass(ni)-b_mass(j,ni)).lt.0.2_dp) then
          scatlencoh(i) = b_coh(j,ni)
          scatleninc(i) = b_inc(j,ni)
        endif
      enddo
    elseif (q_caltype.eq.5) then
      scatlencoh(i) = b_coh(1,ni)
      scatleninc(i) = b_inc(1,ni)
      do j = 1,maxisotopes
        if (abs(atmass(ni)-b_mass(j,ni)).lt.0.2_dp) then
          scatlencoh(i) = b_coh(j,ni)
          scatleninc(i) = b_inc(j,ni)
        endif
      enddo
    endif
  enddo
! 
!  Check size for largest number of k points
!  
  if (n_qpointmax.gt.maxskpt.or.nq_intstep.gt.maxskpt) then
    maxskpt = max(n_qpointmax,nq_intstep)
    call changemaxskpt
  endif

!***********************************************************************************************
!   END OF MEMORY ALLOCATION SECTION
!***********************************************************************************************

  call kvector3D
  qv = kv                                        !kt  kv is recirpocal cell vector with 2Pi
  call matrix_inversion(qv,3_i4,3_i4,wrk,ifail)  !kt  by inversion, qv is used to obtain q in the fractional unit in rso
!  Reciprocal Space Onion Sampling Method
  if (q_runtype .eq. 1) then
    call rso(qv,init_par)
! Write Q vector array to appropriate file using type of data format = 1 !kt this is for mpi parallel calc.
    flabel = 'Q'
    call dump_hold(1_i4,q_filename,74_i4,flabel,init_par)
    flabel = 'smq'
    call dump_hold(1_i4,q_filename,75_i4,flabel,init_par)
    flabel = 'tau'
    call dump_hold(1_i4,q_filename,76_i4,flabel,init_par)
!
!  Linear Sampling Set Method
!
  elseif (q_runtype .eq. 2) then
    call lss(init_par)
! Write Q vector array to appropriate file using type of data format =2
    flabel = 'Q'
    call dump_hold(2_i4,q_filename,74_i4,flabel,init_par)
    flabel = 'smq'
    call dump_hold(2_i4,q_filename,75_i4,flabel,init_par)
    flabel = 'tau'
    call dump_hold(2_i4,q_filename,76_i4,flabel,init_par)
!
!  Reciprocal Space Cylinder Method   -   NONFUNCTIONAL
!
  elseif (q_runtype .eq. 3) then
!     call rsc(qv, Hold_Q, Hold_smq, Hold_Tau,init_par)
! Write Q vector array to appropriate file using type of data format = 1
     
!
!  Monkhorst Pack Sampling Method   -   NONFUNCTIONAL
!
  elseif (q_runtype .eq. 4) then
!     call mps(qv, Hold_Q, Hold_smq, Hold_Tau,init_par)
! Write Q vector array to appropriate file using type of data format = 1
!
    
!
!  Monte Carlo Sampling Method   -   NONFUNCTIONAL
!
  elseif (q_runtype .eq. 5) then
!     call mcs(qv, Hold_Q, Hold_smq, Hold_Tau,init_par)
! Write Q vector array to appropriate file using type of data format = 1
  
!
!  Single Crystal method
!
  elseif (q_runtype .eq. 6) then
!  determine Q, q and Tau vectors for each qpoint
    call lss(init_par)
! Write Q vector array to appropriate file using type of data format =6
    flabel = 'Q'
    call dump_hold(6_i4,q_filename,74_i4,flabel,init_par)
    flabel = 'smq'
    call dump_hold(6_i4,q_filename,75_i4,flabel,init_par)
    flabel = 'tau'
    call dump_hold(6_i4,q_filename,76_i4,flabel,init_par)
  endif

!***********************************************************
!    Calculate dynamic Debye-Waller factor for each atom   *
!    in the lattice.    DUMMY                              *
!***********************************************************   
  do i = 1,numat
    dwaller(i) = 1.0_dp
  enddo
  
!*************************************************************************************************
!       
!      SELECT SCATTERING KERNEL AND CALL THE SPECIFIC SCATTERING KERNEL FOR EACH Q_RUNTYPE
!      THEN WRITE THE INTERMEDIATE VALUES TO FILE AND MAIN MEMORY (IN SUBROUTINES)
!
!*************************************************************************************************  

  if (q_runtype .eq. 1) then
    if (nq_intstep.gt.maxskpt) then  ! kt this was checked before. Why check again?
      maxskpt = nq_intstep
      call changemaxskpt
    endif
!
!  Allocate and initialise endarray
!
    allocate(endarray(nw_step,nq_step),stat = status)
    if (status/=0) call outofmemory('scatter','endarray')      
    endarray(1:nw_step,1:nq_step) = 0.0_dp
!
!  Output: For this option, the file needs to be opened outside the loops over variables
!          so that multiple calls to cscatter can be recorded.
!
    if (lprintloc) then
      outdump_name = trim(adjustl(q_filename))//'.sqo'  ! Onion output, in terms of theta and phi
      open(60,file=outdump_name,status='unknown')
      write(60,'(a100,/)') outdump_name
      write(60,'(''# Reciprocal space - RSO sampling output from Scatter calculation :'')')
!
!  Lead in with line giving dimensions of holding array for future read-in
!
      write(60,'(''# Size of file:  No. of shells, Angular steps in theta and phi, Modes per Q point'')')
      write(60,'(3i6)') nq_step, nq_intstep, 3*numat
!
!    End of output lead-in section...
!
!      write(ioout,'('' nq_step = '',i6, ''  nq_intstep = '',i6)') nq_step, nq_intstep
    endif
    nskpt = nq_intstep
    Qmodu = init_par(1)
    do shellcount = 1,nq_step
      if (lprintloc) then
        write(60,'(/,''# Modulus of Q = '',f10.4,'' Inverse Angstroms.'')') Qmodu
        write(60,'(''#----------------------------------------------'')')
        write(60,'(90a,/)') title_str1
      endif
      theta = theta_initial
      do thetacount = 1,nq_intstep ! kt only phi's are processed at once. 
        xskpt     = Hold_smq(1,:,thetacount,shellcount)  
        yskpt     = Hold_smq(2,:,thetacount,shellcount)
        zskpt     = Hold_smq(3,:,thetacount,shellcount)
        Qvector   = Hold_Q(:,:,thetacount,shellcount)
        tauvector = Hold_tau(:,:,thetacount,shellcount)
		wskpt(1:nq_intstep) = 1.0_dp
!
!  Call energy and phonon to compute frequencies and eigenvectors !! kt for each qvec! so it is time-consuming!!
!
        lscattercall = .true.
        lsymderv2 = .false.
        t1s = g_cpu_time()
        call energy(fc,.true.,.true.) ! kt fc is total energy
        call phonon(.false.,fc,0_i4,0_i4) ! kt phonon calls cscatter which calculates sofomega
<<<<<<< HEAD
		sofomega = sofomega * abs(sin(theta)) ! kt I think this sin should be multiplied with sofomega as the directional weight 
		!print *, "kt, chk",  sin(theta), theta  ! kt I think this sin should be multiplied with sofomega as the directional weight 
=======
>>>>>>> 4ee3e8ac7e3f4d4e86c52737f945af3adff2f7ec
        t2s = g_cpu_time()
        tnonscatter = tnonscatter + t2s - t1s
        lscattercall = .false.
!
        qpt_numin = nq_intstep
        call freqhist(endarray,shellcount,qpt_numin)
               
        if (lprintloc) then
          phi = phi_initial
          do phicount = 1,nq_intstep
            do i = 1,3*numat
              write(60,'(i4,4x,f10.4,3x,f10.4,3x,f10.4,5x,i4,3x,f10.4,f25.8)') &
                         shellcount,Qmodu,theta,phi,i,freq(i,phicount),sofomega(i,phicount)
            enddo
            phi = phi + init_par(4)
          enddo ! over phicount the second time...
        endif
        theta = theta + init_par(4)
      enddo  ! over theta
      Qmodu = Qmodu + init_par(3)
    enddo ! over shells
!
!  In parallel perform reduction of endarray
!
    if (nprocs.gt.1) then
      call sumhist(endarray)
    endif
!
!  Close output file      
!
    if (lprintloc) then
      close(60)
    endif
! 
    lscatsurfaceweight = .false.
    call out_rso(endarray,nw_step,nq_step,nq_step,lscatsurfaceweight) ! kt output sqw to a text file.
!
!  Compute sum of squares of difference in endarray for fitting
!
    if (lfitrun) then
      call sqomegadifference(endarray,nw_step,nq_step,nq_step,fsumsq)
    endif
!
! Deallocate kpoints for phonon call     
!
    deallocate(endarray)
   
  elseif (q_runtype .eq. 2) then
!
!  Linear Sampling Set Method
!
    if (n_qpointmax.gt.maxskpt) then
      maxskpt = n_qpointmax
      call changemaxskpt
    endif
!
!  Output: For this option, the file needs to be opened outside the no. of 
!            q_points loop, so that multiple calls to cscatter can be recorded.
!
    if (lprintloc) then
      outdump_name = trim(adjustl(q_filename))//'.sqv'  ! Set of vectors
      open(71,file=outdump_name,status='unknown')
      write(71,'(a100,/)') outdump_name
      write(71,'(''# Reciprocal space - Linear sampling output from Scatter calculation :'')')
!
!  Lead in with line giving dimensions of holding array for future read-in
!
      write(71,'(''# Please Note: all outvector components are (lattice) fractional'')')
      write(71,'(''# Size of file:  Number of vectors, largest number of increments'')')
      write(71,'(2i6,/)') n_qvector,n_qpointmax
!
!  End of output lead-in section...
!
    endif
!
!   LOOP OVER NUMBER OF VECTORS
!
    do qcount = 1,n_qvector
!
!  Collect from holding arrays and assign to x,y,zskpt(nskpt) and Qvector arrays      
!
      do icount = 1, n_qpoint(qcount)
        nskpt = n_qpoint(qcount)
        xskpt(icount)       = Hold_smq(1,icount,qcount,1)  
        yskpt(icount)       = Hold_smq(2,icount,qcount,1)
        zskpt(icount)       = Hold_smq(3,icount,qcount,1)
        Qvector(1,icount)   = Hold_Q(1,icount,qcount,1)
        Qvector(2,icount)   = Hold_Q(2,icount,qcount,1)
        Qvector(3,icount)   = Hold_Q(3,icount,qcount,1)
        tauvector(1,icount) = Hold_tau(1,icount,qcount,1)
        tauvector(2,icount) = Hold_tau(2,icount,qcount,1)
        tauvector(3,icount) = Hold_tau(3,icount,qcount,1)
!
!  wkpt needs to be present, but is set to unity for this configuration
!
        wskpt(icount) = 1.0_dp
        if (ldebug.and.ioproc) then
!
!  Check to see if the values passed are sensible        
!
          write(ioout,'(3f16.8)') xskpt(icount),yskpt(icount),zskpt(icount)
          write(ioout,'(3f16.8)') Qvector(1,icount),Qvector(2,icount),Qvector(3,icount)
        endif
      enddo
!
!  Call energy and phonon to compute frequencies and eigenvectors
!
      lscatsurfaceweight = .true.
      lscattercall = .true.
      lsymderv2 = .false.
      t1s = g_cpu_time()
      call energy(fc,.true.,.true.)
      call phonon(.false.,fc,0_i4,0_i4)
      t2s = g_cpu_time()
      tnonscatter = tnonscatter + t2s - t1s
      lscattercall = .false.
      if (lprintloc) then
!
!  Write the data in 'data_array' into the file
!
        write(71,'(/,''Vector '',f10.4,2x,f10.4,2x,f10.4,'' to '',f10.4,2x,f10.4,2x,f10.4,'' in '',i6, '' incs'')') &
                       q_ordinate(1,qcount),q_ordinate(2,qcount),q_ordinate(3,qcount), &
                       q_ordinate(4,qcount),q_ordinate(5,qcount),q_ordinate(6,qcount),n_qpoint(qcount)
        write(71,'(''#-----------------------------------------------------------------------------------------------'',/)')
        do i = 1, n_qpoint(qcount)
          write(71,'(/,''# Q point  Mode             S(Q,w)    '',/)')
!          write(71,'(''#--------------------------------------'',/)')
          do j = 1,3*numat
            write(71,'(2i6,f25.8)')i,j, sofomega(j,i)
          enddo
        enddo
      endif
    enddo !n_qvector
    if (lprintloc) then
      close(71)
    endif
!   
!   END OF LOOP OVER NUMBER OF VECTORS
!
  elseif (q_runtype .eq. 3) then
!
!  Reciprocal Space Cylinder Method   -   NONFUNCTIONAL
!   
  elseif (q_runtype .eq. 4) then
!
!  Monkhorst Pack Sampling Method     -   NONFUNCTIONAL
!   
  elseif (q_runtype .eq. 5) then
!
!  Monte Carlo Sampling Method        -   NONFUNCTIONAL
!   
  elseif (q_runtype .eq. 6) then
!
!  Single Crystal method
!

!
!  Check sizes of kpoint arrays
!
    if (n_qpoint(1).gt.maxskpt) then
      maxskpt = n_qpoint(1)
      call changemaxskpt
    endif
!
!  Allocate and initialise endarray
!
    allocate(endarray(nw_step,n_qpoint(1)),stat = status)
    if (status/=0) call outofmemory('scatter','endarray') 
    endarray(1:nw_step,1:n_qpoint(1)) = 0.0_dp
!
!  Collect from holding arrays and assign to x,y,zskpt(nskpt) and Qvector arrays      
!
    do icount = 1,n_qpoint(1)
      xskpt(icount)       = Hold_smq(1,icount,1,1)  
      yskpt(icount)       = Hold_smq(2,icount,1,1)
      zskpt(icount)       = Hold_smq(3,icount,1,1)
      Qvector(1,icount)   = Hold_Q(1,icount,1,1)
      Qvector(2,icount)   = Hold_Q(2,icount,1,1)
      Qvector(3,icount)   = Hold_Q(3,icount,1,1)
      tauvector(1,icount) = Hold_tau(1,icount,1,1)
      tauvector(2,icount) = Hold_tau(2,icount,1,1)
      tauvector(3,icount) = Hold_tau(3,icount,1,1)
!
!  wkpt needs to be present, but is set to unity for this configuration
!
      wskpt(icount) = 1.0_dp
      if (ldebug.and.ioproc) then
!
!  Check to see if the values passed are sensible        
!
        write(ioout,'(3f16.8)') xskpt(icount),yskpt(icount),zskpt(icount)
        write(ioout,'(3f16.8)') Qvector(1,icount),Qvector(2,icount),Qvector(3,icount)
      endif
    enddo
!
!  Block to set globals for phonon prior to call
!
    nskpt = n_qpoint(1)
!
!  Output: For this option, eigenvectors are output in the .sqei file: 
!
!  - initialise file
!  - write header
!  - data written in cscatter
!
    if (lprintloc) then
      outeig_name = trim(adjustl(q_filename))//'.sqei'  ! Single
      open(69,file=outeig_name,status='unknown')
      write(69,'(a100,/)') trim(adjustl(outeig_name))
      write(69,'(''# Single direction Eigenvector Output:'')')
      write(69,'(''# Format:  Q Point, Eigenvectors with associated frequency. '')')
    endif
!
!  Call to energy then phonon to compute frequencies and eigenvectors
!
    lscatsurfaceweight = .true.
    lscattercall = .true.
    lsymderv2 = .false.
    t1s = g_cpu_time()
    call energy(fc,.true.,.true.)
    call phonon(.false.,fc,0_i4,0_i4)
    t2s = g_cpu_time()
    tnonscatter = tnonscatter + t2s - t1s
    lscattercall = .false.
!
    call freqhist(endarray,0_i4,n_qpoint(1))
!
!  In parallel perform reduction of endarray
!
    if (nprocs.gt.1) then
      call sumhist(endarray)
    endif

    if (ldebug.and.ioproc) then
      do i = 1,n_qpoint(1)
        do j = 1,nw_step
          write(ioout,'(f36.10)') endarray(j,i)
        enddo
      enddo
    endif
    call out_rso(endarray,nw_step,nskpt,n_qpoint(1),lscatsurfaceweight)
!
!  Compute sum of squares of difference in endarray for fitting
!
    if (lfitrun) then
      call sqomegadifference(endarray,nw_step,nskpt,n_qpoint(1),fsumsq)
    endif
!
!  Close channel 69 since all eigenvectors have been written
!
    if (lprintloc) then
      close(69) ! for eigenvector output
    endif
!
!  Check output intensities from sofomega
!
    if (ldebug.and.ioproc) then
      do i = 1, n_qpoint(1)
        write(ioout,'(''q point: '',i4)') i
        do j = 1,3*numat
          write(ioout,'(''Freq for mode '',i4,'' is '',f25.8,/)') j,freq(j,i)
          write(ioout,'(''S(Q,w) for mode '',i4,'' is '',f25.8,/)') j,sofomega(j,i)
        enddo
      enddo
    endif ! debug
!     
!  Output for S(Q,w) in standard format
!
    if (lprintloc) then
      outdump_name = trim(adjustl(q_filename))//'.sqs'  ! Single
      open(70,file=outdump_name,status='unknown')
      write(70,'(a100,/)') outdump_name
      write(70,'(''# Reciprocal space - single dir. output from Scatter calculation :'')')
      write(70,'(''# Format =  individual increment, outvector (x,y,z)'')')
!
!  Lead in with line giving dimensions of holding array for future read-in
!
      write(70,'(''# Please Note: all outvector components are (lattice) fractional.'',/)')
!
!  Then write the data in 'data_array' into the file
!
      write(70,'(''# Vector is gamma point  to '',f10.4,''   '',f10.4,''   '',f10.4,''In '',i6, '' increments. '')') &
                      q_ordinate(4,1),q_ordinate(5,1),q_ordinate(6,1),n_qpoint(1)
      write(70,'(''#-------------------------------------------------------------------------------------------------'',/)')
      do i = 1,n_qpoint(1)
        write(70,'(/,''# Q point (x,y,z) =  '', 3f16.8)') &
                      Hold_Q(1,i,1,1),Hold_Q(2,i,1,1),Hold_Q(3,i,1,1)
        write(70,'(''#---------------------------------------------------------------------'',/)')
        write(70,'(''# Mode               Omega        S(Q,w)        '')')
        write(70,'(''#----------------------------------'',/)')
        do j = 1,3*numat
          write(70,'(i6,2x,f24.10,2x,f24.10)')j,freq(j,i),sofomega(j,i)
        enddo
      enddo
      close(70)
    endif
!
!  Deallocate kpoints for phonon call     
!
    deallocate(endarray)
     
  elseif (q_runtype .eq. 7) then
!
!  Read existing sampling from file   -   NONFUCTIONAL
!
  endif ! selection and calling of scattering kernel by sampling type

!*************************************************************************************************
!       
!      SELECT APPROPRIATE SCATTERING OPTIONS, PERFORM CALCULATIONS AND WRITE TO FILE 
!      (IN SUBROUTINES), THEN CONVOLVE WITH MAIN ARRAY VALUES IF FLAGS ARE SET
!
!************************************************************************************************* 

   if (lscatopt(1)) then
   
   elseif (lscatopt(2)) then
   
   elseif (lscatopt(3)) then
   
   elseif (lscatopt(4)) then
   
   elseif (lscatopt(5)) then
   
   elseif (lscatopt(6)) then
   
   elseif (lscatopt(7)) then
   
   elseif (lscatopt(8)) then

   else
   
   endif

!*************************************************************************************************
!       
!      SELECT APPROPRIATE OUTPUT AND PRE-PROCESSING OPTIONS, PERFORM CALCULATIONS AND WRITE  
!      TO FILE (IN SUBROUTINES), THEN CONVOLVE WITH MAIN ARRAY VALUES IF FLAGS ARE SET
!
!************************************************************************************************* 

   if (lscatanopt(1)) then
   
   elseif (lscatanopt(2)) then
   
   elseif (lscatanopt(3)) then
   
   elseif (lscatanopt(4)) then
   
   elseif (lscatanopt(5)) then
   
   elseif (lscatanopt(6)) then
   
   elseif (lscatanopt(7)) then

   else
!***********************************************************
!                 DUMMY put primary output here            *
!***********************************************************   
   endif

!     
!  Free local memory
!     
   if (allocated(w1)) deallocate(w1)
   if (allocated(dwaller)) deallocate(dwaller)
!
!  Restore original value of this logical 
!
  lfreqout = lsavefreqout
!
   t2t = g_cpu_time()
   tscatter = t2t - t1t + tscatter - tnonscatter
   return
   end
!
        

!*********************************************************************
!*                                                                   *
!*                        SUBROUTINES FOR SCATTER.F90                *
!*                                                                   *
!*********************************************************************



!***************************************************************************
!  Subroutine Linear Set Sampling (LSS): Generates a set of Q vectors,
!  and corresponsing q vectors, from Qmodu(min) to Qmodu(max).  Also returns
!  a set of matching Tau (r.l.) vectors for each (Q,q) set, giving Q,q,tau
!  and a sampling area size for each direction selected in reciprocal space.
!
!       NEW FOR f90
!
!***************************************************************************

  subroutine lss(work_initpar)
!
  use datatypes
  use current
  use constants
  use scatterdata
  
  implicit none
!
!  Passed variables
!
  real(dp),     intent(inout)  :: work_initpar(4)
!
!  Local variables  
!
  integer(i4)                  :: i
  integer(i4)                  :: j
  integer(i4)                  :: inc_count
  integer(i4)                  :: increments
  real(dp)                     :: work_Q(3)
  real(dp)                     :: work_smq(3)
  real(dp)                     :: work_tau(3)
  real(dp)                     :: remain_vec(3)
  real(dp)                     :: init_vec(3)
  real(dp)                     :: end_vec(3)
  real(dp)                     :: inc_size(3) 
!
!  For each n_qvector, assign a beginning and end vector to the working vector variables,
!    then cycle through the increments and populate the holding array...
!
  do i = 1,n_qvector
    increments = n_qpoint(i)
!
!  Read in Q vectors from main array, and assign to working arrays
!
    do j = 1,3
      init_vec(j) = q_ordinate(j,i)
      end_vec(j)  = q_ordinate(j+3,i)
      inc_size(j) = (end_vec(j) - init_vec(j))/dble(increments)
!
!  Initial assignment outside increment loop of Working array holding Q vector components
!
      work_Q(j)   = init_vec(j)
    enddo
!
!  Loop over increments, adding individual vector elements to Work_Q, then determining q and tau,
!  storing the x,y and z elements of Q, smq and tau in the holding arrays.
!
    do inc_count = 1,increments
      do j = 1,3
        remain_vec(j) = mod(work_Q(j),1.0_dp)
        if (remain_vec(j).lt.0.5_dp) then
          work_smq(j) = remain_vec(j)
        else
          work_smq(j) = remain_vec(j) - 1.0_dp
        endif
        work_tau(j) = work_Q(j) - work_smq(j)
!
!  Assign to global holding arrays
!
        Hold_Q(j,inc_count,i,1)    = work_Q(j)
        Hold_smq(j,inc_count,i,1)  = work_smq(j)
        Hold_tau(j,inc_count,i,1)  = work_tau(j)
!
!  Increment Q
!
        work_Q(j)   = work_Q(j) + inc_size(j)           
          
      enddo ! over j

    enddo  !  over inc_count
  
  enddo ! over i (number of vectors)

  work_initpar(1) = 1.0_dp
  work_initpar(2) = 1.0_dp
  work_initpar(3) = 1.0_dp
  work_initpar(4) = 1.0_dp

  return
  end

!***************************************************************************
!  Subroutine Reciprocal Space Onion (RSO): Generates a set of Q vectors,
!  and corresponsing q vectors, from Qmodu(min) to Qmodu(max).  Also returns
!  a set of matching Tau (r.l.) vectors for each (Q,q) set, giving Q,q,tau
!  and a sampling area size for each layer of the 'onion'
!
!       NEW FOR f90, but modified from old routine (then in main body of code)
!
!***************************************************************************

  subroutine rso(qv_in, work_initpar)
!
  use datatypes
  use current
  use constants
  use scatterdata
  use iochannels
  implicit none
!
!  Passed variables
!
  real(dp),     intent(inout)  :: work_initpar(4)
  real(dp),     intent(in)     :: qv_in(3,3)
!
!  Local variables
!
  integer(i4)                  :: i
  integer(i4)                  :: layercount
  integer(i4)                  :: thetacount
  integer(i4)                  :: phicount
  integer(i4)                  :: gencount
  integer(i4)                  :: qloop_step
  real(dp)                     :: Qmodulus
  real(dp)                     :: rsteps
  real(dp)                     :: work_vec(3)
  real(dp)                     :: work_Q(3)
  real(dp)                     :: work_smq(3)
  real(dp)                     :: work_tau(3)
  real(dp)                     :: RemainQ(3)
  real(dp)                     :: Qinc
  real(dp)                     :: Q_init
  real(dp)                     :: Q_final
  real(dp)                     :: tinc
  real(dp)                     :: pinc
  real(dp)                     :: theta
  real(dp)                     :: phi
  
!***********************************************************
!*   Angular sampling process for reciprocal space onion   *
!*   in terms of a loop over Q (outer loop), then          *
!*   theta (middle loop) and phi (inner loop).             *
!***********************************************************
  
  rsteps = dble(nq_intstep)
  Qinc = (q_qmax-q_qmin)/dble(nq_step)
  tinc = pi/rsteps
  pinc = pi/rsteps
!
!  Set Qmodulus for zero and non-zero q_qmin  
!
<<<<<<< HEAD
  !kt if (q_qmin .gt. globenull_chk) then
  !kt   Qmodulus =q_qmin
  !kt   Q_init = Qmodulus
  !kt   nq_step = nq_step + 1
  !kt   if (nq_step.gt.maxHold4) then
  !kt     maxHold4 = nq_step
  !kt     call changemaxhold
  !kt   endif
  !kt else
  !kt   Qmodulus = q_qmin + Qinc
  !kt   Q_init = Qmodulus
  !kt endif !if larger than zero 
  Qmodulus = q_qmin + Qinc                 !kt I do not like the modification above.
  Q_init = Qmodulus
=======
  if (q_qmin .gt. globenull_chk) then
    Qmodulus =q_qmin
    Q_init = Qmodulus
    nq_step = nq_step + 1
    if (nq_step.gt.maxHold4) then
      maxHold4 = nq_step
      call changemaxhold
    endif
  else
    Qmodulus = q_qmin + Qinc
    Q_init = Qmodulus
  endif !if larger than zero 
>>>>>>> 4ee3e8ac7e3f4d4e86c52737f945af3adff2f7ec

  qloop_step = nq_step

  do layercount = 1,qloop_step
    theta = theta_initial
    do thetacount = 1,nq_intstep
      phi = phi_initial
      do phicount = 1, nq_intstep
           
        work_vec(1) = Qmodulus*dsin(theta)*dcos(phi)
        work_vec(2) = Qmodulus*dsin(theta)*dsin(phi)
        work_vec(3) = Qmodulus*dcos(theta)
		work_Q = matmul(qv_in, work_vec)
		work_smq = mod(work_Q, 1.0_dp) ! work_smq is now from -1 to 1
		where (mod(work_smq, 1.0_dp)>=0.5_dp) work_smq(:) = work_smq(:) - 1.0_dp  ! work_smq is now from -1 to 0.5. I checked that wherther -0.5 to
			                                                                      ! 0.5 or -1 to 0.5 just causes negligible numeriacl errors in sqw of SiO2.
		work_tau = work_Q - work_smq                                              ! work_tau is a reciprocal lattice vector which  connects work_smq and work_Q
!
!  Store sampling of (Q,w) space in Holding arrays - remember that arrays store fractional x,y,z
!
        Hold_Q(:,phicount,thetacount,layercount)   = work_Q(:)
        Hold_smq(:,phicount,thetacount,layercount) = work_smq(:)
        Hold_Tau(:,phicount,thetacount,layercount) = work_tau(:)
        phi = phi + pinc    
      enddo  ! over phi
      theta = theta + tinc  
    enddo !  End loop over theta and hence entire angular integration 
    if (Qmodulus .lt. q_qmax) then
      Qmodulus = Qmodulus + Qinc 
      Q_final = Qmodulus         
    endif
  enddo
!
!  End loop over modulus of Q
!  output init_par info for dumphold routine
!
  work_initpar(1) = Q_init
  work_initpar(2) = Q_final
  work_initpar(3) = Qinc
  work_initpar(4) = tinc
  return

  end

!**********************************************************
!*   This routine calculates the Debye-Waller Factor W    *
!*   for the crystal system, for use in scatter           *
!*                                                        *
!*      UNDER COMPLETE OVERHAUL - NON-ESSENTIAL           *
!**********************************************************

!      subroutine dwalcalc(
!     &            )
        
!      use constants
!      use scatterdata
!      implicit none
!
!  Passed arguments
!
! 
!
!      return
!      end

!**********************************************************
!*   This routine calculates S(q,w) inc for a given Q, q  *
!*   e(ds) and scattering length b.                       *
!**********************************************************
! kt wrong!!

  subroutine iscatter(numb,realeds,imageds,maxd2,omegas,Qvect,sQwout,phaser,phasei,thetain,&
                      scatlen,dtheta,dphi,Qmodu)
                     
  use constants   
  use current 
  use scatterdata

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)   :: maxd2
  integer(i4), intent(in)   :: numb
  real(dp)                  :: Qvect(3), omegas(3*numb), sQwout(3*numb),realeds(maxd2,*)
  real(dp)                  :: imageds(maxd2,*),phaser(3*numb),phasei(3*numb),thetain
  real(dp)                  :: scatlen(*),dtheta,dphi,Qmodu
!
!  Local variables
!
  integer(i4)               :: i,j,jx,jy,jz
  logical                   :: flatflag
  real(dp)                  :: FofQr,FofQi,qedotr,qedoti,FofQ2,qcdot,partr,parti,cscatcon
!
!  Define constants
!
  cscatcon = 1.0_dp 
  flatflag = .false.
!
!  Calculate atomic phase factors
!
  do i = 1,numb
    qcdot = Qvect(1)*xclat(i) + Qvect(2)*yclat(i) + Qvect(3)*zclat(i)
    phaser(i) = cos(qcdot)
    phasei(i) = sin(qcdot)
  enddo
!
!  Outer loop over modes
!
  do i = 1,3*numb
!
!  Inner loop over components of eigenvector
!
    FofQr = 0.0_dp
    FofQi = 0.0_dp
    jx = -2
    jy = -1
    jz = 0
    do j = 1,numb
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!
!  Eigenvector . qvector for each atom
!
      qedotr = Qvect(1)*realeds(jx,i) + Qvect(2)*realeds(jy,i) + Qvect(3)*realeds(jz,i)
      qedoti = Qvect(1)*imageds(jx,i) + Qvect(2)*imageds(jy,i) + Qvect(3)*imageds(jz,i)
      partr = qedotr*phaser(j) - qedoti*phasei(j)
      parti = qedoti*phaser(j) + qedotr*phasei(j)
      FofQr = FofQr + partr*scatlen(j)
      FofQi = FofQi + parti*scatlen(j)
    enddo
    FofQ2 = FofQr*FofQr + FofQi*FofQi
    if (flatflag) then
      sQwout(i) = (dtheta*dphi*dsin(thetain)*FofQ2*cscatcon)/(3*dble(numb)*Qmodu*Qmodu) 
    else
      sQwout(i) = (dtheta*dphi*dsin(thetain)*FofQ2*cscatcon)/(3*dble(numb)*omegas(i))
    endif
  enddo

  return
  end

!*****************************************************************************
!*   This routine calculates the frequency histogram by boxing the data into *
!*   discrete intervals for return to main code.                             *
!*                                                                           *
!*       UPDATED FOR F90                                                     *
!*                                                                           *
!*****************************************************************************

  subroutine freqhist(outarray,currentQ,qpoint_number)

  use scatterdata
  use frequencies
  use constants   
  use current 
  use parallel,    only : nprocs, procid

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: currentQ
  integer(i4), intent(in)                      :: qpoint_number
  real(dp),    intent(inout)                   :: outarray(nw_step,*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: freqbox
  real(dp)                                     :: freqstep
!
!  Single vector outarray option
!
  if (currentQ .eq. 0_i4) then
    freqstep = (q_wmax-q_wmin)/dble(nw_step)
    do i = procid+1,qpoint_number,nprocs
      do j = 1,3*numat
        freqbox = nint(freq(j,i)/freqstep) + 1
        if (freqbox.le.nw_step.and.freqbox.gt.0) then
          outarray(freqbox,i) = outarray(freqbox,i) + sofomega(j,i)
        endif
      enddo
    enddo      
  else
!
!  For rso output get the frequencies, find the correct 'box' and add
!  the associated scattering value, sofomega, to the contents of that 'box'.    
!
    freqstep = (q_wmax-q_wmin)/dble(nw_step)
    do i = procid+1,qpoint_number,nprocs
      do j = 1,3*numat
        freqbox = nint(freq(j,i)/freqstep) + 1
        if (freqbox.le.nw_step.and.freqbox.gt.0) then
          outarray(freqbox,currentQ) = outarray(freqbox,currentQ) + sofomega(j,i)
        endif
      enddo
    enddo
  endif
!
  return
  end

!*****************************************************************************
!*     Out_rso: Get the array of S(Q,w) from RSO and then box it into energy *
!*     intervals set by the 'resolution' of the x axis.  The boxing is       *
!*     also done for the Q axis, although this is seen to be naturally       *
!*     boxed by the previous routine, Freqhist.                              *
!*                                                                           *
!*       UPDATED FOR F90                                                     *
!*                                                                           *
!*****************************************************************************

  subroutine out_rso(sqwarray, nsqwarraylhs, nsqwarrayrhs, q_vectors_in, sflagin)   
      
  use scatterdata
  use constants
  use iochannels
  use parallel,    only : ioproc
      
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                       :: nsqwarraylhs
  integer(i4), intent(in)                       :: nsqwarrayrhs
  real(dp),    intent(in)                       :: sqwarray(nsqwarraylhs,nsqwarrayrhs)
  logical,     intent(in)                       :: sflagin
  integer(i4), intent(in)                       :: q_vectors_in
!
!  Local variables
!
  character(len=85)                             :: singlefn
  character(len=85)                             :: sqwfname
  integer(i4)                                   :: i
  integer(i4)                                   :: j
  logical                                       :: hertzflag
  real(dp)                                      :: qi
  real(dp)                                      :: qscale
  real(dp)                                      :: vertcon
  real(dp)                                      :: freqscale
  real(dp)                                      :: qmin
  real(dp)                                      :: qmax
!
!  Only write from ioproc
!
  if (.not.ioproc) return
!
  hertzflag = .false.
!
!  name & trim off trailing blanks from both sqwfname and
!  singlefn for polycrystal output and single crystal direction output 
!  data .sqw and .xy file
!
  singlefn = trim(adjustl(q_filename))//'.sxy'
  sqwfname = trim(adjustl(q_filename))//'.sqw'
!
!  nb, if unitflag, then units of energy will be in THz, otherwise
!  they will be in cm-1.  
!
  if (hertzflag) then
    vertcon = 0.02998_dp
  else
    vertcon = 1.0_dp
  endif
!
!  Select between single vector and other sampling outputs
!
! kt I do not understand the writing procedure if (sflagin). It just print q and omega mesh.
  if (sflagin) then
    freqscale = (q_wmax - q_wmin)/dble(nw_step)
    qmin = 0.0_dp
!   qmin = sqrt(q_ordinate(1,1)*q_ordinate(1,1)+q_ordinate(2,1)*q_ordinate(2,1)+ &
!                     q_ordinate(3,1)*q_ordinate(3,1))
    qmax = sqrt(q_ordinate(4,1)*q_ordinate(4,1)+q_ordinate(5,1)*q_ordinate(5,1)+ &
                q_ordinate(6,1)*q_ordinate(6,1))           
    qscale = (qmax - qmin)/dble(q_vectors_in)
    open(72,file=singlefn,status='unknown')
    write(72,'(''# Single direction scattering (.xy file) :'')')
    write(72,'(''# Format = Q, Omega (x,y)'')')
         
    do i = 1,q_vectors_in
      qi = q_qmin + i*qscale
      do j = 1,nw_step
        if (sqwarray(j,i).gt.globenull_chk) then
          write(72,'(2f25.6)') qi,vertcon*dble(j)*freqscale
        endif
      enddo
    enddo
    close(72)
  else
    open(71,file=sqwfname,status='unknown')
    freqscale = (q_wmax - q_wmin)/dble(nw_step)
    qscale    = (q_qmax - q_qmin)/dble(nsqwarrayrhs)
    write(71,'(''# Output from scattering calculation : '',a85)') sqwfname
    write(71,'(''# Q '',12x,''Omega (cm-1)'',10x,''S(Q,omega)'')')
    write(71,'(''#-------------------------------------------------------------'',/)')
    write(71,'(i8,1x,i8)') nw_step, q_vectors_in
    do i = 1,q_vectors_in
      qi = q_qmin + i*qscale
      do j = 1,nw_step
        write(71,'(f12.6,4x,f16.8,4x,f35.8)') qi,vertcon*dble(j)*freqscale,sqwarray(j,i)
      enddo
    enddo      
    close(71)
  endif
  end

!******************************************************************************************
!                              Subroutine dump_hold                                       *
!                                                                                         *
! Given the Q array provided by the sampling method selected, create a file labelled      *
! appropriately, write the Q vectors, along with co-ordinate parameters in the scheme     *
! selected, to the file, then close and save.                                             *
!                                                                                         *
!       NEW FOR f90                                                                       *
!                                                                                         *
!                                                                                         *
!******************************************************************************************

  subroutine dump_hold(data_struct,filename,outchan,data_label,init_in)
!  
  use scatterdata
  use parallel,     only : ioproc

  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)  :: data_struct
  integer(i4),       intent(in)  :: outchan
  real(dp),          intent(in)  :: init_in(4)
  character(len=80), intent(in)  :: filename
  character(len=20), intent(in)  :: data_label
!
!  Local variables
!
  integer(i4)                    :: icounter
  integer(i4)                    :: jcounter
  integer(i4)                    :: kcounter
  integer(i4)                    :: gencount
!
!  Formatting parameters for outputting the data (and subsequent read-ins)
!
  character(len=100), parameter :: lss_header = '(i5)'
  character(len=100), parameter :: sc_header = '()'
  character(len=100), parameter :: rso_header = '(3f16.8,2i4,f16.8,2i4,f16.8)'
  character(len=100), parameter :: lss_title = '(i5)'
  character(len=100), parameter :: sc_title = '()'
  character(len=100), parameter :: rso_title = '()'
   
  character(len=100) :: Qdump_name
  character(len=100) :: title_string
  character(len=70)  :: title_string_1
  character(len=70)  :: title_string_2
  character(len=70)  :: title_string_3
  character(len=70)  :: title_string_4
    
  
  real(dp)           :: Qmodulus
  real(dp)           :: ang_size
  real(dp)           :: curr_Q(3)
  real(dp)           :: curr_smq(3)
  real(dp)           :: curr_tau(3)
  real(dp)           :: curr_out(3)
  real(dp)           :: theta
  real(dp)           :: phi
  real(dp)           :: tinc
  real(dp)           :: Qinc
!
!  Only write from ioproc
!
  if (ioproc) return
!
!  NB => Channels are 74 (Q), 75 (smq) and 76 (tau) for output options
!
!  File header setup
!
  title_string_1 = '  Qx             Qy              Qz'
  title_string_2 = 'qx             qy              qz'
  title_string_3 = 'taux           tauy            tauz'
  title_string_4 = '# Q Step     theta             phi'
!
! work_initpar = [q_qmin,q_qmax,Qinc,tinc]
!  If type of data structure matches that of Q, tau or small q (rso data output):
!
  if (data_struct.eq.1) then
!
!  Name & trim off trailing blanks from input filename, then assign .rss (rec. space sample)
!
    Qdump_name = trim(filename)//'_'//trim(data_label)//'.rso'  ! onion
    open(outchan,file=Qdump_name,status='unknown')
    write(outchan,'(a80,/)') trim(adjustl(Qdump_name))
    write(outchan,'(''# Reciprocal space onion sampling output from Scatter calculation :'',/)')
    write(outchan,'(''# Format = Magnitude of Q (inv. angs) theta(deg) phi(deg) outvector (x,y,z)'',/)')
    write(outchan,'(''# Please Note: all outvector components are in a (lattice) fractional format.'')')
!
!  Lead in with line giving dimensions of holding array for future read-in
!
    write(outchan,'(/,''# Size of file: Minimum Q, Maximum Q, Q increment size,No. of Q shells, '')')
!
!  Calculate ang_size from globals
!
    ang_size = 180.0_dp/dble(nq_intstep)
    write(outchan,'(''# number of angular steps in theta and phi and size of each angular step (degrees): '')')
    write(outchan,rso_header) q_qmin,q_qmax,init_in(1),nq_step,nq_intstep,ang_size
!
!  Set up modulus output
!    
    Qmodulus = init_in(1)
    Qinc = init_in(3)
    tinc = init_in(4)
    if (outchan .eq. 74) then
      title_string_1 = '              ' // trim(title_string_1)
      title_string = trim(title_string_4) // trim(title_string_1)
    elseif (outchan .eq. 75) then
      title_string_2 = '              ' // trim(title_string_2)
      title_string = trim(title_string_4) // trim(title_string_2)
    elseif (outchan .eq. 76) then
      title_string_3 = '              ' // trim(title_string_3)
      title_string = trim(title_string_4) // trim(title_string_3)
    endif  
    do icounter = 1, nq_step
      write(outchan,'(/)')
      write(outchan,'(''# Q modulus = '',1f10.4, '' Inverse Angstroms. '')') Qmodulus
      write(outchan,'(''#-------------------------------------------'',/)')
      write(outchan,'(/)')
      write(outchan,'(90a)') title_string
!
!  Integrate over pi rads
!
      theta = theta_initial
      do jcounter = 1, nq_intstep
        phi = phi_initial
!
!  Integrate over pi rads (centre of inversion-half a sphere)
!
        do kcounter = 1, nq_intstep
!
!  Access holding arrays, and write contents to properly formatted file
!
          do gencount = 1,3 
            curr_Q(gencount)   = Hold_Q  (gencount,kcounter,jcounter,icounter)
            curr_smq(gencount) = Hold_smq(gencount,kcounter,jcounter,icounter)
            curr_tau(gencount) = Hold_Tau(gencount,kcounter,jcounter,icounter)
          enddo
          if (outchan .eq. 74) then
            curr_out(1) = curr_Q(1)
            curr_out(2) = curr_Q(2)
            curr_out(3) = curr_Q(3)
          elseif (outchan .eq. 75) then
            curr_out(1) = curr_smq(1)
            curr_out(2) = curr_smq(2)
            curr_out(3) = curr_smq(3)
          elseif (outchan .eq. 76) then
            curr_out(1) = curr_tau(1)
            curr_out(2) = curr_tau(2)
            curr_out(3) = curr_tau(3)
          endif 
          write(outchan,'(i6,5f16.6,f10.4)') icounter,theta,phi,curr_out(1),curr_out(2),curr_out(3)
!
!  Increment phi        
!
          phi = phi + tinc
        enddo
!
!  Increment theta
!
        theta = theta + tinc
      enddo
!
!  End loop over angular integration 
!
      if (Qmodulus .lt. q_qmax) then
        Qmodulus = Qmodulus + Qinc
      endif
    enddo      !*  End loop over modulus of Q     
!
!  End of RSO format
!

!
!  If type of data structure matches that of LSS:
!
  elseif (data_struct.eq.2) then
    Qdump_name = trim(filename)//'_'//trim(data_label)//'.lss' ! linear set
    open(outchan,file=Qdump_name,status='unknown')
    write(outchan,'(a100,/)') trim(adjustl(Qdump_name))
    write(outchan,'(''# Reciprocal space - Linear sampling output from Scatter calculation :'',/)')
    write(outchan,'(''# Format = Vector (number), individual increment, outvector (x,y,z)'',/)')
!
!  Lead in with line giving dimensions of holding array for future read-in
!
    write(outchan,'(''# Please Note: all outvector components are (lattice) fractional'',/)')
    write(outchan,'(''# Size of file:  Number of vectors, largest number of increments'')')
    write(outchan,'(/,2i6)') n_qvector, n_qpointmax

    do icounter = 1, n_qvector
      write(outchan,'(/,''# Vector '',i4, '' :'',3f16.6)') icounter,q_ordinate(1,icounter), &
            q_ordinate(2,icounter),q_ordinate(3,icounter)
      write(outchan,'(''   to'',10x,3f16.6)') &
            q_ordinate(4,icounter),q_ordinate(5,icounter),q_ordinate(6,icounter)
      write(outchan,'(''   in '',i6, '' increments. '',/)') n_qpoint(icounter)

      if (outchan .eq. 74) then
        title_string = '# Vector   Point       ' // title_string_1
        write(outchan,'(60a,/)') title_string
      elseif (outchan .eq. 75) then
        title_string = '# Vector   Point       ' // title_string_2
        write(outchan,'(60a,/)') title_string
      elseif (outchan .eq. 76) then
        title_string = '# Vector   Point       ' // title_string_3
        write(outchan,'(60a,/)') title_string
      endif  
      do jcounter = 1, n_qpoint(icounter)
        do gencount = 1,3 
          curr_Q(gencount)   = Hold_Q(gencount,jcounter,icounter,1)  
          curr_smq(gencount) = Hold_smq(gencount,jcounter,icounter,1)
          curr_tau(gencount) = Hold_tau(gencount,jcounter,icounter,1)
        enddo
        if (outchan .eq. 74) then
          curr_out(1) = curr_Q(1)
          curr_out(2) = curr_Q(2)
          curr_out(3) = curr_Q(3)
        elseif (outchan .eq. 75) then
          curr_out(1) = curr_smq(1)
          curr_out(2) = curr_smq(2)
          curr_out(3) = curr_smq(3)
        elseif (outchan .eq. 76) then
          curr_out(1) = curr_tau(1)
          curr_out(2) = curr_tau(2)
          curr_out(3) = curr_tau(3)
        endif  
        write(outchan,'(2i6,3f16.6)') icounter,jcounter,curr_out(1),curr_out(2),curr_out(3)
      enddo

    enddo
!
!  End of LSS format
!

!
!  If type of data structure matches that of RSC (Reciprocal Space Cylinder):
!
  elseif (data_struct.eq.3) then
    Qdump_name = trim(filename)//'_'//trim(data_label)//'.rsc'  ! cylinder
    open(outchan,file=Qdump_name,status='unknown')
    write(outchan,'(a100,/)') trim(adjustl(Qdump_name))
!
!  If type of data structure matches that of MPS (Monkhorse Pack Sampling):
!
  elseif (data_struct.eq.4) then
    Qdump_name = trim(filename)//'_'//trim(data_label)//'.mps'  !Monkhorse-Pack
    open(outchan,file=Qdump_name,status='unknown')
    write(outchan,'(a100,/)') trim(adjustl(Qdump_name))
!
!  If type of data structure matches that of MCS (Monte Carlo Sampling):
!
  elseif (data_struct.eq.5) then
    Qdump_name = trim(filename)//'_'//trim(data_label)//'.mcs'  ! Monte Carlo
    open(outchan,file=Qdump_name,status='unknown')
    write(outchan,'(a100,/)') trim(adjustl(Qdump_name))
!
!  If type of data structure matches that of single:
!
  elseif (data_struct.eq.6) then
    Qdump_name = trim(filename)//'_'//trim(data_label)//'.sif'  ! Single
    open(outchan,file=Qdump_name,status='unknown')
    write(outchan,'(a100,/)') trim(adjustl(Qdump_name))
    write(outchan,'(''# Reciprocal space - single dir. output from Scatter calculation :'')')
    write(outchan,'(''# Format =  individual increment, outvector (x,y,z)'')')
    write(outchan,'(''# n_qvector = '', i4)') n_qvector
!
!  Lead in with line giving dimensions of holding array for future read-in
!
    write(outchan,'(''# Please Note: all outvector components are (lattice) fractional.'',/)')
!
!  then write the data in 'data_array' into the file
!
    write(outchan,'(''# Vector is : '',f10.4,''   '',f10.4,''   '',f10.4,''   to '',f10.4,''   '',f10.4,''   '',f10.4)') &
          q_ordinate(1,1),q_ordinate(2,1),q_ordinate(3,1),q_ordinate(4,1),q_ordinate(5,1),q_ordinate(6,1)
    write(outchan,'(''# In '',i6, '' increments. '')') n_qpoint(1)
    write(outchan,'(''#-------------------------------------------------------------------------------------------------'',/)')
!  
    if (outchan .eq. 74) then
      title_string = '# Point          ' // title_string_1
      write(outchan,'(60a)') title_string
    elseif (outchan .eq. 75) then
      title_string = '# Point          ' // title_string_2
      write(outchan,'(60a)') title_string
    elseif (outchan .eq. 76) then
      title_string = '# Point          ' // title_string_3
      write(outchan,'(60a)') title_string
    endif  
    do icounter = 1, n_qpoint(1)
      do gencount = 1,3 
        curr_Q(gencount)   = Hold_Q(gencount,icounter,1,1)  
        curr_smq(gencount) = Hold_smq(gencount,icounter,1,1)
        curr_tau(gencount) = Hold_tau(gencount,icounter,1,1)
      enddo   
      if (outchan .eq. 74) then
        curr_out(1) = curr_Q(1)
        curr_out(2) = curr_Q(2)
        curr_out(3) = curr_Q(3)
      elseif (outchan .eq. 75) then
        curr_out(1) = curr_smq(1)
        curr_out(2) = curr_smq(2)
        curr_out(3) = curr_smq(3)
      elseif (outchan .eq. 76) then
        curr_out(1) = curr_tau(1)
        curr_out(2) = curr_tau(2)
        curr_out(3) = curr_tau(3)
      endif  
      write(outchan,'(i6,3f16.6)') icounter,curr_out(1),curr_out(2),curr_out(3)
    enddo
    
  endif
!
!  End of Single format
! 
  close(outchan)
  return
  end

!***********************************************************************
!*   Subroutine cscatter: Calculation Kernel for coherent scattering   *
!*   This routine calculates S(q,w) coh for a given Q, q, polarisation *
!*   vectors and scattering length b, over a range of q points passed  *
!*                                                                     *
!*       UPDATED FOR F90                                               *
!*                                                                     *
!***********************************************************************

  subroutine cscatter(qpoint,ncfoc,maxd2,eigr,eigi,lprint)
        
  use constants   
  use current 
  use frequencies
  use scatterdata
  use configurations,   only : tempcfg
  use ksample_scatter,  only : xskpt, yskpt, zskpt
  use iochannels
  use parallel,         only : ioproc
  implicit none
!
!  Passed arguments
!
  integer,     intent(in)                      :: qpoint          ! Q point being processed
  integer(i4), intent(in)                      :: maxd2           ! LHS dimension of eigr/eigi
  integer(i4), intent(in)                      :: ncfoc           ! Number of cores in the system
  logical,     intent(in)                      :: lprint          ! Print output or not?
  real(dp),    intent(in)                      :: eigr(maxd2,*)   ! Real eigenvectors              !kt assumed-size array. eigr stores a maxd2-size eigenvector one by one as called in phonon.
  real(dp),    intent(in)                      :: eigi(maxd2,*)   ! Imaginary  eigenvectors        
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: indj
  integer(i4)                                  :: status
  logical                                      :: off_gamma
  real(dp),     dimension(:),     allocatable  :: phaser
  real(dp),     dimension(:),     allocatable  :: phasei
  real(dp)                                     :: cmfact
  real(dp)                                     :: Qmodu
  real(dp)                                     :: FofQr
  real(dp)                                     :: FofQi
  real(dp)                                     :: qedotr
  real(dp)                                     :: qedoti
  real(dp)                                     :: FofQ2
  real(dp)                                     :: qcdot
  real(dp)                                     :: mytemp
  real(dp)                                     :: omsense
  real(dp)                                     :: Qsense
  real(dp)                                     :: partr
  real(dp)                                     :: parti
  real(dp)                                     :: nfac
  real(dp)                                     :: boltzfac
  real(dp)                                     :: tsense
  real(dp)                                     :: Qvect(3)
  real(dp)                                     :: tauvect(3)
  real(dp)                                     :: gamma_cut
!
!  Allocate local memory
!
  allocate(phaser(3*numat),stat=status)
  if (status/=0) call outofmemory('cscatter','phaser')
  allocate(phasei(3*numat),stat=status)
  if (status/=0) call outofmemory('cscatter','phasei')
!
!  Define constants
!
  boltzfac = 1.0_dp
  omsense = 1.0d-4
  tsense  = 1.0d-6
  Qsense  = 1.0d-8


  mytemp = max(tempcfg(ncf),1.0d-8)
  gamma_cut = 1.0d-6
  cmfact = planck*speedl/(boltz*mytemp)
!
  off_gamma = .true.
!
!  Set if flag for gamma point; set gamma intensities to zero (note that S(Q,w) is strictly undefined for gamma
!  also set Qmodu; sqrt of modulus... a test to ensure that Q is sizeable enough to bother with...
!
  Qmodu = (abs(Qvector(1,qpoint)) + abs(Qvector(2,qpoint)) + abs(Qvector(3,qpoint)))
  if (Qmodu .lt. gamma_cut) then
    off_gamma = .false.
  endif
!
!  Convert Q and tau to absolute units (Inverse Angstroms)
!
  Qvect = matmul(kv, Qvector(:, qpoint))
!        
<<<<<<< HEAD
  !tauvect(1) = tauvector(1,qpoint)*kv(1,1) + tauvector(2,qpoint)*kv(1,2) + tauvector(3,qpoint)*kv(1,3)
  !tauvect(2) = tauvector(2,qpoint)*kv(2,1) + tauvector(2,qpoint)*kv(2,2) + tauvector(3,qpoint)*kv(2,3) !kt first term should be tauvector(1,qpoint)*kv(2,1), this is probably a bug.
  !tauvect(3) = tauvector(3,qpoint)*kv(3,1) + tauvector(2,qpoint)*kv(3,2) + tauvector(3,qpoint)*kv(3,3) !kt first term should be tauvector(1,qpoint)*kv(2,1)
=======
  !kt tauvect(1) = tauvector(1,qpoint)*kv(1,1) + tauvector(2,qpoint)*kv(1,2) + tauvector(3,qpoint)*kv(1,3)
  !kt tauvect(2) = tauvector(2,qpoint)*kv(2,1) + tauvector(2,qpoint)*kv(2,2) + tauvector(3,qpoint)*kv(2,3) !kt first term should be tauvector(1,qpoint)*kv(2,1), this is probably a bug.
  !kt tauvect(3) = tauvector(3,qpoint)*kv(3,1) + tauvector(2,qpoint)*kv(3,2) + tauvector(3,qpoint)*kv(3,3) !kt first term should be tauvector(1,qpoint)*kv(2,1)
>>>>>>> 4ee3e8ac7e3f4d4e86c52737f945af3adff2f7ec
  tauvect = matmul(kv, tauvector(:, qpoint))
!
!  Calculate atomic phase factors
!
  if (off_gamma) then
    do i = 1,ncfoc
      qcdot = tauvect(1)*xclat(i) + tauvect(2)*yclat(i) + tauvect(3)*zclat(i)
      phaser(i) = cos(qcdot)
      phasei(i) = sin(qcdot)
    enddo      
!
!  Outer loop over modes
!
    do i = 1,3*ncfoc
         
      if (freq(i,qpoint) .lt. omsense) then
        freq(i,qpoint) = 0.0_dp
      endif
!
!  Inner loop over components of eigenvector
!
      FofQr = 0.0_dp
      FofQi = 0.0_dp

      indj = 0
      do j = 1,ncfoc  ! over cores
!
!  Eigenvector . Qvect for each atom DEBUG 
!
		qedotr = dot_product(Qvect, eigr(indj+1:indj+3,i))
		qedoti = dot_product(Qvect, eigi(indj+1:indj+3,i))
!
!  Inner product of Q.e(ds) with atomic phase factors
!
        partr = qedotr*phaser(j) - qedoti*phasei(j)
        parti = qedoti*phaser(j) + qedotr*phasei(j)
!
!  Scattering function kernel calculation
!
        FofQr = FofQr + partr*scatlencoh(j)
        FofQi = FofQi + parti*scatlencoh(j)
!
        indj = indj + 3
      enddo  ! over cores
        
      if (Qmodu .gt. Qsense) then
        FofQ2 = FofQr*FofQr + FofQi*FofQi
      else
        FofQ2 = 0.0_dp
      endif
!
!  Calculate <n(s) +1>, nfac
!
      if (mytemp .lt. tsense) then     
        nfac = 1.0_dp
      else
        boltzfac = exp(cmfact*freq(i,qpoint))
        nfac = boltzfac/(boltzfac - 1.0_dp)
      endif

	  !if (lscatsurfaceweight) then                     ! kt this if sentence is no meaning presently.
      !  if (freq(i,qpoint) .lt. omsense) then
      !    sofomega(i,qpoint) = 0.0_dp
      !  else
      !    sofomega(i,qpoint) = nfac*FofQ2/(3*dble(numat)*freq(i,qpoint))
      !  endif ! sflag
      !else
        if (freq(i,qpoint) .lt. omsense) then
          sofomega(i,qpoint) = 0.0_dp
        else
          sofomega(i,qpoint) = nfac*FofQ2/(3*dble(numat)*freq(i,qpoint))
!
!  Check that weighting the intensity according to angle theta is valid!  Valid for single crystal, but powder? ! kt for powder
!                                                                                                                    average, sofomega(i, qpoint) should be multiplied by sin(theta).
!          sofomega(i,qpoint) = nfac*(dtheta*dphi*dsin(thetain)*FofQ2)/(3*dble(numat)*freq(i,qpoint))
        endif          
      !endif                                             ! kt this if sentence is no meaning presently.
    enddo      ! over modes
  else 
!
!  For the case of a gamma point, for which S(Q,w) is strictly undefined, we insert zero intensity
!
    do i = 1,3*numat
      sofomega(i,qpoint) = 0.0_dp
    enddo
  endif
!  if (ldebug.and.ioproc) then
!
!  Check returned arrays from scatter
!
!    write(ioout,'(/,''Kpoint: '',i4,/)') qpoint
!    write(ioout,'(''----------'',/)')
!    do i = 1,3*numat    !over modes
!      write(ioout,'(/,''Frequency for mode '',i4,'' is '',f16.8)') j,freq(j,i)
!      write(ioout,'(/,''Atom & ords'',2x,''Real part'',7x,''Imaginary part'',/)')
!      indj = 0
!      do j = 1,numat
!        write(ioout,'(i6,'' x '',2f16.8)') j, eigr(indj+1,i), eigi(indj+1,i)
!        write(ioout,'(i6,'' y '',2f16.8)') j, eigr(indj+2,i), eigi(indj+2,i)
!        write(ioout,'(i6,'' z '',2f16.8)') j, eigr(indj+3,i), eigi(indj+3,i)
!        indj = indj + 3
!      enddo
!    enddo
!  endif 
  if (q_runtype.eq.6.and.ioproc.and.lprint) then
!
!  For q_runtype = 6, write eigenvector information to channel 69
!
    write(69,'(/,''Q point: '',i4,'' is :'',f16.8,2x,f16.8,2x,f16.8)') &
      qpoint,Qvector(1,qpoint),Qvector(2,qpoint),Qvector(3,qpoint)
    write(69,'(/,''q point: '',i4,'' is :'',f16.8,2x,f16.8,2x,f16.8)') &
      qpoint,xskpt(qpoint),yskpt(qpoint),zskpt(qpoint)
    write(69,'(/,''Tau vec: '',i4,'' is :'',f16.8,2x,f16.8,2x,f16.8)') &
      qpoint,tauvector(1,qpoint),tauvector(2,qpoint),tauvector(3,qpoint)
    write(69,'(''----------------------------------------------------------------------'',/)')
    do i = 1,3*ncfoc
      write(69,'(/,''Frequency for mode '',i4,'' is '',f16.8)') i,freq(i,qpoint)
      write(69,'(''S(Q,w) for mode '',i4,'' is '',f25.8,/)') i,sofomega(i,qpoint)
      write(69,'(/,''Atom & ords'',2x,''Real part'',7x,''Imaginary part'',/)')
      indj = 0
      do j = 1,ncfoc
        write(69,'(i6,'' x '',2f16.8)') j, eigr(indj+1,i), eigi(indj+1,i)
        write(69,'(i6,'' y '',2f16.8)') j, eigr(indj+2,i), eigi(indj+2,i)
        write(69,'(i6,'' z '',2f16.8)') j, eigr(indj+3,i), eigi(indj+3,i)
        indj = indj + 3
      enddo 
    enddo 
  endif
!
!  Deallocate local memory
!
  deallocate(phasei,stat=status)
  if (status/=0) call deallocate_error('cscatter','phasei')
  deallocate(phaser,stat=status)
  if (status/=0) call deallocate_error('cscatter','phaser')
!
  return
  end

  subroutine sqomegadifference(sqwarray, nsqwarraylhs, nsqwarrayrhs, q_vectors_in, fsumsq)   
      
  use scatterdata
      
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                       :: nsqwarraylhs
  integer(i4), intent(in)                       :: nsqwarrayrhs
  integer(i4), intent(in)                       :: q_vectors_in
  real(dp),    intent(in)                       :: sqwarray(nsqwarraylhs,nsqwarrayrhs)
  real(dp),    intent(out)                      :: fsumsq
!
!  Local variables
!
  integer(i4)                                   :: i
  integer(i4)                                   :: j
  real(dp)                                      :: rnormcalc
!
!  Check dimensions
!
  if (q_vectors_in.ne.nq_step_fit) then
    call outerror('S(Q,omega) array size does not match in fitting',0_i4)
    call stopnow('fitword')
  endif
  if (nw_step.ne.nw_step_fit) then
    call outerror('S(Q,omega) array size does not match in fitting',0_i4)
    call stopnow('fitword')
  endif
!
!  Compute normalisation factors
!
  rnormcalc = 0.0_dp
  do i = 1,q_vectors_in
    do j = 1,nw_step
      rnormcalc = rnormcalc + (sqwarray(j,i))**2
    enddo
  enddo
  if (rnormcalc.gt.1.0d-12) then
    rnormcalc = 1.0_dp/sqrt(rnormcalc)
  endif
!
  do i = 1,q_vectors_in
    do j = 1,nw_step
      fsumsq = fsumsq + (rnormcalc*sqwarray(j,i) - sofomega_fit(j,i))**2
    enddo
  enddo
!
  fsumsq = sqrt(fsumsq)
!
  end

  subroutine sumhist(outarray)
!
!  Sum up histogrammed data accumulated over multiple nodes
!
!  Julian Gale, Curtin University, September 2010
!
  use scatterdata
  use parallel,    only : nprocs
  use times,       only : tsum

  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout)                   :: outarray(nw_step,*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: status
  real(dp)                                     :: g_cpu_time
  real(dp),     dimension(:,:),   allocatable  :: sumarray
  real(dp)                                     :: tsum0
      
  if (nprocs.gt.0) then
!
!  Allocate workspace array
!
    allocate(sumarray(nw_step,nq_step),stat=status)
    if (status/=0) call outofmemory('sumhist','sumarray')
!
!  Perform parallel reduction of outarray
!
    tsum0 = g_cpu_time()
    call sumall(outarray(1,1),sumarray(1,1),nw_step*nq_step,"sumhist","outarray")
!
!  Copy data back from sum array to endarray
!
    do i = 1,nq_step
      do j = 1,nw_step
        outarray(j,i) = sumarray(j,i)
      enddo
    enddo      
    tsum = tsum + g_cpu_time() - tsum0
!
!  Deallocate workspace array
!
    deallocate(sumarray,stat=status)
    if (status/=0) call deallocate_error('sumhist','sumarray')
  endif
!
  return
  end
