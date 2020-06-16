  subroutine setkeyword
!
!  Sets general flags based on keywords
!
!   2/04 Created from getkeyword
!   7/05 lSandM added
!   7/05 Setting of linten flag improved to remove ambiguity
!   9/06 libdump keyword added
!  11/06 NEB modifications added
!   1/07 Gasteiger charges added
!   3/07 lPrintEAM keyword added
!   3/07 lPrintTwo keyword added
!   3/07 More robust checking of keywords added by ensuring
!        string comes at the start of the word
!   3/07 lpreserveQ added
!   4/07 Conjugate gradients set if running in parallel
!   5/07 Mean KE option added
!   5/07 qbond keyword added
!   7/07 lmeta added
!  10/08 COSMO/COSMIC keywords merged in 
!  10/08 Error in logic for lmodco setting corrected
!  10/08 Error in logic for other keywords also corrected
!  11/08 lPrintFour keyword added
!  11/08 lPrintThree keyword added
!   2/09 lconj only set to be true for > one processor if not LM-BFGS
!   6/09 Module name changed from three to m_three
!   6/09 PDF keywords added
!   7/09 Symmetry turned off for derivatives in reaxFF case
!  12/09 pregionforce keyword added
!   4/10 qtpie keyword added
!   6/10 Hopping keyword added
!   8/10 lconvert, lphase, lcutbelow keywords removed
!   8/10 Keyword settings adjusted for PDF case
!   8/10 lfix1atom added
!  10/10 Symmetry turned off for EDIP potentials
!  12/10 Hiding of shells added
!   1/11 Force minimisation option added
!   2/11 Keyword added to turn off ReaxFF charges
!   3/11 lstressout added
!   9/11 Metadynamics internal code replaced with Plumed
!   9/11 Madelung correction added
!  11/11 eregion keyword added
!   5/12 Atomic stress keyword added
!   6/12 Thermal conductivity added
!   7/12 Oldinten keyword added
!   9/12 Pacha added
!   9/12 Eckart transformation added 
!  10/12 Option added to allow two atoms co-exist on a site with an
!        occupancy greater than one.
!  12/12 site_energy keyword added
!   7/13 Modified so that multiple calls do not overwrite flags set by options
!   8/13 linclude_imaginary added
!   8/13 lscatter added back
!   8/13 lfindsym added
!   8/13 Raman keyword added
!  10/13 lsopt keyword added
!   3/14 Harmonic relaxation keyword added
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
  use bondorderdata,  only : nboQ
  use control
  use cosmo,          only : lcosmic, lsegsmooth
  use distances,      only : lStoreVectors
  use eam,            only : lPrintEAM
  use element
  use fitting
  use four,           only : lPrintFour
  use library,        only : llibsymdump
  use m_pdfneutron,   only : lmakeeigarray, lcoreinfo, lpdf
  use m_pdfneutron,   only : lnowidth, lpartial, lfreqcut, lkeepcut, lnoksym
  use m_pdfneutron,   only : lpdfout
  use m_three,        only : lPrintThree
  use mdlogic
  use molecule
  use neb,            only : lnebclimbingimage, lnebdoublynudged
  use optimisation
  use parallel,       only : nprocs
  use phonout,        only : linclude_imaginary
  use symmetry
  use synchro,        only : lfixtangent
  use two,            only : lPrintTwo
  use wolfcosmo,      only : lPureCoulomb0D
  implicit none
!****************************
!  Set flags from keywords  *
!****************************
!
!  Keywords with default = .false.
!
  if (.not.lallowgt1) lallowgt1 = (index(keyword,' allo').ne.0.or.index(keyword,'allo').eq.1)
  if (.not.langle) langle = (index(keyword,' angl').ne.0.or.index(keyword,'angl').eq.1)
  if (.not.lanneal) lanneal = (index(keyword,' anne').ne.0.or.index(keyword,'anne').eq.1)
  if (.not.latomicstress) latomicstress = (index(keyword,' atom').ne.0.or.index(keyword,'atom').eq.1)
  if (.not.laver) laver = (index(keyword,' aver').ne.0.or.index(keyword,'aver').eq.1)
  if (.not.lbond) lbond = (index(keyword,' bond').ne.0.or.index(keyword,'bond').eq.1)
  if (.not.lbroad) lbroad = (index(keyword,' broa').ne.0.or.index(keyword,'broa').eq.1)
  if (.not.lbulknoopt) lbulknoopt = (index(keyword,' bulk').ne.0.or.index(keyword,'bulk').eq.1)
  if (.not.lc6) lc6 = ((index(keyword,' c6').ne.0.or.index(keyword,'c6').eq.1).and.index(keyword,'noel').eq.0)
  if (.not.lcello) lcello = (index(keyword,' cell').ne.0.or.index(keyword,'cell').eq.1)
  if (.not.lcomp) lcomp = (index(keyword,' comp').ne.0.or.index(keyword,'comp').eq.1)
  if (.not.lconj) lconj = (index(keyword,' conj').ne.0.or.index(keyword,'conj').eq.1)
  if (.not.lconp) lconp = (index(keyword,' conp').ne.0.or.index(keyword,'conp').eq.1)
  if (.not.lconv) lconv = (index(keyword,' conv').ne.0.or.index(keyword,'conv').eq.1)
  if (.not.lcosmic) lcosmic = (index(keyword,'cosmi').ne.0)
  if (.not.lcosmo) lcosmo = (index(keyword,'cosmo').ne.0)
  if (.not.ldcharge) ldcharge = (index(keyword,' dcha').ne.0.or.index(keyword,'dcha').eq.1)
  if (.not.ldebug) ldebug = (index(keyword,' debu').ne.0.or.index(keyword,'debu').eq.1)
  if (.not.ldefect) ldefect = (index(keyword,' defe').ne.0.or.index(keyword,'defe').eq.1)
  if (.not.ldfp) ldfp = ((index(keyword,' dfp').ne.0.or.index(keyword,'dfp').eq.1).and.index(keyword,'bfgs').eq.0)
  if (.not.ldipole) ldipole = (index(keyword,' dipo').ne.0.or.index(keyword,'dipo').eq.1)
  if (.not.lStoreVectors) lStoreVectors = (index(keyword,' stor').ne.0.or.index(keyword,'stor').eq.1)
  if (.not.ldist) ldist = (index(keyword,' dist').ne.0.or.index(keyword,'dist').eq.1)
  if (.not.leckart) leckart = (index(keyword,' eck').ne.0.or.index(keyword,'eck').eq.1)
  if (.not.leem) leem = (index(keyword,' eem').ne.0.or.index(keyword,'eem').eq.1)
  if (.not.lefg) lefg = (index(keyword,' efg').ne.0.or.index(keyword,'efg').eq.1)
  if (.not.leigen) leigen = (index(keyword,' eige').ne.0.or.index(keyword,'eige').eq.1)
  if (.not.leregion) leregion = (index(keyword,' ereg').ne.0.or.index(keyword,'ereg').eq.1)
  if (.not.lfbfgs) lfbfgs = (index(keyword,' fbfg').ne.0.or.index(keyword,'fbfg').eq.1)
  if (.not.lfindsym) lfindsym = (index(keyword,' find').ne.0.or.index(keyword,'find').eq.1)
  if (.not.lfit) lfit = (index(keyword,' fit').ne.0.or.index(keyword,'fit').eq.1)
  if (.not.lforcemin) lforcemin = (index(keyword,' forc').ne.0.or.index(keyword,'forc').eq.1)
  if (.not.lfree) lfree = (index(keyword,' free').ne.0.or.index(keyword,'free').eq.1)
  if (.not.lfreq) lfreq = (index(keyword,' freq').ne.0.or.index(keyword,'freq').eq.1)
  if (.not.lfreqout) lfreqout = (index(keyword,' nofr').eq.0.and.index(keyword,'nofr').ne.1)
  if (.not.lfixtangent) lfixtangent = (index(keyword,' tfix').ne.0.and.index(keyword,'tfix').eq.1)
  if (.not.lga) lga = (index(keyword,' gene').ne.0.or.index(keyword,'gene').eq.1)
  if (.not.lgasteiger) lgasteiger = (index(keyword,' gast').ne.0.or.index(keyword,'gast').eq.1)
  if (.not.lgrad) lgrad = (index(keyword,' grad').ne.0.or.index(keyword,'grad').eq.1)
  if (.not.lharmrelax) lharmrelax = (index(keyword,' ehar').ne.0.or.index(keyword,'ehar').eq.1)
  if (.not.lhex) lhex = (index(keyword,' hex').ne.0.or.index(keyword,'hex').eq.1)
  if (.not.lhideshells) lhideshells = (index(keyword,' hide').ne.0.or.index(keyword,'hide').eq.1)
  if (.not.linclude_imaginary) linclude_imaginary = (index(keyword,' incl').ne.0.or.index(keyword,'incl').eq.1)
  if (.not.linten) linten = (index(keyword,' inte').ne.0.or.index(keyword,'inte').eq.1)
  if (.not.lkfull) lkfull = (index(keyword,' kful').ne.0.or.index(keyword,'kful').eq.1)
  if (.not.llbfgs) llbfgs = (index(keyword,' lbfg').ne.0.or.index(keyword,'lbfg').eq.1)
  if (.not.llibsymdump) llibsymdump = (index(keyword,' libd').ne.0.or.index(keyword,'libd').eq.1)
  if (.not.lmadelung) lmadelung = (index(keyword,' made').ne.0.or.index(keyword,'made').eq.1)
  if (.not.lmarvreg2) lmarvreg2 = (index(keyword,' marv').ne.0.or.index(keyword,'marv').eq.1)
  if (.not.lmc) lmc = (index(keyword,' mont').ne.0.or.index(keyword,'mont').eq.1)
  if (.not.lmd) lmd = (index(keyword,' md').ne.0.or.index(keyword,'md').eq.1)
  if (.not.lmeanke) lmeanke = (index(keyword,' mean').ne.0.or.index(keyword,'mean').eq.1)
  if (.not.lminimage) lminimage = (index(keyword,' mini').ne.0.or.index(keyword,'mini').eq.1)
  if (.not.lmol) lmol = (index(keyword,' mol').ne.0.or.index(keyword,'mol').eq.1)
  if (.not.lmolq) lmolq = (index(keyword,' molq').ne.0.or.index(keyword,'molq').eq.1)
  if (.not.lmolmec) lmolmec = (index(keyword,' molm').ne.0.or.index(keyword,'molm').eq.1)
  if (.not.lmolfix) lmolfix = (index(keyword,' fix').ne.0.or.index(keyword,'fix').eq.1)
  if (.not.lneb) lneb = (index(keyword,' neb').ne.0.or.index(keyword,'neb').eq.1)
  if (.not.lnebclimbingimage) lnebclimbingimage = (index(keyword,' cineb').ne.0.or.index(keyword,'cineb').eq.1)
  if (.not.lnoautobond) lnoautobond = (index(keyword,' noau').ne.0.or.index(keyword,'noau').eq.1)
  if (.not.lnoenergy) lnoenergy = (index(keyword,' noen').ne.0.or.index(keyword,'noen').eq.1)
  if (.not.lnoflags) lnoflags = (index(keyword,' nofl').ne.0.or.index(keyword,'nofl').eq.1)
  if (.not.lnoreal) lnoreal = (index(keyword,' noreal').ne.0.or.index(keyword,'noreal').eq.1)
  if (.not.lnorecip) lnorecip = (index(keyword,' noreci').ne.0.or.index(keyword,'noreci').eq.1)
  if (.not.lnumdiag) lnumdiag = (index(keyword,' numd').ne.0.or.index(keyword,'numd').eq.1)
  if (.not.loldinten) loldinten = (index(keyword,' oldi').ne.0.or.index(keyword,'oldi').eq.1)
  if (.not.lopt) lopt = (index(keyword,' opti').ne.0.or.index(keyword,'opti').eq.1)
  if (.not.loptcellpar) loptcellpar = (index(keyword,' ocel').ne.0.or.index(keyword,'ocel').eq.1)
  if (.not.lpacha) lpacha = (index(keyword,' pac').ne.0.or.index(keyword,'pac').eq.1)
  if (.not.lphon) lphon = (index(keyword,' phon').ne.0.or.index(keyword,'phon').eq.1)
  if (.not.lposidef) lposidef = (index(keyword,' posi').ne.0.or.index(keyword,'posi').eq.1)
  if (.not.lpot) lpot = (index(keyword,' pot').ne.0.or.index(keyword,'pot').eq.1)
  if (.not.lpredict) lpredict = (index(keyword,' pred').ne.0.or.index(keyword,'pred').eq.1)
  if (.not.lpreserveQ) lpreserveQ = (index(keyword,' pres').ne.0.or.index(keyword,'pres').eq.1)
  if (.not.lPrintEAM) lPrintEAM = (index(keyword,' prt_eam').ne.0.or.index(keyword,'ptr_eam').eq.1)
  if (.not.lPrintFour) lPrintFour = (index(keyword,' prt_fo').ne.0.or.index(keyword,'ptr_fo').eq.1)
  if (.not.lPrintThree) lPrintThree = (index(keyword,' prt_th').ne.0.or.index(keyword,'ptr_th').eq.1)
  if (.not.lPrintTwo) lPrintTwo = (index(keyword,' prt_two').ne.0.or.index(keyword,'ptr_two').eq.1)
  if (.not.lprop) lprop = (index(keyword,' prop').ne.0.or.index(keyword,'prop').eq.1)
  if (.not.lPureCoulomb0D) lPureCoulomb0D = (index(keyword,'pure').ne.0)
  if (.not.lqbond) lqbond = (index(keyword,' qbon').ne.0.or.index(keyword,'qbon').eq.1)
  if (.not.lqeq) lqeq = (index(keyword,' qeq').ne.0.or.index(keyword,'qeq').eq.1)
  if (.not.lqtpie) lqtpie = (index(keyword,' qtp').ne.0.or.index(keyword,'qtp').eq.1)
  if (.not.lquicksearch) lquicksearch = (index(keyword,' quic').ne.0.or.index(keyword,'quic').eq.1)
  if (.not.lraman) lraman = (index(keyword,' rama').ne.0.or.index(keyword,'rama').eq.1)
  if (.not.lrelax) lrelax = ((index(keyword,' rela').ne.0.or.index(keyword,'rela').eq.1).and.lfit)
  if (.not.lregionforce) lregionforce = (index(keyword,' preg').ne.0.or.index(keyword,'preg').eq.1)
  if (.not.lrest) lrest = (index(keyword,' rest').ne.0.or.index(keyword,'rest').eq.1)
  if (.not.lrfo) lrfo = (index(keyword,' rfo').ne.0.or.index(keyword,'rfo').eq.1.or.index(keyword,'tran').ne.0)
  if (.not.lSandM) lSandM = (index(keyword,' sm').ne.0.or.index(keyword,'sm ').eq.1)
  if (.not.lsave) lsave = (index(keyword,' save').ne.0.or.index(keyword,'save').eq.1)
  if (.not.lscatter) lscatter = (index(keyword,' scat').ne.0.or.index(keyword,'scat').eq.1)
  if (.not.lshello) lshello = (index(keyword,' shel').ne.0.or.index(keyword,'shel').eq.1)
  if (.not.lsiteenergy) lsiteenergy = (index(keyword,' site').ne.0.or.index(keyword,'site').eq.1)
  if (.not.lsopt) lsopt = (index(keyword,' sopt').ne.0.or.index(keyword,'sopt').eq.1)
  if (.not.lspatial) lspatial = (index(keyword,' spat').ne.0.or.index(keyword,'spat').eq.1)
  if (.not.lstaticfirst) lstaticfirst = (index(keyword,' stat').ne.0.or.index(keyword,'stat').eq.1)
  if (.not.lstressout) lstressout = (index(keyword,' stre').ne.0.or.index(keyword,'stre').eq.1)
  if (.not.lsymoff) lsymoff = (index(keyword,' symoff').ne.0.or.index(keyword,'symoff').eq.1)
  if (.not.lthermal) lthermal = (index(keyword,' ther').ne.0.or.index(keyword,'ther').eq.1)
  if (.not.ltors) ltors = (index(keyword,' tors').ne.0.or.index(keyword,'tors').eq.1)
  if (.not.ltran) ltran = (index(keyword,' tran').ne.0.or.index(keyword,'tran').eq.1)
  if (.not.lunit) lunit = (index(keyword,' unit').ne.0.or.index(keyword,'unit').eq.1)
  if (.not.lzsisa) lzsisa = (index(keyword,' zsis').ne.0.or.index(keyword,'zsis').eq.1)
!
!  Keywords with default = .true.
!
  if (ldsym) ldsym = (index(keyword,' nodsym').eq.0.and.index(keyword,'nodsym').ne.1)
  if (lfix1atom) lfix1atom = (index(keyword,' unfi').eq.0.and.index(keyword,'unfi').ne.1)
  if (lDoElectrostatics) lDoElectrostatics = (index(keyword,' noel').eq.0.and.index(keyword,'noel').ne.1)
  if (lmodco) lmodco = (index(keyword,' nomod').eq.0.and.index(keyword,'nomod').ne.1)
  if (lnebdoublynudged) lnebdoublynudged = (index(keyword,' nodn').eq.0.and.index(keyword,'nodn').ne.1)
  if (lreaxFFQ) lreaxFFQ = (index(keyword,' norx').eq.0.and.index(keyword,'norx').ne.1)
  if (lSandMnozz) lSandMnozz = (index(keyword,' smzz').eq.0.and.index(keyword,'smzz').ne.1)
  if (lsegsmooth) lsegsmooth = (index(keyword,' nosmo').eq.0.and.index(keyword,'nosmo').ne.1)
  if (lsym) lsym = (index(keyword,' nosy').eq.0.and.index(keyword,'nosy').ne.1)
  if (lsymdok) lsymdok = (index(keyword,' nosd').eq.0.and.index(keyword,'nosd').ne.1)
!
!  Neutron Keywords (ers29)
!
  if (.not.lmakeeigarray) lmakeeigarray = (index(keyword,' make').ne.0.or.index(keyword,'make').eq.1)
  if (.not.lcoreinfo) lcoreinfo = (index(keyword,' core').ne.0.or.index(keyword,'core').eq.1)
  if (.not.lnoksym) lnoksym = (index(keyword,' noks').ne.0.or.index(keyword,'noks').eq.1)
  if (.not.lpdf) lpdf = (index(keyword,' pdf').ne.0.or.index(keyword,'pdf').eq.1)
  if (.not.lfreqcut) lfreqcut = (index(keyword,' pdfc').ne.0.or.index(keyword,'pdfc').eq.1)
  if (.not.lkeepcut) lkeepcut = (index(keyword,' pdfk').ne.0.or.index(keyword,'pdfk').eq.1)
  if (.not.lnowidth) lnowidth = (index(keyword,' nowi').ne.0.or.index(keyword,'nowi').eq.1)
  if (.not.lpartial) lpartial = (index(keyword,' nopa').eq.0.or.index(keyword,'nopa').eq.1)
!*******************************
!  Handle keyword dependances  *
!*******************************
!
!  lflags depends on other keywords
!
  if (lnoflags) then
    lflags = .false.
  else
    if (lopt.or.lgrad.or.lharmrelax.or.lfit.or.lrfo.or.lmc.or.lmd.or.lneb) then
      if ((.not.lconp).and.(.not.lconv).and.(.not.lcello)) lflags = .true.
    endif
  endif
!
!  PDF dependances: (ers29)
!
  if (lkeepcut) lfreqcut = .true.
  if (lfreqcut.or.lkeepcut) lpdf = .true.
  if (lpdf) then
    lmakeeigarray = .true.
    lsym = .false.
    if (.not.(index(keyword,' full').ne.0.or.index(keyword,'full').eq.1)) then
      write(keyword,*) trim(adjustl(keyword)), " full"
    endif
  endif
  if (lcoreinfo) lphon = .true.
  if (lmakeeigarray) then
    lphon = .true.
    leigen = .true.
    lnoksym = .true.
    lpdfout = .true.
  endif
!
!  If CI-NEB is requested then turn on NEB too
!
  if (lnebclimbingimage) lneb = .true.
!
!  If defect calc, then we need properties
!
  if (ldefect.and..not.lrest) lprop = .true.
!
!  If prop calc, then strain derivatives are needed
!
  if (lprop) then
    lstr = .true.
  endif
!
!  If prop calc, set flag for Born effective charges to be true
!
  if (lprop) then
    lborn = .true.
  endif
!
!  If atomic stresses are requested then strain derivatives are needed
!
  if (latomicstress) then
    lstr = .true.
  endif
!
!  If symmetry is to be lowered from imaginary modes
!  then we need nosym and phonon calc
!
  if (index(keyword,'lowe').ne.0) then
    lsym = .false.
    lphon = .true.
  endif
!
!  If IR intensities or eigenvectors requested, then must be a phonon calc
!
  if (linten.or.leigen) lphon = .true.
!
!  If Raman susceptibility requested, then properties must be computed
!
  if (lraman) lprop = .true.
!
!  If this is a thermal conductivity calculation then it also must be a phonon calc
!
  if (lthermal) lphon = .true.
!
!  If transition state calc, then lopt and lrfo must be true
!
  if (ltran) then
    morder = 1
    lopt = .true.
    lrfo = .true.
  endif
!
!  If number of processors is greater than 1 set conjugate gradients as optimiser
!
  if (nprocs.gt.1.and..not.llbfgs) lconj = .true.
!
!  Change default update for unit hessian
!
  if (lunit.and.lfit) then
    nfupdate = 100
  endif
!
!  If QEq or SM or Pacha, then leem must be true
!
  if (lqeq.or.lSandM.or.lqtpie.or.lpacha) leem = .true.
!
!  If variable charge then there are no third derivatives
!  but we do need charge second derivatives
!
  if (leem) then
    lnoanald3 = .true.
    lDoQDeriv1 = .false.
    lDoQDeriv2 = .true.
  endif
!
!  If bond order charge potentials or ReaxFF are present turn off symmetry for derivatives
!
  if (nboQ.gt.0.or.lreaxFF.or.lEDIP) then
    lsymdok = .false.
  endif
!
!  If COSMIC then set cosmo to be true
!
  if (lcosmic) lcosmo = .true.
!
!  If COSMO then there are no third derivatives and symmetry
!  must not be enabled for derivatives.
!
  if (lcosmo) then
    lnoanald3 = .true.
    lsymdok = .false.
  endif
!
!  If option to store vectors has been set then reset maxat arrays
!
  if (lStoreVectors) then
    call changemaxat
  endif
!
  return
  end
