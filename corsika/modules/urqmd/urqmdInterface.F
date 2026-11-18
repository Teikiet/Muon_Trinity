c 18.11.2011 Link routines between UrQMD 1.3 and CONEX.
c author T. Pierog based on CORSIKA and EPOS link to UrQMD

c adapted by M. Reininghaus for linking UrQMD to CORSIKA 8 (Apr 2019)

#ifdef __STD__
#define __GHEISHA__
#define __QGSJET__
#define __ANALYSIS__
#endif


c-----------------------------------------------------------------------
      subroutine IniUrQMDC8
c-----------------------------------------------------------------------
c Primary initialization for UrQMD 1.31
c-----------------------------------------------------------------------
      implicit none
c CONEX includes
c~ #include "conex.h"
c~ #include "conex.incnex"
#ifndef __CXCORSIKA__
c~       character*500 furqdat
c~       integer ifurqdat, nfurqdat
c~       common/urqfname/  furqdat, ifurqdat, nfurqdat

      include 'boxinc.f'
      include 'inputs.f'
      include 'options.f'
      
c commons from coms.f
      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl

c local
      INTEGER          i,io,ia,ie,id
      CHARACTER        CTPStrg(numctp)*60, CTOStrng(numcto)*60
      integer mxie,mxid,mxia
      parameter (mxie=41,mxid=10,mxia=3)
      character adum
      double precision sig_u1,ekdummy
      integer iamaxu,idmaxu,iemaxu
c~       common /cxs_u1/ sig_u1(mxie,mxid,mxia),iamaxu,idmaxu,iemaxu
c~       double precision xs(3),bim(3)
c~ c M.R.: bim added to cxs_u2
c~       common /cxs_u2/ xs,bim
      integer iudebug
c~       data bim/6.d0,6.d0,7.d0/
      integer init
      data init/0/
      SAVE

      if(init.ge.1)return
      init=init+1
#ifdef __CXDEBUG__
      call utisx1('iniurqmd  ',4)
      write(*,'(a)')'initialize URQMD ...'
#endif

C-----------------------------------------------------------------------

c~       IF ( isx.ge.2 ) THEN
c~         IUDEBUG = isx-1
c~       ELSE
c~         IUDEBUG = 0
c~       ENDIF

      WRITE (*,*)
     $   '############################################################'
      WRITE (*,*)
     $   '##                                                        ##'
      WRITE (*,*)
     $   '##     UrQMD 1.3.1  University of Frankfurt               ##'
      WRITE (*,*)
     $   '##                  urqmd@th.physik.uni-frankfurt.de      ##'
      WRITE (*,*)
     $   '##                                                        ##'
      WRITE (*,*)
     $   '############################################################'
      WRITE (*,*)
     $   '##                                                        ##'
      WRITE (*,*)
     $   '##     please cite when using this model:                 ##'
      WRITE (*,*)
     $   '##     S.A.Bass et al. Prog.Part.Nucl.Phys. 41 (1998) 225 ##'
      WRITE (*,*)
     $   '##     M.Bleicher et al. J.Phys. G25  (1999) 1859         ##'
      WRITE (*,*)
     $   '##                                                        ##'
      WRITE (*,*)
     $   '############################################################'

C  SET THE 'LARGE' CROSS-SECTIONS FOR ALL 3 TARGET ELEMENTS
c~       DO  I = 1, 3
c~         XS(I) = 10.D0 * PI * BIM(I)**2
c~       ENDDO

C  SET NMAX TO DEFAULT VALUE
      call set0
      call params_urqmd

C  THIS IS THE SUBSTITUE FOR THE URQMD INPUT ROUTINE
C  INITIALIZE COUNTERS
      boxflag = 0
      mbflag  = 0
      edens   = 0.d0
      para    = 0
      solid   = 0
      mbox    = 0
      io      = 0

C  THE FOLLOWING FLAGS CHECK, WHETHER ALL NECESSARY INPUT IS GIVEN
C  PROJECTILE
      prspflg = 0
C  TARGET
      trspflg = 0
C
      srtflag = 0
      firstev = 0
C  EXCITATION FUNCTION
      nsrt    = 1
      npb     = 1
      efuncflag = 0
C  DEFAULT NUMBER OF EVENTS
      nevents = 1
C  DEFAULT NUMBER OF TIMESTEPS
      nsteps  = 1000

C  SKIP CONDITIONS ON UNIT 13, 14, 15, 16 & 18
C  SUPPRESS ALL OUTPUT
      bf13 = .true.
      bf14 = .true.
      bf15 = .true.
      bf16 = .true.
      bf18 = .true.
      bf19 = .true.
      bf20 = .true.
C  SET DEBUG OUTPUT DEPENDING ON CHOSEN DEBUG LEVEL
C  SET THE OUTPUT OF UNITS 13, 14, 15 TO THE DEBUG OUTPUT UNIT
c~       IF     ( IUDEBUG .EQ. 1 ) THEN
c~         bf13 = .true.
c~         bf14 = .false.
c~         call uounit(14,IFCK)
c~         bf15 = .true.
c~       ELSEIF ( IUDEBUG .EQ. 2 ) THEN
c~         bf13 = .false.
c~         call uounit(13,IFCK)
c~         bf14 = .true.
c~         bf15 = .true.
c~       ELSEIF ( IUDEBUG .GT. 2 ) THEN
c~         bf13 = .true.
c~         bf14 = .true.
c~         bf15 = .false.
c~         call uounit(15,IFCK)
c~       ENDIF
      do  i = 1, numcto
         CTOdc(i) = '  '
      enddo
      do  i = 1, numctp
         CTPdc(i) = '  '
      enddo
      do  i = 1, maxstables
         stabvec(i) = 0
      enddo
      nstable = 0

C  DEFAULT SETTINGS FOR CTParam AND CTOption
C  DEFAULT SETTINGS FOR CTParam
      CTParam(1)=1.d0
      CTPStrg(1)='scaling factor for decay-width'
      CTParam(2)=0.52d0
      CTPStrg(2)='used for minimal stringmass & el/inel cut in makestr'
      CTParam(3)=2.d0
      CTPStrg(3)='velocity exponent for modified AQM'
      CTParam(4)=0.3d0
      CTPStrg(4)='transverse pion mass, used in make22 & strexct'
      CTParam(5)=0.d0
      CTPStrg(5)='probabil. for quark rearrangement in cluster'
      CTParam(6)=0.37d0
      CTPstrg(6)='strangeness probability'
      CTParam(7)=0.d0
      CTPStrg(7)='charm probability (not yet implemented in UQMD)'
      CTParam(8)=0.093d0
      CTPStrg(8)='probability to create a diquark'
      CTParam(9)=0.35d0
      CTPStrg(9)='kinetic energy cut off for last string break'
      CTParam(10)=0.25d0
      CTPStrg(10)='min. kinetic energy for hadron in string'
      CTParam(11)=0.d0
      CTPStrg(11)='fraction of non groundstate resonances'
      CTParam(12)=.5d0
      CTPStrg(12)='probability for rho 770 in String'
      CTParam(13)=.27d0
      CTPStrg(13)='probability for rho 1450 (rest->rho1700)'
      CTParam(14)=.49d0
      CTPStrg(14)='probability for omega 782'
      CTParam(15)=.27d0
      CTPStrg(15)='probability for omega 1420(rest->om1600)'
      CTParam(16)=1.0d0
      CTPStrg(16)='mass cut betw. rho770 and rho 1450'
      CTParam(17)=1.6d0
      CTPSTRG(17)='mass cut betw. rho1450 and rho1700'
      CTParam(18)=.85d0
      CTPStrg(18)='mass cut betw. om 782 and om1420'
      CTParam(19)=1.55d0
      CTPStrg(19)='mass cut betw. om1420 and om1600'
      CTParam(20)=0.0d0
      CTPStrg(20)=' distance for second projectile'
      CTParam(21)=0.0d0
      CTPStrg(21)=' deformation parameter'

      CTParam(25)=.9d0
      CTPStrg(25)=' probability for diquark not to break'
      CTParam(26)=50.d0
      CTPStrg(26)=' maximum trials to get string masses'
      CTParam(27)=1.d0
      CTPStrg(27)=' scaling factor for xmin in string excitation'
      CTParam(28)=1.d0
      CTPStrg(28)=' scaling factor for transverse fermi motion'
      CTParam(29)=0.4d0
      CTPStrg(29)=' single strange di-quark suppression factor '
      CTParam(30)=1.5d0
      CTPStrg(30)=' radius offset for initialization  '
      CTParam(31)=1.6d0
      CTPStrg(31)=' sigma of gaussian for tranverse momentum tranfer '
      CTParam(32)=0.d0
      CTPStrg(32)=' alpha-1 for valence quark distribution  '
      CTParam(33)=2.5d0
      CTPStrg(33)=' betav for valence quark distribution  (DPM)'
      CTParam(34)=0.1d0
      CTPStrg(34)=' minimal x multiplied with ecm  '
      CTParam(35)=3.0d0
      CTPStrg(35)=' offset for cut for the FSM '
      CTParam(36)=0.275d0
      CTPStrg(36)=' fragmentation function parameter a  '
      CTParam(37)=0.42d0
      CTPStrg(37)=' fragmentation function parameter b  '
      CTParam(38)=1.08d0
      CTPStrg(38)=' diquark pt scaling factor '
      CTParam(39)=0.8d0
      CTPStrg(39)=' strange quark pt scaling factor '
      CTParam(40)=0.5d0
      CTPStrg(40)=' betas-1 for valence quark distribution (LEM)'
      CTParam(41)=0.d0
      CTPStrg(41)=' distance of initialization'
      CTParam(42)=0.55d0
      CTPStrg(42)=' width of gaussian -> pt in string-fragmentation '
      CTParam(43)=5.d0
      CTPStrg(43)=' maximum kinetic energy in mesonic clustr '
      CTParam(44)=0.8d0
      CTPStrg(44)=' prob. of double vs. single excitation for AQM inel.'
      CTParam(45)=0.5d0
      CTPStrg(45)=' offset for minimal mass generation of strings'
      CTParam(46)=800000.d0
      CTPStrg(46)=' maximal number of rejections for initialization'
      CTParam(47)=1.0d0
      CTPStrg(47)=' field feynman fragmentation funct. param. a'
      CTParam(48)=2.0d0
      CTPStrg(48)=' field feynman fragmentation funct. param. b'

      CTParam(50)=1.d0
      CTPStrg(50)=' enhancement factor for 0- mesons'
      CTParam(51)=1.d0
      CTPStrg(51)=' enhancement factor for 1- mesons'
      CTParam(52)=1.d0
      CTPStrg(52)=' enhancement factor for 0+ mesons'
      CTParam(53)=1.d0
      CTPStrg(53)=' enhancement factor for 1+ mesons'
      CTParam(54)=1.d0
      CTPStrg(54)=' enhancement factor for 2+ mesons'
      CTParam(55)=1.d0
      CTPStrg(55)=' enhancement factor for 1+-mesons'
      CTParam(56)=1.d0
      CTPStrg(56)=' enhancement factor for 1-*mesons'
      CTParam(57)=1.d0
      CTPStrg(57)=' enhancement factor for 1-*mesons'
      CTParam(58)=1.d0
      CTPStrg(58)=' scaling factor for DP time-delay'

C  DEFAULT SETTINGS FOR CTOption
      CTOption(1)=1                  ! hjd1
      CTOStrng(1)=' resonance widths are mass dependent '
      CTOption(2)=0
      CTOStrng(2)=' conservation of scattering plane'
      CTOption(3)=0
      CTOStrng(3)=' use modified detailed balance'
      CTOption(4)=0
      CTOStrng(4)=' no initial conf. output '
      CTOption(5)=0
      CTOStrng(5)=' fixed random impact parameter'
      CTOption(6)=0
      CTOStrng(6)=' no first collisions inside proj/target'
      CTOption(7)=0
      CTOStrng(7)=' elastic cross-section enabled (<>0:total=inelast)'
      CTOption(8)=0
      CTOStrng(8)=' extrapolate branching ratios '
      CTOption(9)=0
      CTOStrng(9)=' use tabulated pp cross-sections '
      CTOption(10)=0
      CTOStrng(10)=' enable Pauli Blocker'
      CTOption(11)=0
      CTOStrng(11)=' mass reduction for cascade initialization'
      CTOption(12)=0
      CTOStrng(12)=' string condition =0 (.ne.0 no strings)'
      CTOption(13)=0
      CTOStrng(13)=' enhanced file16 output '
      CTOption(14)=0
      CTOStrng(14)=' cos(the) is distributet between -1..1 '
      CTOption(15)=0
      CTOStrng(15)=' allow mm&mb-scattering'
      CTOption(16)=0
      CTOStrng(16)=' propagate without collisions'
      CTOption(17)=0
      CTOStrng(17)=' colload after every timestep '
      CTOption(18)=0
      CTOStrng(18)=' final decay of unstable particles'
      CTOption(19)=0
      CTOStrng(19)=' allow bbar annihilaion'
      CTOption(20)=0
      CTOStrng(20)=' dont generate e+e- instead of bbar'
      CTOption(21)=0
      CTOStrng(21)=' use field feynman frgm. function'
      CTOption(22)=1
      CTOStrng(22)=' use lund excitation function'
      CTOption(23)=0
      CTOStrng(23)=' lorentz contraction of projectile & targed'
      CTOption(24)=2      ! 1 is default    2 means fast method
      CTOStrng(24)=' Wood-Saxon initialization'
      CTOption(25)=0
      CTOStrng(25)=' phase space corrections for resonance mass'
      CTOption(26)=0
      CTOStrng(26)=' use z -> 1-z for diquark-pairs'
      CTOption(27)=1             ! hjd1
      CTOStrng(27)=' reference frame (1=target, 2=projectile, else=cms)'
      CTOption(28)=1             ! M.R. 2019-04-15
      CTOStrng(28)=' propagate spectators also '
      CTOption(29)=2
      CTOStrng(29)=' no transverse momentum in clustr '
      CTOption(30)=1
      CTOStrng(30)=' frozen fermi motion '
      CTOption(31)=0
      CTOStrng(31)='  reduced mass spectrum in string'
      CTOption(32)=0
      CTOStrng(32)=' masses are distributed acc. to m-dep. widths'
      CTOption(33)=0
      CTOStrng(33)=' use tables & m-dep. for pmean in fprwdt & fwidth'
      CTOption(34)=1
      CTOStrng(34)=' lifetme according to m-dep. width'
      CTOption(35)=1
      CTOStrng(35)=' generate high precision tables'
      CTOption(36)=0
      CTOStrng(36)=' normalize Breit-Wigners with m.dep. widths '
      CTOption(37)=0
      CTOStrng(37)=' heavy quarks form di-quark clusters'
      CTOption(38)=0
      CTOStrng(38)=' scale p-pbar to b-bbar with equal p_lab '
      CTOption(39)=0
      CTOStrng(39)=' dont call pauliblocker'
      CTOption(40)=0
      CTOStrng(40)=' read old fort.14 file '
      CTOption(41)=0
      CTOStrng(41)=' generate extended output for cto40'
      CTOption(42)=0
      CTOStrng(42)=' hadrons now have color fluctuations'
      CTOption(43)=0
      CTOStrng(43)=' dont generate dimuon intead of dielectron output'
      CTOption(44)=0
      CTOStrng(44)=' not used at the moment'
      CTOption(45)=0
      CTOStrng(45)=' not used at the moment'

C  INITIALIZE ARRAYS FOR SPECIAL PRO/TAR COMBINATIONS
      do  i = 1, 2
         spityp(i) = 0
         spiso3(i) = 0
      enddo

C  INITIALIZE ARRAYS FOR SPECIAL PARTICLES
      EoS = 0

C  READ CROSS-SECTION FILES
Cdh   CALL URQREC()

C  INITIALIZES SOME ARRAYS
      call strini      ! initialize mixing angles for meson-multipletts
      call loginit

      IF ( CTOption(33) .EQ. 0  .OR.  CTOption(9) .EQ. 0 ) THEN
        call loadwtab(io)
c~         IF ( IUDEBUG .GT. 0 ) WRITE(IFCK,*) 'URQINI: AFTER LOADWTAB'
      ENDIF

C READ URQMD TOTAL CROSS SECTION TABLE
c
c   ie=1..41   E=10.0**(float(ie)/10-1.0-0.05)  (bin-middle)
c   id=1..9    p,ap,n,an,pi+,pi-,K+,K-,KS
c   ia=1..3    N,O,Ar
c
c~       if(ifurqdat.eq.1)then
c~         OPEN(UNIT=76,FILE=furqdat(1:nfurqdat),STATUS='OLD')
c~       else
c~         OPEN(UNIT=76,FILE='UrQMD-1.3.1-xs.dat',STATUS='OLD')
c~       endif
c~       read(76,*) adum,iamaxu,idmaxu,iemaxu
c~       do ia=1,iamaxu
c~         do id=1,idmaxu
c~           do ie=1,iemaxu
c~             read(76,*) ekdummy,sig_u1(ie,id,ia)
c~           enddo
c~           read(76,*)
c~           read(76,*)
c~         enddo
c~       enddo
c~       close(76)

C  IN CASE OF CASCADE MODE, THE POTENTIALS NEED NOT BE CALCULATED

C  CALCULATE NORMALIZATION OF RESONANCES DISTRIBUTION...
      call norm_init
#endif


c~       xsegymin=0.25d0

#ifdef __CXDEBUG__
      call utisx2
#endif

      end

