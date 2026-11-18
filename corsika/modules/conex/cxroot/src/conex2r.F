#include "conexConfig.h"
c M. Unger, T. Pierog 03.2005
c 31.03.2005 PDG code input format for primary particle (T. Pierog)
c 28.02.2005 implicit none for conex.inc (T. Pierog)
c 14.02.2005 added get_shower_muon (T. Pierog)
c 17.08.2004 added get_shower_edep (T. Pierog)
c 08.09.2004 removed hardcoded ".param" file names (M.U.)
c 18.05.2005 added get_shower_gamma (RU)
c 18.05.2005 added get_shower_electron (RU)
c 02.09.2005 added QGSJETII
c 25.07.2006 added EPOS and link to conex3
c 13.10.2007 added low energy MC parameters and output (KASCADE)
c 01.07.2008 added MuMia (rho600) for Mia detector
c 11.11.2008 added S1000 for Auger
c 14.02.2013 correct interface to new models in 3D
c 07.07.2017 added DPMJETIII
c 07.07.2020 add CORSIKA8 specific command (TP)
c 04.04.2022 stop wit exit code 3 (M. Reininghaus)
c 02.02.2025 add QGSJETIII+++++++++++++++++++++++++++++++++ (TP)
c 18.02.2025 replace "stop 3" by cxstop (TP)

#ifdef LEADING_INTERACTIONS_TREE
#ifndef __SORT_FOR_ENERGY__
#error It makes no sense to use LEADING_INTERACTIONS_TREE without __SORT_FOR_ENERGY__
#endif
#endif

      subroutine InitConex(nShow,seedVal,model,nMaxDetail,
#ifdef CONEX_EXTENSIONS
     &                     partListModePar, ! RU Mon Oct 23 11:11:45 CEST 2006
#endif 
     &                     conex_root)
***************************************************************
*
*  initialization of Conex
*
*   input: nShow      - number of showers to generate
*          iseed      - random seed
*          model      - HE model switch
*          conex_root - path of directory containing data
*                       card files
*
***************************************************************
*
      implicit none
#include "conex.h"
      integer nShow,iRndmSeed(3),seedVal,model
      character*(*) conex_root
*
      integer inLun,outLun
      parameter (inLun=40,outLun=41)
      character*800 fileName
*
#ifdef CONEX_EXTENSIONS
      integer partListModePar         ! RU Mon Oct 23 11:12:16 CEST 2006
#endif
      integer nMaxDetail

#ifdef CONEX_EXTENSIONS
      particleListMode=partListModePar
      writeFirstIntPart=-1
      interactionCounter=0
      isFirstInt=.false.
      if(particleListMode.ge.2)then
        writeFirstIntPart=particleListMode-2
      endif
c     CONEX_EXTENSIONS
#endif 
      
      if ( model .eq. 2 ) then
         fileName=conex_root//'/conex_qgsjet.param'
      elseif ( model .eq. 1 ) then
         fileName=conex_root//'/conex_nexus.param'
      elseif ( model .eq. 4 ) then
         fileName=conex_root//'/conex_epos.param'
      elseif ( model .eq. 5 ) then
         fileName=conex_root//'/conex_sibyll.param'
      elseif ( model .eq. 6 ) then
         fileName=conex_root//'/conex_qgsjetII.param'
      elseif ( model .eq. 9 ) then
         fileName=conex_root//'/conex_dpmjet.param'
      elseif ( model .eq. 10 ) then
         fileName=conex_root//'/conex_qgsjetIII.param'
      else
         write(6,*)' conex_f: Unknown model : ',model
         call cxstop('          exit ...')
      endif

      iRndmSeed(1) = seedVal
      iRndmSeed(2) = 0
      iRndmSeed(3) = 0

      open(unit=inLun,file=fileName(1:INDEX(filename,' ')-1)
     .    ,status='old',err=999)
      open(unit=outLun,file='/dev/null',status='unknown',err=999)
c RU: not sure if this is still needed, but it doesn't work like this anymore:
c      call rmmaq(iStartSeed,"S")
      call ConexInit(nShow,nMaxDetail,inLun,outLun,iRndmSeed)

      return

 999  write(6,*)' InitConex: Error opening parameter file ',fileName
      call cxstop('          exit ...')
*                                                      end InitConex
      end

      subroutine GetSeed(iRndmSeed)
***************************************************************
*
*  interface to conex subroutine
*
*   output: iRndmSeed      - random seed
***************************************************************
*
      implicit none
      common/cxransto/diu0(100),iiseed(3,2)
      integer iRndmSeed(3),iiseed,i
      double precision seed,diu0

      call ranfgt(seed)         !get seed before shower
      do i=1,3
        iRndmSeed(i)=iiseed(i,1)
      enddo

      return
*                                                                end RunConex
      end

      subroutine RunConex(iRndmSeed,e,theta,phi,impact,ipart)
***************************************************************
*
*  interface to conex subroutine
*
*   input: iseed      - random seed
*          e          - particle energy in GeV
*          theta      - zenith angle (deg)
*          phi        - azimuth angle (deg) (AUGER def : 0 = East)
*          impact     - impact parameter
*          ipart      - primary particle type
*
***************************************************************
*
      implicit none
      integer ipart,iRndmSeed(3),idum(3)
      double precision e, theta, phi, impact
*
#ifdef CONEX_EXTENSIONS
#include "conex.h"
      interactionCounter=0
      isFirstInt=.false.
#endif
      call ConexRun(ipart, e, theta, phi, 0d0, impact, idum)

      return
*                                                                end RunConex
      end


      subroutine get_shower_data(icut,ipart,ndepth,depth,profile,fitres
     >                           ,alt,dis)
**************************************************************************
*
*  interface to retrieve data from conex
*  common blocks
*
*   input: cut number "icut"
*          secondary particle type "ipart"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: number of slant depth points "ndepth"
*           slant depth "depth"
*           particle profile "profile"
*           altitude above see level (m) "alt"
*           slant distance to obs. point (or closer point to Earth) "dis"
*             dis (m) > 0 : shower goes to the obs. point
*             dis (m) < 0 : shower comes from the obs. point
*
***************************************************************************
      implicit none
#include "conex.h"
*
      integer ndepth,icut,ipart,ntmp,i,j
      real depth(*),profile(*),fitres(*),alt(*),dis(*)
      double precision heightt,distL,distance0
      external heightt,distance0

      ntmp=nmaxX-nminX+1
      if(ntmp.gt.ndepth) then
         print *," Error in getShowerData - too many slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop('          exit ...')
       endif
       ndepth=ntmp
       do i=nminX,nmaxX 
         j=i-nminX+1
         depth(j)=sngl(zha(i))
         profile(j)=sngl(XProf(i,icut,ipart)) 
         distL=distz(i)                 
         alt(j)=sngl(heightt(abs(distL),radtr0))
         dis(j)=sngl(distL)
       enddo
       fitres(1)=sngl(log10(eprima))+9.
       fitres(2)=sngl(thetas)
       fitres(3)=sngl(xfirst)           !X1 (first interaction point)
       fitres(4)=sngl(XmaxShow(4,1))    !Nmax
       fitres(5)=sngl(XmaxShow(4,2))    !X0
       fitres(6)=sngl(XmaxShow(4,4))    !P1
       fitres(7)=sngl(XmaxShow(4,5))    !P2
       fitres(8)=sngl(XmaxShow(4,6))    !P3
       fitres(9)=sngl(XmaxShow(3,0))    !chi2/sqrt(Nmax)
       fitres(10)=sngl(XmaxShow(4,3))   !Xmax
       fitres(11)=sngl(phisho)          !azimuth angle
       fitres(12)=sngl(xfirstin)        !inelasticity of 1st interaction
       distL=distance0(xfirst)
       fitres(13)=sngl(heightt(abs(distL),radtr0))
*                                                    end get_shower_data
      end


      subroutine get_data_card1(icut,hmCut,egCut,icutm,mCut,icutg,gCut
     & ,icute,eCut,icuth,hcut,hThr,mThr,eThr,HiLowE,leModel,cVersion)
****************************************************************
*
*  interface to retrieve data card options (1)
*  from conex  common blocks
*
*   input: cut number "icut (min in simu), icutm, icutg, icute, icuth"
*
*   output: hadron/muon and emag cutoffs hmCut and egCut for simu,
*           muon, hadron, electron and gamma cutoffs mCut, hCut, eCut and gCut for profiles,
*           hadron, muon and emag threshold (MC->CE) hThr, mThr, eThr
*           low energy model number leModel and conex and bin size dX
*           version cVersion
*
****************************************************************
      implicit none
#include "conex.h"   
      real mCut,egCut,hmCut,hCut,eCut,gCut,hThr,mThr,eThr
     *         ,cVersion,HiLowE
      integer icut, icutm, icute, icutg, icuth, leModel

      hmCut=sngl(HaCutP(icut))
      egCut=sngl(EMCutP(icut))
      gCut=sngl(EMCutP(icutg))
      eCut=sngl(EMCutP(icute))
      mCut=sngl(HaCutP(icutm))
      hCut=sngl(HaCutP(icuth))
      hThr=sngl(fehcut)
      mThr=sngl(femcut)
      eThr=sngl(feecut)
      HiLowE=sngl(EgyHiLoLim)
      leModel=MCleModel
      cVersion=float(ivers)/100.

*                                                end get_data_card1
      end 

#if __MC3D__

      subroutine get_data_card2(WHmx,WMmx,WTmx,hLow,mLow,egLow,zsLow
     &                         ,iThn,Thn,hgro)
****************************************************************
*
*  interface to retrieve data card options (2)
*  from conex  common blocks
*
*   output: hadron/muon and emag maximum weight in MC,
*           hadron, muon and emag minimum energy for CE (CE->MC) hLow, mLow, 
*           egLow and zslow
*           thinning flag and threshold for EM MC
*           ground height hgrd and height of shower axis, 
*           at (x=0,y=0) above ground
*
****************************************************************
      implicit none
#include "conex.h"   
      real WHmx,WMmx,WTmx,hLow,mLow,egLow,zsLow,Thn,hgro,alt
      integer iThn
      WHmx=sngl(fwhmax)
      WMmx=sngl(fwmmax)
      WTmx=sngl(fwemax)
      hLow=sngl(ehlowi)
      mLow=sngl(emlowi)
      egLow=sngl(eelowi)
      zsLow=sngl(zshlow)
      iThn=iothin
      Thn=sngl(thin)
      hgro=sngl(hground)

*                                                end get_data_card2
      end 

#ifndef __CORSIKA8__
      
      subroutine get_shower_lat(icuth,sEh,Ehmax,icutm,MuTr,MuMi,S1000)
**********************************************************
*
*  interface to retrieve muon profiles 
*  from conex  common blocks
*
*   input: cut number "icuth" for hadrons and  "icutm" for muons
*
*   output: sum of hadron energy
*           at ground and Maximum energy of the hadron at ground
*           Truncated number of muons at ground and S1000
*
**********************************************************
      implicit none
#include "conex.h"   
*
      integer icuth,icutm,i
      real    sEh,Ehmax,Mutr,MuMi,S1000(0:4)  
*    

       sEh=sngl(EHaSum(icuth))
       Ehmax=sngl(EHaMax)
       MuTr=sngl(MuTrunc(icutm))
       MuMi=sngl(MuMia)
       do i=0,4
         S1000(i)=sngl(SD1000(i))
       enddo
*                                                   end get_shower_lat
       end

#endif

#endif

      subroutine get_shower_muon(icut,ndepth,muProf,muProd)
**********************************************************
*
*  interface to retrieve muon profiles
*  from conex  common blocks
*
*   input: cut number "icut"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: Number of muons for cutoff icut and production rate
*
**********************************************************
      implicit none
#include "conex.h"   
*
      integer ndepth,icut,ntmp,i,j
      real muProf(*),muProd(*)  
*    

      ntmp=nmaxX-nminX+1
      if(ntmp.ne.ndepth) then
         print *," Error in get_shower_muon - inconsistent",
     >       " slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop(' exit ...')
       endif

       do i=nminX,nmaxX 
         j=i-nminX+1
         muProf(j)=sngl(XProf(i,icut,3)) 
         muProd(j)=sngl(XdMu(i)*delzsh) 
       enddo
*                                                            end get_shower_muon
       end

      subroutine get_shower_edep(icut,ndepth,edepcx,eground)
**********************************************************
*
*  interface to retrieve energy deposit profile
*  from conex  common blocks
*
*   input: cut number "icut"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: Total energy deposit profile "Edep" in GeV/g.cm2
*
**********************************************************
      implicit none
#include "conex.h"
*
      integer ndepth,icut,ntmp,i,j
      real EdepCx(*),eground(*)
*

      ntmp=nmaxX-nminX+1
      if(ntmp.ne.ndepth) then
         print *," Error in get_shower_edep - inconsistent",
     >       " slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop(' exit ...')
       endif

       do i=nminX,nmaxX 
         j=i-nminX+1
         EdepCx(j)=sngl(max(0.d0,Edepo(i,0))) 
       enddo
*
** retrieve rest energy at end of shower
*
       eground(1)=sngl(Ebalan(nmaxX,1))  ! emag
       eground(2)=sngl(Ebalan(nmaxX,2))  ! hadrons
       eground(3)=sngl(Ebalan(nmaxX,4))  ! muons
*                                                            end get_shower_edep
       end



      subroutine get_shower_gamma(icut,ndepth,gaProf)
**********************************************************
*
*  interface to retrieve gamma profile 
*  from conex  common blocks
*
*   input: cut number "icut"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: Number of gammas for cutoff icut
*
**********************************************************
      implicit none
#include "conex.h"   
*
      integer ndepth,icut,ntmp,i,j
      real gaProf(*)    
*    

      ntmp=nmaxX-nminX+1
      if(ntmp.ne.ndepth) then
         print *," Error in get_shower_gamma - inconsistent",
     >       " slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop(' exit ...')
       endif

       do i=nminX,nmaxX 
         j=i-nminX+1
         gaProf(j)=sngl(XProf(i,icut,1)) 
       enddo
*                                                            end get_shower_gamma
       end

      subroutine get_shower_electron(icut,ndepth,elProf)
**********************************************************
*
#ifdef __CORSIKA8__
*     interface to retrieve electron profile
#else
*     interface to retrieve electron+positron profile
#endif
*  from conex  common blocks
*
*   input: cut number "icut"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: Number of electrons for cutoff icut
*
**********************************************************
      implicit none
#include "conex.h"   
*
      integer ndepth,icut,ntmp,i,j
      real elProf(*)    
*    

      ntmp=nmaxX-nminX+1
      if(ntmp.ne.ndepth) then
         print *," Error in get_shower_electron - inconsistent",
     >       " slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop(' exit ...')
       endif

       do i=nminX,nmaxX 
         j=i-nminX+1
#ifdef __CORSIKA8__
         elProf(j)=sngl(XProf(i,icut,-1))
#else
         elProf(j)=sngl(XProf(i,icut,2))
#endif
       enddo
*                                                            end get_shower_electron
       end

#ifdef __CORSIKA8__

      subroutine get_shower_positron(icut,ndepth,elProf)
**********************************************************
*
*  interface to retrieve positron profile 
*  from conex  common blocks
*
*   input: cut number "icut"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: Number of electrons for cutoff icut
*
**********************************************************
      implicit none
#include "conex.h"   
*
      integer ndepth,icut,ntmp,i,j
      real elProf(*)    
*    

      ntmp=nmaxX-nminX+1
      if(ntmp.ne.ndepth) then
         print *," Error in get_shower_positron - inconsistent",
     >       " slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop(' exit ...')
       endif

       do i=nminX,nmaxX 
         j=i-nminX+1
         elProf(j)=sngl(XProf(i,icut,2)) 
       enddo
*                                                            end get_shower_electron
       end

#endif

      subroutine get_shower_hadron(icut,ndepth,haProf)
**********************************************************
*
*  interface to retrieve hadron profile 
*  from conex  common blocks
*
*   input: cut number "icut"
*          dimension "ndepth" of arrays "depth"/"profile"
*
*   output: Number of hadrons for cutoff icut
*
**********************************************************
      implicit none
#include "conex.h"   
*
      integer ndepth,icut,ntmp,i,j
      real haProf(*)    
*    

      ntmp=nmaxX-nminX+1
      if(ntmp.ne.ndepth) then
         print *," Error in get_shower_hadron - inconsistent",
     >       " slant depth bins!"
         print *," nmaxX,nminX,ndepth",nmaxX,nminX,ndepth
         call cxstop(' exit ...')
       endif

       do i=nminX,nmaxX 
         j=i-nminX+1
         haProf(j)=sngl(XProf(i,icut,4)) 
       enddo
*                                                            end get_shower_hadron
       end

      integer function get_number_of_depth_bins()
**************************************************************************
*
*  returns current number of depth bins from CONEX common blocks
*
***************************************************************************
      implicit none
#include "conex.h"
*
      get_number_of_depth_bins = nmaxX-nminX+1
      return
*                                               end get_number_of_depth_bins
      end
