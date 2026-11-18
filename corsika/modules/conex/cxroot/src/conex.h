c-----------------------------------------------------------------------------
c                       Main include file for conex
c
c Last modification 16.05.2006 by T. Pierog (ver 3.06)
c Author R. Engel, N. Kalmykov, S. Ostaptchenko, T. Pierog, K. Werner
c-----------------------------------------------------------------------------


      integer mxblk,mxstk,mxstkj,maximE,maximZ
      parameter (mxblk=16,mxstkj=(mxblk+1)*50*2)
#ifdef __SAVEMEMO__
      parameter (maximE=281,maximZ=201)
      parameter (mxstk=mxblk*50*2)
#else
      parameter (maximE=281,maximZ=3701)
      parameter (mxstk=mxblk*100*2)
#endif
      integer maximEd

      parameter (maximEd=261)     !do not change without updating dky and dedx tables

      integer iseed,lseq
      common /cxseed/iseed(3,2),lseq

      common /cxrrandpa/facrdm,u1rdm,u2rdm,knordm
      double precision facrdm,u1rdm,u2rdm
      logical          knordm

c shower basics

      double precision dptl,stack,stackj,eecut,epcut,ehcut,emcut
     &,feecut,fehcut,femcut
      double precision eelow,ehlow,emlow,zshlow,zmclow,zmchlow
     &,eelowi,ehlowi,emlowi,fwhmax,fwemax,fwmmax
     &,wshmax,wsmmax,wsemax
      integer nshower,mshow,mshowEGS,istack,irec,ifsa
     &,jstack,jrec,jfsa
      logical cx2corsha,cx2corsmu,cx2corsem

      common/cxsho1/nshower,mshow,mshowEGS
      common /cxstack/stack(mxstk),istack,irec,ifsa
      common /cxstackj/stackj(mxstkj),jstack,jrec,jfsa
      common /cxoptl/dptl(mxblk)
        ! 1-5 ... 5-momentum
        ! 6,7 ... lateral coordinates x,y (m) in obs frame
        ! 8 ..... altitude
        ! 9 ..... time (m)
        ! 10 .... id
        ! 11 .... weight
        ! 12 .... generation
        ! 13 .... slant depth passed (g/cm^2)
        ! 14,15 . lateral coordinates x,y (m) in shower frame
      common /cxcut/ eecut,epcut,ehcut,emcut  !cutoff el,phot,had,mu also in egs4
      common /cxsubcut/ feecut,fehcut,femcut !factor for cutoff em and had for subroutine mode
c cutoff el,had,mu,z and weight for low energy MC
      common /cxlow/ eelow,ehlow,emlow,eelowi,ehlowi,emlowi
     &    ,zshlow,zmclow,zmchlow,fwhmax,fwemax,fwmmax
     &    ,wshmax,wsemax,wsmmax
     &           ,cx2corsha,cx2corsmu,cx2corsem

c basics

      double precision estck,zshmin,zshmax,delzsh,hground,enymin
     &,enymax,decade,emin,exmin
      double precision eprima,thetas,costhet,phisho,c2bas
     &,sinthet,cosphi,sinphi,dphmin0,XminSlant,HGrd,distMaxi
      double precision altitude,RadGrd,DistALt,dphmaxi0,dphlim0
     &,latitude,longitude
      real year
      integer irdelz,muse,musz,ipreshow
#ifdef __MC3D__
     &,musLh,musLm
#endif
      logical goOutGrd

      common /cxstene/estck(2)
      common/cxbas3/zshmin,zshmax,delzsh,hground,enymin,enymax,decade
     &,emin,irdelz
      common/cxbas4/eprima,thetas,costhet,phisho,muse,musz,c2bas
     &,sinthet,cosphi,sinphi,XminSlant,HGrd,distMaxi
#ifdef __MC3D__
     &,musLh,musLm
#endif
      common/cxbas5/exmin
      common/cxbas6/altitude,RadGrd,DistALt,dphmaxi0,dphlim0,dphmin0
     &,goOutGrd
      common/cxbas7/latitude,longitude,year,ipreshow

c source

      double precision hsource,EgyHiLoLim,p2ha
      integer isoumax,jsoumin,MCModel,MCleModel,ilowegy
      character*500 nexusdir

      common/cxhsource/hsource(8,maxime,maximz),p2ha(8,maxime)
     &                ,isoumax,jsoumin
        ! 1 ... protons
        ! 2 ... charged pions
        ! 3 ... charged kaons
        ! 4 ... neutral kaon Klong
        ! 5 ... neutral kaon Kshort
        ! 6 ... neutrons
        ! 7 ... muons
        ! 8 ... gamma from photonuclear effect
      common/cnexdir/nexusdir
      common/cxmodel/EgyHiLoLim,MCModel,MCleModel,ilowegy
        ! 1 ... NeXuS
        ! 2 ... QGSJet
        ! 3 ... Gheisha
        ! 4 ... Sibyll

#if __MC3D__ || __CXLATCE__
      double precision hpt2source
      common/cxhpt2source/hpt2source(8,maxime,maximz)
#endif

#ifdef __ANALYSIS__
c plots

      double precision hadspec
      common/cxhadcas/hadspec(7,maxime,maximz)
#if __MC3D__ || __CXLATCE__
      double precision ptspec
      common/cxhadcaspt/ptspec(2,7,maxime,maximz)

#endif
#endif

c utilities

      double precision airz,aira,airw,airavz,airava,AATM,BATM,CATM
     &,DATM,EATM,pi,radlth,radearth,cxlight,xgauss7,wgauss7,airi,falpha
     &,fialpha,xgauss10,wgauss10
      double precision bdeca,xsaxis,ysaxis,zsaxis,avog,bbatm,ccatm
      integer mxatm,mnatm,mode,iwrt,i1DMC,iphonu


      common/cxair/airz(3),aira(3),airw(3),airavz,airava,airi(3)!also in geisha_conex and nexus_conex and egs4_conex
      parameter(mxatm=5)                   !so241103
      common/cxatmos/AATM(mxatm),BATM(mxatm),CATM(mxatm),DATM(mxatm)
     &,EATM(mxatm+1)    !so120903
      common/cxatmos2/bbatm(mxatm),ccatm(mxatm),mnatm
      common/cxdeca/ bdeca(18)                                   !tp110414
      common/cxaxis/ xsaxis,ysaxis,zsaxis
      common/cxetc/mode,iwrt,i1DMC,iphonu                      !also in egs4
      common/cxconst/pi,radlth,radearth,avog,falpha
      common/cxconst2/cxlight,fialpha
      common/cxgauss/xgauss7(7),wgauss7(7)
      common/cxdga20/xgauss10(10),wgauss10(10)

c files

      integer nfnho,nfnck,nfnwle,nfnwhe,nfndkz,nfndkl,nfndks
     &,nfndkm,nfnilo,nfndke,nfndkn,nfndkg,nfnwgh,nfnwgl
     &,nfnda,nfnemcs,nfninput,nfnrt
      integer ifho,ifck,ifwle,ifwhe
     &,ifdkz,ifdkl,ifdks,ifdkm,ifilo
     &,ifdke,ifdkn,ifdkg,ifwgh
     &,ifwgl,ifda,ifemcs,ifinput,ifrt
      character*500 fnho,fnck,fnwle,fnwhe,fndkz,fndkl,fndks,fndkm
     &,fnilo,fndke,fndkn,fndkg,fnwgh,fnwgl      !also in nexus/epos_conex and gheisha_conex
      character*500 fnda,fnemcs,fninput,fnrt

      common /cxfiles/fnho,ifho,fnck,ifck,fnwle,ifwle,fnwhe,ifwhe
     &,fndkz,ifdkz,fndkl,ifdkl,fndks,ifdks,fndkm,ifdkm,fnilo,ifilo
     &,fndke,ifdke,fndkn,ifdkn,fndkg,ifdkg,fnwgh,ifwgh
     &,fnwgl,ifwgl
      common /cxnfiles/nfnho,nfnck,nfnwle,nfnwhe,nfndkz,nfndkl,nfndks
     &,nfndkm,nfnilo,nfndke,nfndkn,nfndkg,nfnwgh,nfnwgl
      common /cxfiles2/fnda,ifda,fnemcs,ifemcs
      common /cxnfiles2/nfnda,nfnemcs
      common /cxfiles3/fninput,ifinput,fnrt,ifrt
      common /cxnfiles3/nfninput,nfnrt

      integer nfnp2le,nfnp2he,nfnp2d,ifp2le,ifp2he,ifp2d
c     &,nfnp4le,nfnp4he,nfnp4d,ifp4le,ifp4he,ifp4d
      character*500 fnp2le,fnp2he,fnp2d
c     &,fnp4le,fnp4he,fnp4d

      common /cxfiles4/fnp2le,ifp2le,fnp2he,ifp2he,fnp2d,ifp2d
c     &,fnp4le,ifp4le,fnp4he,ifp4he,fnp4d,ifp4d
      common /cxnfiles4/nfnp2le,nfnp2he,nfnp2d
c     &,nfnp4le,nfnp4he,nfnp4d


c debug

      integer mxisx,isx,nisx,isxsave,isxxsave,isxsub
      character*500 subisx,textisx

      parameter (mxisx=200)
      common/cxisx/isx,nisx,subisx(mxisx),isxsub(mxisx)
     &      ,isxsave,isxxsave,textisx                    !also in gheisha_nexus

c effects

      double precision thin,ethin,wtmax,rthmax,hthin,ehthin,whmax
     &,radtr0,depthmaxi0,pmass,sterncor,etotlost,etotsource,etotsta
      integer ilpmeffect,iothin,ihthin,istern,i1DEM,imscat,ionloss
     &,isubin,ifragm,iMuInt,iLatCE
#ifdef __CXCORSIKA__
      logical lFlat
#endif

      common /cxlpm/ ilpmeffect
      common /cxthin/ thin,ethin,wtmax,rthmax,iothin   !egs4
      common /cxHadthin/ hthin,ehthin,whmax,ihthin
      common /cxcoord/ radtr0,depthmaxi0
      common /cxmass/ pmass(10) !mass (proton, pion ch, kaon ch, kaon 0, pion 0, neutron, nucleon,lambda,muon,electron)
      common/cxstern/sterncor,istern        !also in egs4 and conex-eph
      common/cx1dem/i1DEM,imscat,ionloss    !also in egs4
      common/cxegybal/etotlost,etotsource,etotsta
      common/cxsubro/isubin                 !also in nexus_conex.fpp
      common/cxfragm/ifragm
      common/cxmuint/iMuInt
      common/cxlatce/iLatCE
#ifdef __CXCORSIKA__
      common/cxflatg/lFlat
#endif

#ifdef __ANALYSIS__
#if __MC3D__ || __CXLATCE__
      integer maxo,maxmm,muso,musmm,iimom,i1mom,i2mom,i3mom
#ifdef __CXLATCE__
      parameter(maxo=3,maxmm=(maxo+1)*(maxo**2+5*maxo+6)/6-1)
#else
      parameter(maxo=2,maxmm=maxo)
#endif
      common/cxmom1/muso,musmm,iimom(0:maxo,0:maxo,0:maxo)
     &           ,i1mom(0:maxmm),i2mom(0:maxmm),i3mom(0:maxmm)
#else
      integer maxo,maxmm
      parameter(maxo=0,maxmm=0)
#endif
#endif

c  shower analysis


#ifdef __ANALYSIS__
      integer maxiz,numiz
      integer ngenmx
      parameter (ngenmx=100)
      double precision cntgen(0:ngenmx)
      common/countgen/cntgen
      double precision zamin,zamax,delza,yieldz,yiex,spec
     &,tamin,tamax,ctime,eamin,eamax
      integer maxjz,maxin,maxie,izfirst,modz,numie
     &,kfirsth,modkh

      common/cxtime/ tamin,tamax,ctime
      common/cxlimitsh/kfirsth,modkh
#if __MC3D__ || __CXLATCE__
      double precision ramin,ramax,xamin,xamax,yieldr,yieldx
     &,specr,spex,speca,yieldr2,yieldr1
#ifdef __CXLATCE__
     &,yieldrt,yieldrt1,yieldrt2
#endif
      integer maxir,maxjr,maxiex,maxix,numix,numir,irfirst,modr
     &,iefirst,moden

      parameter(maxiz=401,maxjz=3,maxir=101,maxjr=3
     &          ,maxie=151,maxix=101,maxin=12,maxiex=5 )
      common/cxlimits/zamin,zamax,delza,ramin,ramax
     & ,eamin(2),eamax(2),xamin(maxiex),xamax(maxiex),numix(maxiex)
     & ,numie,numir,irfirst,modr,numiz,izfirst,modz,iefirst,moden
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz)
     & ,yieldr(maxin,maxjz,maxir),yieldr1(maxin,maxjz,maxir)
     & ,yieldr2(maxin,maxjz,maxir)
     & ,yieldx(maxiex,maxin,maxjz,maxjr,maxix)
     & ,spec(0:maxmm,maxin,maxie,maxjz),spex(0:maxmm,maxin,maxie,maxjz)
     & ,specr(maxin,maxie,maxjz,maxjr),speca(maxin,maxjz,maxjr)
#ifdef __CXLATCE__
     & ,yieldrt(maxin,maxjz,maxir),yieldrt2(maxin,maxjz,maxir)
     & ,yieldrt1(maxin,maxjz,maxir)
#endif
#else
      parameter(maxiz=401,maxie=151,maxin=12,maxjz=3 )
      common/cxlimits/zamin,zamax,delza
     & ,eamin(2),eamax(2),numie,numiz,izfirst,modz
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz)
     & ,spec(0:maxmm,maxin,maxie,maxjz)
#endif
#endif

c shower output

      double precision XProf,XmeanP,XmeanP2,EMCutP,HaCutP,XminP,XmaxP
     &,Ebalan,Ebalan1,Ebalan2,Einit,XmaxShow,XmaxMean,XmaxProf,Xfirst
     &,Edepo,Edepo1,Edepo2,XdMu,XdMuMean,XdMuMean2,XfirstIn
      integer mxExpro,mxPxpro,nminX,nmaxX,mZEMHa,ifout,ivers,iXmax
      logical lheader,lxfirst,lxfirstIn,lwrite

#if __CXCORSIKA__ || __CORSIKA8__
      parameter(mxExpro=1,mxPxpro=10)
      common/cxoutput3/XmaxShow(4,-1:mxPxpro),XmaxMean(5,-1:mxPxpro)
     &,XmaxProf(maximz,-1:mxPxpro),iXmax
      common/cxoutput1/XProf(maximz,mxExpro,-1:mxPxpro)
     &,XmeanP(maximz,mxExpro,-1:mxPxpro)
     &,XmeanP2(maximz,mxExpro,-1:mxPxpro)
     &,EMCutP(mxExpro),HaCutP(mxExpro),XminP,XmaxP
     &,nminX,nmaxX,mZEMHa,ifout,ivers,lheader
#else
      double precision MuTrunc,EHaSum,EHaMax,MuMia,SD1000,SD1000m
      parameter(mxExpro=3,mxPxpro=10)
      common/cxoutput3/XmaxShow(4,0:mxPxpro),XmaxMean(5,0:mxPxpro)
     &,XmaxProf(maximz,0:mxPxpro),iXmax
      common/cxoutput1/XProf(maximz,mxExpro,0:mxPxpro)
     &,XmeanP(maximz,mxExpro,0:mxPxpro)
     &,XmeanP2(maximz,mxExpro,0:mxPxpro)
     &,EMCutP(mxExpro),HaCutP(mxExpro),XminP,XmaxP
     &,nminX,nmaxX,mZEMHa,ifout,ivers,lheader
#endif
      common/cxoutput2/Ebalan(maximz,4),Ebalan1(maximz,4)
     &,Ebalan2(maximz,4),Einit
      common/cxoutput4/Xfirst,XfirstIn,lxfirst,lxfirstIn          !also in egs4
      common/cxoutput5/Edepo(maximz,0:mxPxpro)
     &,Edepo1(maximz,0:mxPxpro)
     &,Edepo2(maximz,0:mxPxpro)
      common/cxoutput6/XdMu(maximz),XdMuMean(maximz),XdMuMean2(maximz)
#if !__CXCORSIKA__ && !__CORSIKA8__
     &,MuTrunc(mxExpro),EHaSum(mxExpro),EHaMax,MuMia,SD1000(0:4)
     &,SD1000m(0:4)
#endif
      common/cxoutput0/lwrite

c Hadronic cascade equations

      double precision  dzHa,Cha,c2ha,zha,eeha,rnHa,ppHa,rpHa,rkz,rkl
     *,rks,HaMu,pi0,gamHa,dnHa,distz,rhoz,distmid,rhora,delLe
     *,z2mid,rlamti,wwHa,wwdec,allldec,antibars
      integer n1maxi,n2maxi,nmaxHa,iemin,iemax,iehlim,nminHa
     *,minEhad,mzHa,iehmc,iemmc,mzmc,ienmn

      parameter (n1maxi=7,n2maxi=12)
      common /cxarea2/ dzHa,Cha,c2ha,zha(maximz),eeha(maxime)
      common /cxarea3/ rnHa(maxime,2),ppHa(maxime,2)
     *,rpHa(maxime,2),rkz(maxime,2),rkl(maxime,2)
     *,rks(maxime,2),HaMu(maxime,2),pi0(maxime,2),gamHa(maxime)
      common /cxarea7/ dnHa,nmaxHa,iemin,iemax,iehlim,nminHa,ienmn
      common /cxarea8/ distz(maximz),rhoz(maximz),distmid,rhora,delLe
     *,z2mid
      common /cxinecs/ rlamti(n1maxi,maxime)
     *,wwHa(maxime,maxime,n1maxi,n2maxi)
      common /cxdecsp/wwdec(maxime,maxime,7),allldec(maxime,n2maxi)
      common /cxlinkemha/minEhad
      common/cxantibarsou/antibars(maximZ)
      common /cxarea5/ mzHa
      common /cxarea6/ iehmc,iemmc,mzmc


c decay tables

      double precision akz0,akl0,aks0,akzm,aklm,apim,akze,akle,ap0g
     *,ap0e,amue,akz,akl,aks,akzn,akln,apin,amun,agpi,agpr,agne,agmu
      integer ndecmax,ndecmin

      common /cxarea0dec/ akz0(maximEd,maximEd),akl0(maximEd,maximEd)
     *,aks0(maximEd,maximEd)
      common /cxarea1dec/ akzm(maximEd,maximEd),aklm(maximEd,maximEd)
     *,apim(maximEd,maximEd)
      common /cxarea2dec/ akze(maximEd,maximEd),akle(maximEd,maximEd)
     *,ap0g(maximEd,maximEd),ap0e(maximEd,maximEd)
     *,amue(maximEd,maximEd)
      common /cxarea3dec/ akz(maximEd,maximEd),akl(maximEd,maximEd)
     *,aks(maximEd,maximEd)
      common /cxarea4dec/ akzn(maximEd,maximEd),akln(maximEd,maximEd)
     *,apin(maximEd,maximEd),amun(maximEd,maximEd)
      common /cxarea5dec/ agpi(maximEd,maximEd),agpr(maximEd,maximEd)
     *,agne(maximEd,maximEd),agmu(maximEd,maximEd)
      common/cxarea6dec/ndecmin(maximEd,n1maxi)
     *                 ,ndecmax(maximEd,n2maxi)

C  Ionization loss (part,energy)

      double precision dedxion,dedxionmu
      common/cxdedx/dedxion(4,maximEd),dedxionmu(maximEd,maximZ)
        ! 1 ... proton
        ! 2 ... charged pions
        ! 3 ... charged kaons
        ! 4 ... muons (for interaction part of Eloss)

#if __MC3D__ || __CXLATCE__
      double precision pt2mu
      common/cxarpt2/pt2mu(maxime,maxime,4)

      double precision pt2w,pt2pi
      common/cxarpt3/pt2w(maxime,maxime,n1maxi,n2maxi)     !so050207
     *,pt2pi(3,maxime,maxime)
#endif
#if __MC3D__ || __CXLATCE__
      double precision the2ha
      common/cxarprod/the2ha(8,maxime)
#endif


c Magnetic Field
      double precision Bfield(4,maximz)
      integer iMagne
      common /cxmagnef/Bfield,iMagne

c --------------/cxmumult/-----------------------
c  chc         = constant chi_c   for muomn multiple scattering
c  omc         = constant omega_c for muomn multiple scattering
c  phisct      = azimutal angle of muon multiple scattering
c  cphisct     = cosine of azimutal angle of muon multiple scattering
c  sphisct     = sine   of azimutal angle of muon multiple scattering
c  stepl       = step length for muon transport step
c  vscat       = polar angle of muon multiple scattering
c  iMuScat     = flag for muon multiple scattering
      double precision chc,omc,phisct,stepl,vmscat,eneper1
      integer iMuScat
      common /cxmumult/chc,omc,phisct,stepl,vmscat,eneper1,iMuScat

#ifdef LEADING_INTERACTIONS_TREE
      integer maxDetail, countInt
      logical leadingParticle
      common /firstinttree/maxDetail,countInt,leadingParticle
#endif


#ifdef CONEX_EXTENSIONS
c     RU Mon Oct 23 09:05:39 CEST 2006
!      parameter (nTmp=16) 
!      double precision ruTmp
!      double precision factCrossSection
!      double precision factMultiplicity
!      logical doResampling
      integer writeFirstIntPart
      integer particleListMode
      integer interactionCounter
      logical isFirstInt
      common/cxruext/writeFirstIntPart,particleListMode,
     +               interactionCounter,isFirstInt
#endif 
