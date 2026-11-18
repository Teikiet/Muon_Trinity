c-----------------------------------------------------------------------------
c                       Include file for conex e/m part
c
c Last modification 10.10.2005 by T. Pierog (ver 3.059)
c Author R. Engel, N. Kalmykov, S. Ostaptchenko, T. Pierog, K. Werner
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------------
c     ee(j) ......... discrete energy
c     zz(k) ......... discrete depth
c     aeo(j,k) ...... initial number of electrons at energy j and depth k
c     ago(j,k) ...... initial number of gammas at energy j and depth k
c     ae00(j,k) ..... number of electrons at energy j and depth k
c     ag00(j,k) ..... number of gammas at energy j and depth k
c     ele0(k) ....... number of electrons at depth k
c     gam0(k) ....... number of gammas at  depth k
c     wgg,wge,sg .... cross sections (gammas)
c     weg,wee,se .... cross section tables (electrons)
c-----------------------------------------------------------------------

      INTEGER minE,maxE,maxZ,imaxE0,maximumE,maximumZ,jminZ0
     *,kfirst,modk,lowE,LowZ
      DOUBLE PRECISION AG,AE,AP,AEM,AGM,APM,DLT,DLTP
     *,Suu,Suv,Svu,Svv,Sww,Svw,Swu,Swv,WGG,WGE,WEG,WEE,WGP,WPP,WPE,WPG
     *,amc2,Eo,radle,Cem,ZZo,dZZ,emdecade,c2em,SFE,SFG,SFP,sf2had
     *,ZZEM,EEEM,dethe,dethg,distzem,wgh,wgm,Svh,Svm,wsf2noint,dethp
     *,betheb


#ifdef __SAVEMEMO__
      PARAMETER (maximumE=281,maximumZ=201)
#else
      PARAMETER (maximumE=281,maximumZ=3701)
#endif

      common/cxtempo/AG(maximumE,maximumZ),AE(maximumE,maximumZ)
     *,AP(maximumE,maximumZ)
      common/cxutil/DLT(maximumE,maximumZ),DLTP(maximumE,maximumZ)
      COMMON /CXAREA12/Suu(maximumE,maximumE),Suv(maximumE,maximumE)
     *,Svu(maximumE,maximumE),Svv(maximumE,maximumE)
     *,Svw(maximumE,maximumE),Sww(maximumE,maximumE)
     *,Swu(maximumE,maximumE),Swv(maximumE,maximumE)
      COMMON /CXAREA13/ WGG(maximumE,maximumE),WGE(maximumE,maximumE)
      COMMON /CXAREA14/ WEG(maximumE,maximumE),WEE(maximumE,maximumE)
     *,WGP(maximumE,maximumE),WPP(maximumE,maximumE)
     *,WPE(maximumE,maximumE),WPG(maximumE,maximumE)
      COMMON /CXAREA15/ amc2,Eo,radle
      common /CXAREA16/ Cem,ZZo,dZZ,emdecade,c2em
     *,minE,maxE,maxZ,imaxE0,jminZ0
      common /CXAREA17/ lowE,lowZ
      common/cxsource/SFE(maximumE,maximumZ),SFG(maximumE,maximumZ)
     *,SFP(maximumE,maximumZ)

      common/cxdiscret/ZZEM(maximumZ),EEEM(maximumE)
      common/cxcoutp/kfirst,modk
      common/cxdethe/dethe(maximume),dethg(maximume),distzem(maximumZ)
     *,dethp(maximume),betheb(2,maximume)      !so010304  !tp070503
      common/cxphotmu/wgh(maximume),wgm(maximume)      !tp071204
     *,Svh(maximumE),Svm(maximumE)
      common/cxhalfsource/sf2had(3,maximume),wsf2noint(3,maximume)



#if __MC3D__ || __CXLATCE__

c moment part

      integer lleng,maxoep,maximom,n4m,ISI,ISIr,n4mreal,i4A
      double precision ASC,ASCr,AF4m,AAEM
     *,AAGM,AAPM,Drho,source3d,sef
      logical lcorrp


      PARAMETER(lleng=2)
#if __CXLATCE__
      PARAMETER(maxoep=2)
      PARAMETER (maximom=maxoep+2,n4m=9)
#else
      PARAMETER(maxoep=0)
      PARAMETER (maximom=2,n4m=3)
#endif

      COMMON /CXAREA20/ sef(3,3,maximumE,maximumE)

      common/eqsource2/
     *   ASC(0:n4m, 0:lleng),
     *   ASCr(0:n4m, 0:lleng),
     *   ISI(0:n4m, 0:lleng),
     *   ISIr(0:n4m, 0:lleng,3)
     *        ,n4mreal

      common/cx4moments/ AF4m(0:n4m,maximume,0:3,3),i4A(0:n4m,5)


      common/cxmoments/AEM(0:maximom,maximumE,maximumZ)
     *,AGM(0:maximom,maximumE,maximumZ)
     *,APM(0:maximom,maximumE,maximumZ)
     *,AAEM(0:maximom,maximumE,maximumZ)
     *,AAGM(0:maximom,maximumE,maximumZ)
     *,AAPM(0:maximom,maximumE,maximumZ), Drho

       common/cxprectrl/lcorrp

       common/cxsource2/source3d(0:n4m,maximume,maximumZ,3)

#else

      common/cxmoments/AEM(0:0,maximumE,maximumZ)
     *,AGM(0:0,maximumE,maximumZ),APM(0:0,maximumE,maximumZ)

#endif
