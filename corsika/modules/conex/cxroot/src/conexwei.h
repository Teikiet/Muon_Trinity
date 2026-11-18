      integer mxn1,mxn2,mxppj,iemin0,iemax0,iemin1,iemax1,n1max,n2max
      parameter(mxn1=5) 
      parameter(mxn2=11) 
      parameter(mxppj=262) 
      real ppj,exmin0
      common /cppj1/ ppj(mxppj,mxppj,mxn2,mxn1)
      common /cppj2/ exmin0,iemin0,iemax0,iemin1,iemax1,n1max,n2max
#if __MC3D__ || __CXLATCE__
      real p2j!,p4j
      common /cppj3/ p2j(mxppj,mxppj,mxn2,mxn1)
c      common /cppj4/ p4j(mxppj,mxppj,mxn2,mxn1)
#endif
      character ppjver*20

