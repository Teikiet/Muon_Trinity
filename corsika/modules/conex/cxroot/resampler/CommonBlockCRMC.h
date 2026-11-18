#ifndef _CommonBlockCRMC_h_
#define _CommonBlockCRMC_h_

/*
  Definition and mapping of the CRMC (also NEXUS) particle 
  common block/stack.
*/

namespace resample {
  class CommonBlockParticleCRMC;
  class CommonBlockWrapperCRMC;
  class CommonBlockWrapperSIBYLL_LAB;
}

extern "C" {
  
  static const int gMxpltlxs = 300000; /* from epos.inc */
  
  class CommonBlockCRMC /* : public DataAccess (this is tricky) */ {
    
  public:
    friend class resample::CommonBlockParticleCRMC;
    friend class resample::CommonBlockWrapperCRMC;
    
  private:
    CommonBlockCRMC();
    CommonBlockCRMC(const CommonBlockCRMC&);
    const CommonBlockCRMC& operator=(const CommonBlockCRMC&);
    ~CommonBlockCRMC();
    
  public:
    int GetNumberEntries() const {return number;}
    void SetNumberEntries(int n) {number=n;}
    
    void Copy(const int from, const int to) {
      idptl[to] = idptl[from]; //  ...... particle id      
      pptl[to][0] = pptl[from][0]; 
      pptl[to][1] = pptl[from][1];  
      pptl[to][2] = pptl[from][2]; 
      pptl[to][3] = pptl[from][3]; 
      pptl[to][4] = pptl[from][4]; 
      iorptl[to] = iorptl[from];
      jorptl[to] = jorptl[from]; 
      istptl[to] = istptl[from]; 
      xorptl[to][0] = xorptl[from][0];
      xorptl[to][1] = xorptl[from][1];
      xorptl[to][2] = xorptl[from][2];
      xorptl[to][3] = xorptl[from][3];
      tivptl[to][0] = tivptl[from][0]; 
      tivptl[to][1] = tivptl[from][1]; 
      ityptl[to] = ityptl[from];
      itsptl[to] = itsptl[from];
    }

  protected:
    //public:
    
    /* 
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list
      
      real        pptl,tivptl,xorptl
      integer     nptl,iorptl,idptl,istptl,ifrptl,jorptl,ibptl,ityptl
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real         ekievt
      integer      itsptl
      common/c1ptl/ekievt,itsptl(mxptl)
     
c     nptl .......... current particle index (=number of ptls stored)
c     idptl(i) ...... particle id
c     pptl(1,i) ..... x-component of particle momentum 
c     pptl(2,i) ..... y-component of particle momentum 
c     pptl(3,i) ..... z-component of particle momentum 
c     pptl(4,i) ..... particle energy 
c     pptl(5,i) ..... particle mass 
c     iorptl(i) ..... particle number of father (if .le. 0 : no father) 
c     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
c     istptl(i) ..... status: 40 and 41 : Remnant
c                             30 and 31 : Pomeron
c                             20 and 21 : Parton
c                             10 and 11 : Droplet
c                             00 and 01 : Particle
c                            last digit = 0 : last generation
c                            last digit = 1 : not last generation
c     xorptl(1,i) ... x-component of formation point
c     xorptl(2,i) ... y-component of formation point
c     xorptl(3,i) ... z-component of formation point
c     xorptl(4,i) ... formation time
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)
c     ityptl(i)  .... type of particles origin:
c                         10-19: target
c                         20-29: soft Pom
c                         30-39: hard Pom 
c                         40-49: projectile 
c                         50: string, droplet
c     itsptl(i) ..... string type of particles origin (if string)  
       
    */

    // nptl
    int number;                     // number of entries in block

    int idptl[gMxpltlxs]; //  ...... particle id
    /*
c     pptl(1,i) ..... x-component of particle momentum 
c     pptl(2,i) ..... y-component of particle momentum 
c     pptl(3,i) ..... z-component of particle momentum 
c     pptl(4,i) ..... particle energy 
c     pptl(5,i) ..... particle mass 
    */
    double pptl[gMxpltlxs][5]; 
    int iorptl[gMxpltlxs]; // particle number of father (if .le. 0 : no father) 
    int jorptl[gMxpltlxs]; // particle number of mother (if .le. 0 : no mother)
    int istptl[gMxpltlxs]; // status: 40 and 41 : Remnant
    //  c                             30 and 31 : Pomeron
    //  c                             20 and 21 : Parton
    //  c                             10 and 11 : Droplet
    //  c                             00 and 01 : Particle
    //  c                            last digit = 0 : last generation
    //  c                            last digit = 1 : not last generation
    /*
c     xorptl(1,i) ... x-component of formation point
c     xorptl(2,i) ... y-component of formation point
c     xorptl(3,i) ... z-component of formation point
c     xorptl(4,i) ... formation time
    */
    double xorptl[gMxpltlxs][4];
    /*
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)
    */
    double tivptl[gMxpltlxs][2]; 
    int ityptl[gMxpltlxs]; // type of particles origin:
    //c                         10-19: target
    //c                         20-29: soft Pom
    //c                         30-39: hard Pom 
    //c                         40-49: projectile 
    //c                         50: string, droplet
    int itsptl[gMxpltlxs]; // string type of particles origin (if string)     
  };
  
} // extern "C"

#endif
