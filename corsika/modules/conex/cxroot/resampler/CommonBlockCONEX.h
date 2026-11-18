#ifndef _CommonBlockCONEX_h_
#define _CommonBlockCONEX_h_

/*
  Definition and mapping of the CONEX (also NEXUS) particle 
  common block/stack.
*/

namespace resample {
  class CommonBlockParticleCONEX;
  class CommonBlockWrapperCONEX;
  class CommonBlockWrapperSIBYLL_LAB;
}

extern "C" {
  
#warning Check this size with CONEX !
  static const int gMxpltlxs = 100000; /* from conex.incnex */
  
  class CommonBlockCONEX {
    
  public:
    friend class resample::CommonBlockParticleCONEX;
    friend class resample::CommonBlockWrapperSIBYLL_LAB;
    friend class resample::CommonBlockWrapperCONEX;
    
  private:
    CommonBlockCONEX();
    CommonBlockCONEX(const CommonBlockCONEX&);
    const CommonBlockCONEX& operator=(const CommonBlockCONEX&);
    ~CommonBlockCONEX();
    
  public:
    int GetNumberEntries() const {return number;}
    void SetNumberEntries(int n) {number=n;}
    
  protected:
    //public:
    
    /* 
       parameter (mxptlxs=100000)   !max nr of particles in conex MC particle list
       common/xscptl/xsptl(5,mxptlxs),xstivptl(2,mxptlxs)
       *,xsorptl(4,mxptlxs),ibptlxs(4,mxptlxs),ifrptlxs(2,mxptlxs)
       *,istptlxs(mxptlxs),jorptlxs(mxptlxs),idptlxs(mxptlxs)
       *,iorptlxs(mxptlxs),ityptlxs(mxptlxs),nptlxs
       
       c     nptlxs .......... current particle index (=particle number)
       c     idptlxs(i) ...... particle id
       c     xsptl(1,i) ..... x-component of particle momentum 
       c     xsptl(2,i) ..... y-component of particle momentum 
       c     xsptl(3,i) ..... z-component of particle momentum 
       c     xsptl(4,i) ..... particle energy 
       c     xsptl(5,i) ..... particle mass 
       c     iorptlxs(i) ..... particle number of father (if .le. 0 : no father) 
       c     jorptlxs(i) ..... particle number of mother (if .le. 0 : no mother)
       c     istptlxs(i) ..... status: 40 and 41 : Remnant
       c                             30 and 31 : Pomeron
       c                             20 and 21 : Parton
       c                             10 and 11 : Droplet
       c                             00 and 01 : Particle
       c                            last digit = 0 : last generation
       c                            last digit = 1 : not last generation
       c     xsorptl(1,i) ... x-component of formation point
       c     xsorptl(2,i) ... y-component of formation point
       c     xsorptl(3,i) ... z-component of formation point
       c     xsorptl(4,i) ... formation time
       c     xstivptl(1,i) ... formation time (always in the pp-cms!)
       c     xstivptl(2,i) ... destruction time (always in the pp-cms!)
       c     ityptlxs(i)  .... type of particles origin:
       c                         10-19: target
       c                         20-29: soft Pom
       c                         30-39: hard Pom 
       c                         40-49: projectile 
       c                         50: string, droplet
       
    */
    
    // xsptl(5,mxptlxs)
    double momentum[gMxpltlxs][5];  // px py pz E M
    // xstivptl(2,mxptlxs)
    double time[gMxpltlxs][2];      // formation desctuction (pp-cms)
    // xsorptl(4,mxptlxs)
    double formation[gMxpltlxs][4]; // x y z t
    // ibptlxs(4,mxptlxs)
    int ibptlxs[gMxpltlxs][4];      // unknown
    // ifrptlxs(2,mxptlxs)
    int ifrptlxs[gMxpltlxs][2];     // unknown
    // istptlxs(mxptlxs)
    int status[gMxpltlxs];          // status
    // jorptlxs(mxptlxs)
    int mother[gMxpltlxs];          // mother id
    // idptlxs(mxptlxs)            
    int id[gMxpltlxs];              // id
    // iorptlxs(mxptlxs)            
    int father[gMxpltlxs];          // father
    // ityptlxs(mxptlxs)
    int origin[gMxpltlxs];          // origin
    // nptlxs
    int number;                     // number of entries in block
    
  };
  
} // extern "C"

#endif
