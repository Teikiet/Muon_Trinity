#ifndef _CommonBlockSIBYLL_LAB_h_
#define _CommonBlockSIBYLL_LAB_h_

#include <CommonBlockCONEX.h>

/*
  Definition and mapping of CONEX secondary particle stack
  to individual nucelon-air interactions of CR nuclei
*/
namespace resample {
  class CommonBlockWrapperSIBYLL_LAB;
}

extern "C" {
  
  const static int gMaxNucleons = 200; // from conex.incnex

  class CommonBlockSIBYLL_LAB {
    
  public:
    friend class resample::CommonBlockWrapperSIBYLL_LAB;
    
  private:
    CommonBlockSIBYLL_LAB();
    CommonBlockSIBYLL_LAB(const CommonBlockSIBYLL_LAB&);
    const CommonBlockSIBYLL_LAB& operator=(const CommonBlockSIBYLL_LAB&);
    ~CommonBlockSIBYLL_LAB();
    
  public:
    int GetNumberOfInteractions() const {return numberOfInteractions;}
    int GetNumberEntries(const int id) const {return numSecondariesOfInteraction[id];}
    
  protected:
    /*
      common/ruNucCnt/idnucrct(mxptlxs),npnucs(200),imaxnuc
    */
    int id[gMxpltlxs];
    int numSecondariesOfInteraction[gMaxNucleons];
    int numberOfInteractions;
  };
} // extern "C"

#endif
