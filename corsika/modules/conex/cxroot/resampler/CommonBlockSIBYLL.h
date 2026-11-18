#ifndef _CommonBlockSIBYLL_h_
#define _CommonBlockSIBYLL_h_

/*
  Definition and mapping of the SIBYLL particle 
  common block/stack.
*/

namespace resample {
  class CommonBlockParticleSIBYLL;
  class CommonBlockWrapperSIBYLL;
}

extern "C" {
  
  static const int gMaxParticlesSIBYLL = 8000; /* from conex_sub.f  */
  
  class CommonBlockSIBYLL {
    
  public:
    friend class resample::CommonBlockParticleSIBYLL;
    friend class resample::CommonBlockWrapperSIBYLL;
    
  private:
    CommonBlockSIBYLL();
    CommonBlockSIBYLL(const CommonBlockSIBYLL&);
    const CommonBlockSIBYLL& operator=(const CommonBlockSIBYLL&);
    ~CommonBlockSIBYLL();
    
  public:
    int GetNumberEntries() const {return number;}
    void SetNumberEntries(int n) {number=n;}
    
  protected:
    /*
      COMMON /S_PLIST/ NP, P(8000,5), LLIST(8000)
    */
    int number;
    float momentum[5][gMaxParticlesSIBYLL];  // px py pz E M
    int id[gMaxParticlesSIBYLL];
        
  };
  
} // extern "C"

#endif
