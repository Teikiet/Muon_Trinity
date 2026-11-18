#ifndef _resample_CommonBlockWrapperSIBYLL_h_
#define _resample_CommonBlockWrapperSIBYLL_h_

#include <ParticleBlock.h>
#include <CommonBlockSIBYLL.h>

namespace resample {
  
  class CommonBlockParticleSIBYLL;
  
  class CommonBlockWrapperSIBYLL : public ParticleBlock {
    
  public:
    friend class CommonBlockSIBYLL;

  private:
    CommonBlockWrapperSIBYLL();
    CommonBlockWrapperSIBYLL(const CommonBlockWrapperSIBYLL&);
    const CommonBlockWrapperSIBYLL& operator=(const CommonBlockWrapperSIBYLL&);
    
  public:
    CommonBlockWrapperSIBYLL(CommonBlockSIBYLL* data, bool cms);
    CommonBlockWrapperSIBYLL(CommonBlockSIBYLL& data, bool cms);
    virtual ~CommonBlockWrapperSIBYLL() {}
    
  public: 
    // 'RUParticleBlock' interface implementation
    
    virtual ParticleBlockEntry& GetEntry(int i);
    virtual ParticleBlockEntry& operator[](int i) {return GetEntry(i);}
    virtual const ParticleBlockEntry& GetEntry(int i) const;
    virtual const ParticleBlockEntry& operator[](int i) const {return GetEntry(i);}
    virtual int Size() const {return fData->GetNumberEntries();}
    virtual void IncrementSize() {fData->number++;}
    virtual void DecrementSize() {fData->number--;}
    
    void Copy(const int i1, const int i2);
  
    CommonBlockSIBYLL& GetData() {return *fData;}
    const CommonBlockSIBYLL& GetData() const {return *fData;}
  
  protected:
    CommonBlockSIBYLL* fData;
  };
}

#endif
