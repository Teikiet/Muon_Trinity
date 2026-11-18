#ifndef _resample_CommonBlockWrapperCONEX_h_
#define _resample_CommonBlockWrapperCONEX_h_

#include <ParticleBlock.h>
#include <CommonBlockCONEX.h>

namespace resample {
  
  class CommonBlockParticleCONEX;
  
  class CommonBlockWrapperCONEX : public ParticleBlock {
    
  public:
    friend class CommonBlockCONEX;
    
  private:
    CommonBlockWrapperCONEX();
    CommonBlockWrapperCONEX(const CommonBlockWrapperCONEX&);
    const CommonBlockWrapperCONEX& operator=(const CommonBlockWrapperCONEX&);

  public:
    CommonBlockWrapperCONEX(CommonBlockCONEX* data, bool cms);
    CommonBlockWrapperCONEX(CommonBlockCONEX& data, bool cms);
    virtual ~CommonBlockWrapperCONEX() {}
    
  public: 
    // 'RUParticleBlock' interface implementation
    virtual ParticleBlockEntry& GetEntry(int i);
    virtual ParticleBlockEntry& operator[](int i) {return GetEntry(i);}
    virtual const ParticleBlockEntry& GetEntry(int i) const;
    virtual const ParticleBlockEntry& operator[](int i) const {return GetEntry(i);}
    virtual int Size() const {return fData->GetNumberEntries();}
    virtual void IncrementSize() {fData->number++; }
    virtual void DecrementSize() {fData->number--; }
    
    void Copy(int i1, int i2);
    
    CommonBlockCONEX& GetData() {return *fData;}
    const CommonBlockCONEX& GetData() const {return *fData;}
    
  protected:
    CommonBlockCONEX* fData;
  };

}

#endif
