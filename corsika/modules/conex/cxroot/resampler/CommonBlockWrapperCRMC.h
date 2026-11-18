#ifndef _resample_CommonBlockWrapperCRMC_h_
#define _resample_CommonBlockWrapperCRMC_h_

#include <ParticleBlock.h>
#include <CommonBlockCRMC.h>

namespace resample {
  
  class CommonBlockParticleCRMC;
  
  class CommonBlockWrapperCRMC : public ParticleBlock {
    
  public:
    friend class CommonBlockCRMC;
    
  private:
    CommonBlockWrapperCRMC();
    CommonBlockWrapperCRMC(const CommonBlockWrapperCRMC&);
    const CommonBlockWrapperCRMC& operator=(const CommonBlockWrapperCRMC&);

  public:
    CommonBlockWrapperCRMC(CommonBlockCRMC* data, bool cms);
    CommonBlockWrapperCRMC(CommonBlockCRMC& data, bool cms);
    virtual ~CommonBlockWrapperCRMC() {}
    
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
    
    CommonBlockCRMC& GetData() {return *fData;}
    const CommonBlockCRMC& GetData() const {return *fData;}
    
  protected:
    CommonBlockCRMC* fData;
  };

}

#endif
