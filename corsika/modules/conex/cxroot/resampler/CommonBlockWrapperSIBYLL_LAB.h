#ifndef _resample_CommonBlockWrapperSIBYLL_LAB_h_
#define _resample_CommonBlockWrapperSIBYLL_LAB_h_

#include <ParticleBlock.h>
#include <CommonBlockWrapperCONEX.h>

#include <CommonBlockSIBYLL_LAB.h>
#include <CommonBlockCONEX.h>

#include <vector>

namespace resample {
  
  class CommonBlockWrapperSIBYLL_LAB : public CommonBlockWrapperCONEX {
    
  public:
    friend class CommonBlockSIBYLL_LAB;
    friend class CommonBlockCONEX;
    
  private:
    CommonBlockWrapperSIBYLL_LAB();
    CommonBlockWrapperSIBYLL_LAB(const CommonBlockWrapperSIBYLL_LAB&);
    const CommonBlockWrapperSIBYLL_LAB& operator=(const CommonBlockWrapperSIBYLL_LAB&);

  public:
    CommonBlockWrapperSIBYLL_LAB(CommonBlockCONEX& data, 
				 CommonBlockSIBYLL_LAB& mapping,
				 const int idInteraction,
				 const bool cms);
    virtual ~CommonBlockWrapperSIBYLL_LAB() {}
    
  public: 
    // 'ParticleBlock' interface implementation
    virtual ParticleBlockEntry& GetEntry(int i);
    virtual ParticleBlockEntry& operator[](int i) {return GetEntry(i);}
    virtual const ParticleBlockEntry& GetEntry(int i) const;
    virtual const ParticleBlockEntry& operator[](int i) const {return GetEntry(i);}
    virtual int Size() const;
    virtual void IncrementSize();
    virtual void DecrementSize();

    virtual void Delete(int i);
    virtual bool Duplicate(int i);
    
    virtual void Dump() const;

  protected:
    virtual void Copy(const int i1, const int i2);
        
  private:
    int FindIndex(const int i) const;
    int FindIndexReverse(const int i) const;

  protected:
    int fIdInteraction;
    CommonBlockSIBYLL_LAB* fMapping;
    
    mutable std::vector<int> fIndexMap;
  };

}

#endif
