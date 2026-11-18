#ifndef _resample_CommonBlockParticleSIBYLL_h_
#define _resample_CommonBlockParticleSIBYLL_h_

#include <ParticleBlockEntry.h>
#include <CommonBlockSIBYLL.h>

namespace resample {
  
  class CommonBlockWrapperSIBYLL;

  class CommonBlockParticleSIBYLL : public ParticleBlockEntry {
    
  public:
    friend class CommonBlockSIBYLL;
    friend class CommonBlockWrapperSIBYLL;

  private:
    CommonBlockParticleSIBYLL(const CommonBlockParticleSIBYLL&);
    CommonBlockParticleSIBYLL& operator=(const CommonBlockParticleSIBYLL&);
  
  public:
    virtual ~CommonBlockParticleSIBYLL() {}

    // interface implementation

    static double GetCharge(const int id);
    static std::string GetName(const int id);

    virtual bool IsGood() const;
    
    virtual double GetEnergy() const {return GetMomentum(3);}
    virtual double GetMass() const {return GetMomentum(4);}
    virtual double GetPx() const {return GetMomentum(0);}
    virtual double GetPy() const {return GetMomentum(1);}
    virtual double GetPz() const {return GetMomentum(2);}
    virtual int GetId() const {return fBlock->id[fIndex];}  

    virtual double GetCharge() const {return GetCharge(GetId());}  
    virtual std::string GetName() const {return GetName(GetId());}  

    virtual void SetEnergy(const double v) {GetMomentum(3)=v;}
    virtual void SetMass(const double v) {GetMomentum(4)=v;}
    virtual void SetPx(const double v) {GetMomentum(0)=v;}
    virtual void SetPy(const double v) {GetMomentum(1)=v;}
    virtual void SetPz(const double v) {GetMomentum(2)=v;}
    virtual void SetId(const int id) {fBlock->id[fIndex]=id; CheckId(id);}  
    void CheckId(const int id);	       
        
    CommonBlockParticleSIBYLL(CommonBlockSIBYLL* block, const int index); //, bool isLeading);
    
    virtual int GetIndex() const {return fIndex;}    
    virtual float& GetMomentum(int i) {return fBlock->momentum[i][fIndex];}    
    virtual const float& GetMomentum(int i)  const {return fBlock->momentum[i][fIndex];}
    void SetMomentum(int i, float v) {fBlock->momentum[i][fIndex]=v;}
    
    
  protected:
    mutable CommonBlockSIBYLL* fBlock;
    int fIndex;
  };

}

#endif
