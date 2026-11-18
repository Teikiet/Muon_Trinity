#ifndef _resample_ParticleBlock_h_
#define _resample_ParticleBlock_h_

#include <SummaryInfo.h>
#include <map>

namespace resample {
  
  class ParticleBlockEntry;
  
  class ParticleBlock {
  public:
    enum ELeading {
      eBackward,
      eForward,
      eBoth
    };
    
    enum EExcludeLeading {
      eExcludeLeading,
      eIncludeLeading
    };
    
    enum EIdMode {
      eExact,
      eTypes
    };
    
  public:
    ParticleBlock(int capacity, bool isCMS) 
      : fMaxParticleCapacity(capacity), fIsCMS(isCMS), fLeadingParticleForwardIndex(-1), 
	fLeadingParticleBackwardIndex(-1), fParticleBlockSizeOverflowCounter(0) {}
    virtual ~ParticleBlock();
    
    virtual std::map<int, SummaryInfo> MakeSummary(const EIdMode idMode,
						   const EExcludeLeading excludeLeading) const;
    virtual bool IsLeading(int i, ELeading where=eForward) const;
    virtual int GetLeadingParticleIndex(ELeading where=eForward) const;
    
    virtual ParticleBlockEntry& GetEntry(int i) = 0;
    virtual ParticleBlockEntry& operator[](int i) = 0;
    virtual const ParticleBlockEntry& GetEntry(int i) const = 0;
    virtual const ParticleBlockEntry& operator[](int i) const = 0;
    virtual int Size() const = 0;
    virtual void IncrementSize() = 0;
    virtual void DecrementSize() = 0;
    
    virtual void Delete(int i);
    virtual bool Duplicate(int i);

    virtual bool IsCenterOfMassSystem() const {return fIsCMS;}
    virtual bool IsLaboratorySystem() const {return !fIsCMS;}
    virtual int GetMaxCapacity() const {return fMaxParticleCapacity;}

    virtual void IdentifyLeadingParticle(ELeading where=eForward) const;
    
    virtual void Dump() const; // for debug
    
  protected:
    virtual void Copy(const int i1, const int i2) = 0;
    void SetLeadingParticleIndex(const int index, ELeading where=eForward) const;
    
  protected:
    int fMaxParticleCapacity;
    bool fIsCMS;
    
    mutable int fLeadingParticleForwardIndex; // cached access
    mutable int fLeadingParticleBackwardIndex; // cached access

  protected:
    mutable std::map<int, ParticleBlockEntry*> fParticles; // cached access

    int fParticleBlockSizeOverflowCounter;
  };
}

#endif
