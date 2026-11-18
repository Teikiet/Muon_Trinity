#ifndef _resample_CommonBlockParticleCRMC_h_
#define _resample_CommonBlockParticleCRMC_h_

#include <ParticleBlockEntry.h>
#include <CommonBlockCRMC.h>

namespace resample {
  
  class CommonBlockWrapperCRMC;
  
  class CommonBlockParticleCRMC : public ParticleBlockEntry {
    
  public:
    enum EParticleId {
      eGamma = 10,
      eElectron = 12,
      ePositron = -12,
      eProton = 1120,
      eProtonBar = -1120,
      eNeutron = 1220,
      eNeutronBar = -1220,
      ePi0 = 110,
      ePiP = 120,
      ePiM = -120,
      eKP = 130,
      eKM = -130,
      eKshort = -20, 
      eKlong = 20, 
      eLambda = 2130,
      eLambdaBar = -2130
    };
    
  public:
    friend class CommonBlockCRMC;
    friend class CommonBlockWrapperCRMC;
  
  private:
    CommonBlockParticleCRMC(const CommonBlockParticleCRMC&);
    CommonBlockParticleCRMC& operator=(const CommonBlockParticleCRMC&);
  
  public:
    virtual ~CommonBlockParticleCRMC() {}
    
    // interface implementation
    
    static double GetCharge(const int id);
    static std::string GetName(const int id);
    
    virtual bool IsGood() const;
    virtual double GetCharge() const {return GetCharge(GetId());}  
    virtual std::string GetName() const {return GetName(GetId());}
    
    virtual double GetEnergy() const {return GetMomentum(3);}
    virtual double GetMass() const {return GetMomentum(4);}
    virtual double GetPx() const {return GetMomentum(0);}
    virtual double GetPy() const {return GetMomentum(1);}
    virtual double GetPz() const {return GetMomentum(2);}
    virtual int GetId() const {return fBlock->idptl[fIndex];}
    
    virtual void SetEnergy(const double v) {GetMomentum(3) = v;}
    virtual void SetMass(const double v) {GetMomentum(4) = v;}
    virtual void SetPx(const double v) {GetMomentum(0) = v;}
    virtual void SetPy(const double v) {GetMomentum(1) = v;}
    virtual void SetPz(const double v) {GetMomentum(2) = v;}
    virtual void SetId(const int id) {fBlock->idptl[fIndex] = id; CheckId(id);}
    
    void CheckId(const int id);	       
    
    CommonBlockParticleCRMC(CommonBlockCRMC* block, int index); //, bool isLeading);
    
    virtual int GetIndex() const {return fIndex;}
    
    virtual double& GetMomentum(int i) {return fBlock->pptl[fIndex][i];}
    virtual double& GetFormationTime() {return fBlock->tivptl[fIndex][0];}
    virtual double& GetDestructionTime() {return fBlock->tivptl[fIndex][1];}
    virtual double& GetFormationPoint(int i) {return fBlock->xorptl[fIndex][i];}
    //virtual int& GetIbtlxs(int i) {return fBlock->ibptlxs[fIndex][i];}
    //virtual int& GetIfrptlxs(int i) {return fBlock->ifrptlxs[fIndex][i];}
    virtual int& GetStatus() {return fBlock->istptl[fIndex];}
    virtual int& GetMotherId() {return fBlock->jorptl[fIndex];}
    virtual int& GetFatherId() {return fBlock->iorptl[fIndex];}
    //virtual int& GetId() {std::cout << " GetId: block& " << fBlock << std::endl; return fBlock->id[fIndex];}
    virtual int& GetId() {return fBlock->idptl[fIndex];}
    //virtual int& GetOrigin() {return fBlock->origin[fIndex];}
    
    virtual const double& GetMomentum(int i)  const {return fBlock->pptl[fIndex][i];}
    virtual const double& GetFormationTime()  const {return fBlock->tivptl[fIndex][0];}
    virtual const double& GetDestructionTime()  const {return fBlock->tivptl[fIndex][1];}
    virtual const double& GetFormationPoint(int i)  const {return fBlock->xorptl[fIndex][i];}
    //virtual const int& GetIbtlxs(int i)  const {return fBlock->ibptlxs[fIndex][i];}
    //virtual const int& GetIfrptlxs(int i)  const {return fBlock->ifrptlxs[fIndex][i];}
    virtual const int& GetStatus()  const {return fBlock->istptl[fIndex];}
    virtual const int& GetMotherId()  const {return fBlock->jorptl[fIndex];}
    virtual const int& GetFatherId()  const {return fBlock->iorptl[fIndex];}
    //virtual const int& GetOrigin()  const {return fBlock->origin[fIndex];}
    
    void SetMomentum(int i, double v) {fBlock->pptl[fIndex][i]=v;}
    void SetFormationTime(double v) {fBlock->tivptl[fIndex][0]=v;}
    void SetDestructionTime(double v) {fBlock->tivptl[fIndex][1]=v;}
    void SetFormationPoint(int i, double v) {fBlock->xorptl[fIndex][i]=v;}
    //void SetIbtlxs(int i, int v) {fBlock->ibptlxs[fIndex][i]=v;}
    //void SetIfrptlxs(int i, int v) {fBlock->ifrptlxs[fIndex][i]=v;}
    void SetStatus(int v) {fBlock->istptl[fIndex]=v;}
    void SetMotherId(int v) {fBlock->jorptl[fIndex]=v;}
    void SetFatherId(int v) {fBlock->iorptl[fIndex]=v;}
    //void SetOrigin(int v) {fBlock->origin[fIndex]=v;}
    
    
  protected:
    mutable CommonBlockCRMC* fBlock;
    int fIndex;
  };

}

#endif
