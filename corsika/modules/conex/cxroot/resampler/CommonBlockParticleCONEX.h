#ifndef _resample_CommonBlockParticleCONEX_h_
#define _resample_CommonBlockParticleCONEX_h_

#include <ParticleBlockEntry.h>
#include <CommonBlockCONEX.h>

namespace resample {
  
  class CommonBlockWrapperCONEX;
  class CommonBlockWrapperSIBYLL_LAB;

  class CommonBlockParticleCONEX : public ParticleBlockEntry {

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
      eLambdaBar = -2130,
      eRho0 = 111,
      eRhoP = 121,
      eRhoM = -121
    };
    
  public:
    friend class CommonBlockCONEX;
    friend class CommonBlockWrapperCONEX;
    friend class CommonBlockWrapperSIBYLL_LAB;

  private:
    CommonBlockParticleCONEX(const CommonBlockParticleCONEX&);
    CommonBlockParticleCONEX& operator=(const CommonBlockParticleCONEX&);
  
  public:
    virtual ~CommonBlockParticleCONEX() {}
    
    // interface implementation

    static double GetCharge(const int id);
    static std::string GetName(const int id);
    static double GetMass(const int id);
    static resample::CommonBlockParticleCONEX::EParticleId GetIdOfAntiParticle(const resample::CommonBlockParticleCONEX::EParticleId id);
    
    virtual bool IsGood() const;
    virtual double GetCharge() const {return GetCharge(GetId());}  
    virtual std::string GetName() const {return GetName(GetId());}
    virtual double GetMassPID() const {return GetMass(GetId());}
    
    virtual double GetEnergy() const {return GetMomentum(3);}
    virtual double GetMass() const {return GetMomentum(4);}
    virtual double GetPx() const {return GetMomentum(0);}
    virtual double GetPy() const {return GetMomentum(1);}
    virtual double GetPz() const {return GetMomentum(2);}
    virtual int GetId() const {return fBlock->id[fIndex];}
    
    virtual void SetEnergy(const double v) {GetMomentum(3) = v;}
    virtual void SetMass(const double v) {GetMomentum(4) = v;}
    virtual void SetPx(const double v) {GetMomentum(0) = v;}
    virtual void SetPy(const double v) {GetMomentum(1) = v;}
    virtual void SetPz(const double v) {GetMomentum(2) = v;}
    virtual void SetId(const int id) {fBlock->id[fIndex] = id; CheckId(id);}
    void CheckId(const int id);	       
    
    CommonBlockParticleCONEX(CommonBlockCONEX* block, int index); //, bool isLeading);
    
    virtual int GetIndex() const {return fIndex;}
    
    virtual double& GetMomentum(int i) {return fBlock->momentum[fIndex][i];}
    virtual double& GetFormationTime() {return fBlock->time[fIndex][0];}
    virtual double& GetDestructionTime() {return fBlock->time[fIndex][1];}
    virtual double& GetFormationPoint(int i) {return fBlock->formation[fIndex][i];}
    virtual int& GetIbtlxs(int i) {return fBlock->ibptlxs[fIndex][i];}
    virtual int& GetIfrptlxs(int i) {return fBlock->ifrptlxs[fIndex][i];}
    virtual int& GetStatus() {return fBlock->status[fIndex];}
    virtual int& GetMotherId() {return fBlock->mother[fIndex];}
    //virtual int& GetId() {std::cout << " GetId: block& " << fBlock << std::endl; return fBlock->id[fIndex];}
    virtual int& GetId() {return fBlock->id[fIndex];}
    virtual int& GetFatherId() {return fBlock->father[fIndex];}
    virtual int& GetOrigin() {return fBlock->origin[fIndex];}
    
    virtual const double& GetMomentum(int i)  const {return fBlock->momentum[fIndex][i];}
    virtual const double& GetFormationTime()  const {return fBlock->time[fIndex][0];}
    virtual const double& GetDestructionTime()  const {return fBlock->time[fIndex][1];}
    virtual const double& GetFormationPoint(int i)  const {return fBlock->formation[fIndex][i];}
    virtual const int& GetIbtlxs(int i)  const {return fBlock->ibptlxs[fIndex][i];}
    virtual const int& GetIfrptlxs(int i)  const {return fBlock->ifrptlxs[fIndex][i];}
    virtual const int& GetStatus()  const {return fBlock->status[fIndex];}
    virtual const int& GetMotherId()  const {return fBlock->mother[fIndex];}
    virtual const int& GetFatherId()  const {return fBlock->father[fIndex];}
    virtual const int& GetOrigin()  const {return fBlock->origin[fIndex];}
    
    void SetMomentum(int i, double v) {fBlock->momentum[fIndex][i]=v;}
    void SetFormationTime(double v) {fBlock->time[fIndex][0]=v;}
    void SetDestructionTime(double v) {fBlock->time[fIndex][1]=v;}
    void SetFormationPoint(int i, double v) {fBlock->formation[fIndex][i]=v;}
    void SetIbtlxs(int i, int v) {fBlock->ibptlxs[fIndex][i]=v;}
    void SetIfrptlxs(int i, int v) {fBlock->ifrptlxs[fIndex][i]=v;}
    void SetStatus(int v) {fBlock->status[fIndex]=v;}
    void SetMotherId(int v) {fBlock->mother[fIndex]=v;}
    void SetFatherId(int v) {fBlock->father[fIndex]=v;}
    void SetOrigin(int v) {fBlock->origin[fIndex]=v;}


    static const double fgMassProton;
    static const double fgMassNeutron;
    static const double fgMassLambda;
    static const double fgMassPionPM;
    static const double fgMassPion0;
    static const double fgMassKaon0;
    static const double fgMassKaonPM;
    static const double fgMassEta0;
    static const double fgMassElectron;
    static const double fgMassMuon;
    static const double fgMassTau;

    
    
  protected:
    mutable CommonBlockCONEX* fBlock;
    int fIndex;
  };

}

#endif
