#ifndef _resample_ParticleBlockEntry_h_
#define _resample_ParticleBlockEntry_h_

#include <string>
#include <iostream> // for dump
#include <iomanip> // for dump

namespace resample {

  class ParticleBlockEntry {
  public:
    
    ParticleBlockEntry() : fIdOfAntiParticle(0) {}
    virtual ~ParticleBlockEntry() {}
    
    virtual bool IsGood() const = 0;
    
    virtual double GetEnergy() const = 0;
    virtual double GetMass() const = 0;
    virtual double GetPx() const = 0;
    virtual double GetPy() const = 0;
    virtual double GetPz() const = 0;
    virtual int GetId() const = 0;

    virtual bool IsForward() const {return (GetPz()>0);}

    virtual double GetCharge() const = 0;
    virtual std::string GetName() const = 0;
    virtual int GetIdOfAntiParticle() const { return fIdOfAntiParticle; }
    
    virtual void SetEnergy(const double v)  = 0;
    virtual void SetMass(const double v)  = 0;
    virtual void SetPx(const double v)  = 0;
    virtual void SetPy(const double v)  = 0;
    virtual void SetPz(const double v)  = 0;
    virtual void SetId(const int v)  = 0;
    //    virtual void SetCharge(const int v)  = 0; // determined by Id
    
    virtual void Dump() const {
      const double E = GetEnergy();
      //const double P = sqrt(pow(p.GetMomentum(0),2)+pow(p.GetMomentum(1),2)+pow(p.GetMomentum(2),2));    
      std::cout
	<< " id=" << std::setw(5) << GetId()
	<< " (" << std::setw(11) << GetName() << ")"
	//<< " status=" << setw(2) << p.GetStatus()
	<< " E=" << std::setw(12) << E // p.GetMomentum(3)
	//<< " E_lab=" << setw(12) << Elab
	<< " m=" << std::setw(12) << std::setprecision(6) << GetMass()
	//<< " p=" << setprecision(2) << P
	<< std::endl;
    }
    
  protected:
    int fIdOfAntiParticle;
  };

}

#endif
