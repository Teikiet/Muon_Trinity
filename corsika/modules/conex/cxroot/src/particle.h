#ifndef _cxroot_particle_h_
#define _cxroot_particle_h_

#include <TObject.h>
#include <TString.h>
#include <string>
#include <cmath>

class particle : public TObject {
public: 
  particle()
    : id(0), name("unknow"), px(0), py(0), pz(0), E(0), m(0) {}
  particle(int ID, const std::string& n, double Px, double Py, double Pz, double e, double M)
    : id(ID), name(n.c_str()), px(Px), py(Py), pz(Pz), E(e), m(M) {}
  double GetP() const { return sqrt(px*px + py*py + pz*pz); }
  int id;
  TString name;
  double px;
  double py;
  double pz;
  double E;
  double m;
  
  particle GetBoosted(const double yboost, const bool toLabSystem) const {
    /*
c-----------------------------------------------------------------------
      subroutine cxutlob5(yboost,x1,x2,x3,x4,x5)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      amt=x5**2+x1**2+x2**2
      if(amt.gt.0.d0)then
        amt=sqrt(x5**2+x1**2+x2**2)
      else
        amt=max(1d-15,sqrt(abs((x4+x3)*(x4-x3))))
      endif
      y=sign(1.d0,x3)*log((x4+abs(x3))/amt)
      y=y-yboost
      x4=amt*cosh(y)
      x3=amt*sinh(y)
      return
      end
     */    
    const bool cms = toLabSystem;
    double amt = m*m + px*px + py*py;
    if (amt>0) 
      amt = sqrt(amt);
    else 
      amt = std::max(1.e-15, sqrt(std::abs((E+pz)*(E-pz))));
    const double y = ( pz<0 ? -1.0 : 1.0) * log((E+std::abs(pz))/amt) - 
      (cms ? -yboost : yboost );
    const double pz_boost = amt*sinh(y);
    const double E_boost = amt*cosh(y);
    return particle(id, name.Data(), px, py, pz_boost, E_boost, m);
  }
  
  ClassDef(particle, 1);
};

#endif
