// Settle the historical discrepancy: which calculation does each number represent?
// At guard_b's worst ray, compare the three candidate definitions of the
// "table-equivalent" tau against the converged reference.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static double alphaCB(double lam,double h){ double const dh=50.0; auto const hr=g_ab.getHeightRange();
  double const lo=std::max(hr.first,h-dh*0.5), hi=std::min(hr.second,h+dh*0.5);
  double const a=hi-lo; if(a<=0.0) return 0.0;
  return std::max(0.0,(g_ab.getOpticalDepth(lam,hi)-g_ab.getOpticalDepth(lam,lo))/a); }
static double nP(double r,double f){ return 1.0+N0*std::exp(-(r-R)/HSC)*f; }
static double Aint(double b,double f,double lam,double rTop,double dr){ double acc=0,r=R;
  while(r<rTop-1e-12){ double const st=std::min(dr,rTop-r); double const rm=r+0.5*st;
    double const s=b/(nP(rm,f)*rm); if(s>=1.0) return acc;
    acc+=alphaCB(lam,rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st; r+=st; } return acc; }
int main(int argc,char** argv){ g_ab.load(argv[1]);
  double const h=1500.0, ratio=0.9995, lam=672.5;
  std::printf("guard_b worst ray: h=%.0f ratio=%.4f lam=%.1f\n", h, ratio, lam);
  for(double fTrue : {1.0347, 0.9713, 1.0}){
    double const bCritTrue = nP(R,fTrue)*R;
    double const bTrue = ratio*bCritTrue;                       // photon's true impact parameter
    double const bCritTab = nP(R,1.0)*R;                        // table built at f=1
    double const muObs = bTrue/bCritTrue;                       // the photon's observation-level mu
    double const bTab = muObs*bCritTab;                         // what a mu-parameterised table READS
    double const tR  = Aint(bTrue, fTrue, lam, R+h, 0.5);       // converged reference
    double const tA1 = Aint(bTrue, 1.0,   lam, R+h, 5.0);       // guard_b: true b + f=1 profile
    double const tA2 = Aint(bTab,  1.0,   lam, R+h, 5.0);       // realistic: table's own b + f=1 profile
    std::printf("  f=%.4f  tR=%.6f  tA(trueb)=%.6f |dT|=%.3e   tA(tabb)=%.6f |dT|=%.3e   b: true=%.3f tab=%.3f\n",
      fTrue, tR, tA1, std::abs(std::exp(-tA1)-std::exp(-tR)), tA2, std::abs(std::exp(-tA2)-std::exp(-tR)),
      bTrue, bTab);
  }
  return 0; }
