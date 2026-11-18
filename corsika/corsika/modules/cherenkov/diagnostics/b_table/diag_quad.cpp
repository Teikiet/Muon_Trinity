// Decisive check: was guard_b's 2.2e-3 its OWN quadrature error rather than the table's?
// Compare, at the failing ray, the UNSplit dr=5 integral against the split one.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static std::vector<double> g_hk;
static double alphaCB(double lam, double h) {
  double const dh = 50.0; auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h-dh*0.5), hi = std::min(hr.second, h+dh*0.5);
  double const a = hi-lo; if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam,hi) - g_ab.getOpticalDepth(lam,lo))/a);
}
static double nP(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
static double unsplit(double b,double f,double lam,double rTop,double dr){ double acc=0,r=R;
  while(r<rTop-1e-12){double const st=std::min(dr,rTop-r);double const rm=r+0.5*st;
    double const s=b/(nP(rm,f)*rm); if(s>=1.0)return acc;
    acc+=alphaCB(lam,rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st; r+=st;} return acc; }
static double split(double b,double f,double lam,double rTop,double dr){
  std::vector<double> cuts; for(double h:g_hk){double c1=R+h-25.0,c2=R+h+25.0;
    if(c1>R&&c1<rTop)cuts.push_back(c1); if(c2>R&&c2<rTop)cuts.push_back(c2);} cuts.push_back(rTop);
  std::sort(cuts.begin(),cuts.end()); double acc=0,r=R;
  for(double c:cuts){ while(r<c-1e-12){double const st=std::min(dr,c-r);double const rm=r+0.5*st;
    double const s=b/(nP(rm,f)*rm); if(s>=1.0)return acc;
    acc+=alphaCB(lam,rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st; r+=st;} } return acc; }
int main(int argc,char** argv){ g_ab.load(argv[1]);
  auto const hr=g_ab.getHeightRange(); for(double h=hr.first;h<=hr.second+1e-9;h+=1000.0) g_hk.push_back(h);
  struct C{double h,ratio,f,lam;}; C const cs[]={{1500.0,0.9995,1.0347,672.5},{1500.0,0.9950,1.0347,672.5},
    {1500.0,0.9995,1.0,672.5},{20000.0,0.9995,1.0347,672.5}};
  std::printf("ray                        unsplit(5)   split(5)    unsplit(0.01)  |dT unsplit5 vs split5|"); puts("");
  for(C const& c:cs){ double const b=c.ratio*(nP(R,c.f)*R);
    double const u5=unsplit(b,c.f,c.lam,R+c.h,5.0), s5=split(b,c.f,c.lam,R+c.h,5.0);
    double const u01=unsplit(b,c.f,c.lam,R+c.h,0.01);
    std::printf("h=%5.0f ratio=%.4f f=%.4f lam=%.1f  %.9f  %.9f  %.9f   %.3e",
      c.h,c.ratio,c.f,c.lam,u5,s5,u01,std::abs(std::exp(-u5)-std::exp(-s5))); puts(""); }
  return 0; }
