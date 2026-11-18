// Remaining coverage diagnostics for the cumulative-b absorption design.
// Production untouched. Decisions mirror the PROPOSED eligibility, not wired code.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static double nPhase(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
// proposed eligibility decision for the absorption (cumulative-b) table
enum class D { Accept, FallbackLegacy, RejectUpgoing, RejectBeyondCrit, RejectBelowObs,
               RejectAboveTop, RejectLambdaOut, RejectFOut, AbsorptionOff };
static char const* dn(D d) {
  switch (d) {
    case D::Accept: return "accept";
    case D::FallbackLegacy: return "fallback-legacy";
    case D::RejectUpgoing: return "reject-upgoing";
    case D::RejectBeyondCrit: return "reject-beyond-critical";
    case D::RejectBelowObs: return "reject-below-obs";
    case D::RejectAboveTop: return "reject-above-top";
    case D::RejectLambdaOut: return "reject-lambda-out";
    case D::RejectFOut: return "reject-f-out";
    default: return "absorption-off";
  }
}
static D decide(double h, double b, double lam, double f, bool absorptionEnabled,
                bool tableAvailable, bool upgoing) {
  if (!absorptionEnabled) return D::AbsorptionOff;                  // no draw at all
  if (!tableAvailable) return D::FallbackLegacy;                    // enabled but no table
  if (upgoing) return D::RejectUpgoing;
  if (h < 0.0) return D::RejectBelowObs;
  if (h > 50000.0) return D::RejectAboveTop;
  if (lam < 300.0 || lam > 700.0) return D::RejectLambdaOut;
  if (f < 0.9713 || f > 1.0347) return D::RejectFOut;               // proposed range
  double const rE = R + h;
  if (b > nPhase(rE, f)*rE) return D::RejectBeyondCrit;             // cannot launch
  double const bCrit = nPhase(R, f)*R;
  if (b > bCrit) return D::RejectBeyondCrit;                        // launchable but turns
  return D::Accept;
}
// trajectory: does the ray reach r_obs? returns the minimum r reached
static double minR(double hEmit, double b, double f, double ds) {
  double const rE = R + hEmit;
  double const sE = b/(nPhase(rE,f)*rE);
  if (sE >= 1.0) return rE;
  double r = rE, th = std::acos(-std::sqrt(1.0 - sE*sE)), rmin = rE;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const Nv = N0*std::exp(-(rr-R)/HSC)*f, np = -Nv/HSC, nn = 1.0 + Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c; dt = -s*(1.0/rr + np/nn);
  };
  for (long i = 0; i < 200000000L; ++i) {
    double a1,b1,a2,b2,a3,b3,a4,b4;
    der(r,th,a1,b1); der(r+0.5*ds*a1, th+0.5*ds*b1, a2,b2);
    der(r+0.5*ds*a2, th+0.5*ds*b2, a3,b3); der(r+ds*a3, th+ds*b3, a4,b4);
    r += ds/6.0*(a1+2*a2+2*a3+a4); th += ds/6.0*(b1+2*b2+2*b3+b4);
    if (r < rmin) rmin = r;
    if (r <= R) return R;                       // reached the observation sphere
    if (rmin < r - 1.0 && r > rmin + 1.0) return rmin;   // turned and rising
  }
  return rmin;
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  double const f = 1.0, bCrit = nPhase(R,f)*R;
  // ---- (1) below observation, including the h = 0 boundary ----------------
  std::printf("(1) below observation / h=0 boundary:"); puts("");
  for (double h : {-100.0, -1.0, 0.0, 1.0, 100.0}) {
    D const d = decide(h, 0.5*bCrit, 400.0, f, true, true, false);
    std::printf("    h=%8.1f m -> %-22s", h, dn(d)); puts("");
  }
  // ---- (2) table unavailable while absorption ENABLED ---------------------
  std::printf("(2) absorption enabled, table unavailable vs absorption disabled:"); puts("");
  std::printf("    enabled + table unavailable -> %s", dn(decide(20000.0, 0.5*bCrit, 400.0, f, true, false, false))); puts("");
  std::printf("    enabled + table available   -> %s", dn(decide(20000.0, 0.5*bCrit, 400.0, f, true, true, false))); puts("");
  std::printf("    disabled (any table state)  -> %s", dn(decide(20000.0, 0.5*bCrit, 400.0, f, false, true, false))); puts("");
  // ---- (3) launchable but turns above the ground (b > b_crit) -------------
  std::printf("(3) launchable-but-turning rays (b between b_crit and n(r_e) r_e):"); puts("");
  for (double hE : {10000.0, 30000.0, 50000.0}) {
    double const rE = R + hE, bTop = nPhase(rE,f)*rE;
    for (double frac : {0.5, 0.9, 0.999}) {
      double const b = bCrit + frac*(bTop - bCrit);
      D const d = decide(hE, b, 400.0, f, true, true, false);
      double const rmin = minR(hE, b, f, 2.0);
      std::printf("    h=%6.0f b/b_crit=%.4f -> %-22s min_r - R = %9.1f m", hE, b/bCrit, dn(d), rmin - R); puts("");
    }
  }
  // ---- (4) f range endpoints and just outside -----------------------------
  std::printf("(4) f_phase range [0.9713, 1.0347] as proposed eligibility:"); puts("");
  for (double fv : {0.9709, 0.9713, 1.0000, 1.0347, 1.0351}) {
    D const d = decide(20000.0, 0.5*bCrit, 400.0, fv, true, true, false);
    std::printf("    f=%.4f -> %-22s", fv, dn(d)); puts("");
  }
  // ---- decision counts over a matrix --------------------------------------
  int nA=0, nFL=0, nU=0, nBC=0, nBO=0, nAT=0, nLO=0, nFO=0, nOff=0;
  for (double h : {-100.0, 0.0, 1000.0, 20000.0, 50000.0, 60000.0})
    for (double bb : {0.5*bCrit, 0.999*bCrit, 1.0001*bCrit, 1.05*bCrit})
      for (double lam : {290.0, 300.0, 700.0, 710.0})
        for (double fv : {0.9709, 0.9713, 1.0347, 1.0351}) {
          D const d = decide(h, bb, lam, fv, true, true, false);
          switch (d) { case D::Accept: ++nA; break; case D::FallbackLegacy: ++nFL; break;
            case D::RejectUpgoing: ++nU; break; case D::RejectBeyondCrit: ++nBC; break;
            case D::RejectBelowObs: ++nBO; break; case D::RejectAboveTop: ++nAT; break;
            case D::RejectLambdaOut: ++nLO; break; case D::RejectFOut: ++nFO; break;
            default: ++nOff; }
        }
  std::printf("matrix: accept %d | fallback %d | upgoing %d | beyond-crit %d | below-obs %d | above-top %d | lambda-out %d | f-out %d",
              nA, nFL, nU, nBC, nBO, nAT, nLO, nFO); puts("");
  return 0;
}
