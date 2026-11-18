// ============================================================================
// cp2fb_domain_b.cpp -- Checkpoint 2f-b (T2): domain contract and altitude
// coordinates.  DIAGNOSTIC ONLY.  DEVELOPMENT PROBES ONLY.
//
//   ./cp2fb_domain_b <real profile .dat>
//
// CONTRACT NOTES (ruling of 2f-b/0, corrections 1-4):
//   * Reachability uses the GENERAL definition  p_crit = min over the route of
//     n(r) r, located per cell (endpoints + the root of g' = 0).  The
//     observation-sphere value g(r_obs) is only the minimum when g is strictly
//     increasing on the route; that is verified, not assumed.
//   * An INTERIOR minimum of g is not itself a turning point.  Monotonicity is
//     reported separately from reachability; the turning point of a ray with
//     p > g_min is the radius where g(r) = p, and that is checked.
//   * Equality  p == g_min  has its own status and is NOT accepted.
//   * Physical reachability and EVALUATION-DOMAIN eligibility (wavelength) are
//     separate: an out-of-domain wavelength never makes a ray unreachable.
//   * STABLE-KERNEL SEMANTICS: extinction = native profile slopes x overlap
//     weights of the production +-25 m window (T1).  The T1 FLOATING-POINT
//     equivalence gate remains FAILED (211/10621107 queries at h = H+25-1e-6);
//     nothing here waives it or asserts source-rounding equivalence.
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R_E = 6371000.0;  // constants::EarthRadius::Mean [m]
static constexpr double H_TOP = 50000.0;  // profile top (last native knot) [m]
static constexpr double kHalf = 25.0;     // Cherenkov.inl:97, dh/2
static constexpr double kDh = 50.0;

// --- refractivity models ------------------------------------------------------
static double Nshape(double h) {
  return 2.20e-4 * std::exp(-h / 7500.0) + 0.55e-4 * std::exp(-h / 1600.0);
}
// SAME function as every earlier checkpoint; tabulated to hTopKM so the 120 km
// probes are inside the table (input representation, not an evaluation claim)
static RefractiveIndexTable makeRI(double hTopKM) {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= hTopKM * 1000.0 + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(std::log(Nshape(h)));
  }
  return RefractiveIndexTable(alt, lnN, R_E);
}
// exact vacuum: lnN = -1000 -> N = 0 at every node (RefractiveIndexTable.hpp:141)
static RefractiveIndexTable makeVacuumRI(double hTopKM) {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= hTopKM * 1000.0 + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(-1000.0);
  }
  return RefractiveIndexTable(alt, lnN, R_E);
}
// synthetic profile whose g(r) = (1 + f N) r has an INTERIOR minimum:
// N = 1e-2 exp(-h/7500) gives g'(r) < 0 near the ground and > 0 aloft.
static RefractiveIndexTable makeNonMonotoneRI(double hTopKM) {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= hTopKM * 1000.0 + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(std::log(1.0e-2) - h / 7500.0);
  }
  return RefractiveIndexTable(alt, lnN, R_E);
}

// --- profile -----------------------------------------------------------------
struct Profile {
  std::vector<double> lam, hgt;
  std::vector<std::vector<double> > tau;
};
static bool readProfile(std::string const& fn, Profile& p) {
  std::ifstream f(fn.c_str());
  if (!f) return false;
  int nw = 0, nh = 0;
  f >> nw >> nh;
  if (nw <= 0 || nh <= 0) return false;
  p.lam.resize(nw);
  p.hgt.resize(nh);
  for (int i = 0; i < nw; ++i) f >> p.lam[i];
  for (int i = 0; i < nh; ++i) f >> p.hgt[i];
  p.tau.assign(nw, std::vector<double>(nh, 0.0));
  for (int i = 0; i < nw; ++i)
    for (int j = 0; j < nh; ++j) f >> p.tau[i][j];
  return true;
}
static void slopesAt(Profile const& p, double lam, std::vector<double>& s) {
  std::size_t const nw = p.lam.size();
  double const l = std::max(p.lam.front(), std::min(p.lam.back(), lam));
  std::size_t iw = 0;
  for (std::size_t i = 0; i + 1 < nw; ++i) {
    if (p.lam[i + 1] >= l) {
      iw = i;
      break;
    }
    iw = i;
  }
  std::size_t const iw1 = std::min(iw + 1, nw - 1);
  double const fw = (iw1 > iw) ? (l - p.lam[iw]) / (p.lam[iw1] - p.lam[iw]) : 0.0;
  s.assign(p.hgt.size() - 1, 0.0);
  for (std::size_t j = 0; j + 1 < p.hgt.size(); ++j) {
    double const dh = p.hgt[j + 1] - p.hgt[j];
    double const s0 = (p.tau[iw][j + 1] - p.tau[iw][j]) / dh;
    double const s1 = (p.tau[iw1][j + 1] - p.tau[iw1][j]) / dh;
    s[j] = s0 * (1.0 - fw) + s1 * fw;
  }
}

// --- T1 STABLE kernel at PROFILE altitude z (native slopes x overlap weights) --
static double alphaKernelZ(std::vector<double> const& s, std::vector<double> const& hKnots,
                           double z) {
  double const lo = std::max(0.0, z - kHalf);
  double const hi = std::min(H_TOP, z + kHalf);
  double const d = hi - lo;
  if (!(d > 0.0)) return 0.0;
  double acc = 0.0;
  for (std::size_t j = 0; j + 1 < hKnots.size(); ++j) {
    double const a = std::max(hKnots[j], lo);
    double const b = std::min(hKnots[j + 1], hi);
    if (b > a) acc += s[j] * ((b - a) / d);
  }
  return acc;
}
// source-faithful rule (diagnostic cross-check only; T1's FP gate stays FAILED)
static double alphaSrcZ(AtmosphericAbsorption const& abs, double lam, double z) {
  double const lo = std::max(0.0, z - kHalf);
  double const hi = std::min(H_TOP, z + kHalf);
  double const d = hi - lo;
  if (!(d > 0.0)) return 0.0;
  double const t0 = abs.getOpticalDepth(lam, lo);
  double const t1 = abs.getOpticalDepth(lam, hi);
  return std::max(0.0, (t1 - t0) / d);
}

// ============================================================================
// Status contract.  Reachability and monotonicity are SEPARATE observables.
// ============================================================================
enum class Status {
  Eligible,
  BelowObservation,
  Upgoing,
  CriticalOrBeyond,
  CriticalEquality,      // p == g_min exactly: explicit, NOT accepted
  UnsupportedNonMonotone,
  InvalidInput
};
static char const* statusName(Status s) {
  switch (s) {
    case Status::Eligible: return "Eligible";
    case Status::BelowObservation: return "BelowObservation";
    case Status::Upgoing: return "Upgoing";
    case Status::CriticalOrBeyond: return "CriticalOrBeyond";
    case Status::CriticalEquality: return "CriticalEquality";
    case Status::UnsupportedNonMonotone: return "UnsupportedNonMonotone";
    case Status::InvalidInput: return "InvalidInput";
  }
  return "?";
}
enum class Reach { Reachable, Equality, Unreachable };
static char const* reachName(Reach r) {
  switch (r) {
    case Reach::Reachable: return "reachable";
    case Reach::Equality: return "equality";
    case Reach::Unreachable: return "unreachable";
  }
  return "?";
}

// exact scan of the route [rLo, rHi]: min of g' (monotonicity) and min of g
// (reachability), per cell, using the piecewise-exponential N model.
struct RouteScan {
  bool ok{false};               // cell conditioning satisfied (g' monotone in cell)
  double gPrimeMin{0.0}, gPrimeMinAt{0.0};
  double gMin{0.0}, gMinAt{0.0};
};

static RouteScan scanRoute(RefractiveIndexTable const& ri, double f, double rLo, double rHi) {
  RouteScan s;
  if (!ri.valid() || ri.size() < 2) return s;
  s.gPrimeMin = 1e300;
  s.gMin = 1e300;
  if (!(rHi > rLo + 1e-12)) {  // zero-length route
    s.ok = true;
    s.gPrimeMin = 1.0;
    s.gPrimeMinAt = rLo;
    s.gMin = (1.0 + f * ri.N(rLo)) * rLo;
    s.gMinAt = rLo;
    return s;
  }
  double const dz = (ri.zMax() - ri.zMin()) / double(ri.size() - 1);
  double const lastB = ri.logSlope(ri.size() - 2);
  bool cond = true;
  double a = rLo;
  for (int guard = 0; guard < 200000 && a < rHi; ++guard) {
    double bLog, cellHi;
    if (a >= ri.zMax()) {
      bLog = lastB;  // exponential extrapolation, same log-slope
      cellHi = rHi;
    } else {
      std::size_t const i = ri.cellIndex(a);
      bLog = ri.logSlope(i);
      cellHi = std::min(rHi, ri.zMin() + double(i + 1) * dz);
    }
    if (!(cellHi > a)) cellHi = rHi;
    auto gAt = [&](double r) { return (1.0 + f * ri.N(r)) * r; };
    auto gpAt = [&](double r) { return 1.0 + f * ri.N(r) * (1.0 + r * bLog); };
    double const ga = gAt(a), gb = gAt(cellHi);
    double const gpa = gpAt(a), gpb = gpAt(cellHi);
    // g' monotone in the cell requires g'' = f N b (2 + b r) >= 0
    if (bLog < 0.0 && a < -2.0 / bLog) cond = false;
    if (gpa < s.gPrimeMin) {
      s.gPrimeMin = gpa;
      s.gPrimeMinAt = a;
    }
    if (gpb < s.gPrimeMin) {
      s.gPrimeMin = gpb;
      s.gPrimeMinAt = cellHi;
    }
    double candR;
    if (gpa >= 0.0) {
      candR = a;                      // g increasing on the cell
    } else if (gpb <= 0.0) {
      candR = cellHi;                 // g decreasing
    } else {
      double lo = a, hi = cellHi;     // one interior root of g' (g' increasing)
      for (int it = 0; it < 200; ++it) {
        double const mid = 0.5 * (lo + hi);
        if (gpAt(mid) > 0.0) {
          hi = mid;
        } else {
          lo = mid;
        }
      }
      candR = 0.5 * (lo + hi);
    }
    if (ga < s.gMin) {
      s.gMin = ga;
      s.gMinAt = a;
    }
    if (gb < s.gMin) {
      s.gMin = gb;
      s.gMinAt = cellHi;
    }
    double const gc = gAt(candR);
    if (gc < s.gMin) {
      s.gMin = gc;
      s.gMinAt = candR;
    }
    if (!(cellHi > a)) break;
    a = cellHi;
  }
  s.ok = cond && std::isfinite(s.gMin) && std::isfinite(s.gPrimeMin);
  return s;
}

struct Query {
  Status status{Status::InvalidInput};
  Reach reach{Reach::Unreachable};
  double rObs{0.0}, rE{0.0}, zObs{0.0}, zE{0.0};
  double p{0.0}, gObs{0.0}, gMin{0.0}, rAtMin{0.0}, muObs{0.0}, f{0.0};
  double gPrimeMin{0.0}, gPrimeMinAt{0.0};
  bool routeOk{false}, monotone{false}, launchable{false};
  char const* note{""};
};


// ============================================================================
// Classification.  Order: InvalidInput -> BelowObservation -> Upgoing ->
// reachability (Equality / Beyond) -> monotonicity -> Eligible.
// ============================================================================
static Query classifyByP(RefractiveIndexTable const& ri, double zObs, double hEmit, double p,
                         double f) {
  Query q;
  q.zObs = zObs;
  q.rObs = R_E + zObs;
  q.zE = zObs + hEmit;
  q.rE = q.rObs + hEmit;
  q.f = f;
  q.p = p;
  if (!std::isfinite(p) || !(p >= 0.0) || !std::isfinite(hEmit) || !std::isfinite(zObs) ||
      !(f > 0.0) || !std::isfinite(f) || !(q.rE > 0.0)) {
    q.status = Status::InvalidInput;
    q.note = "non-finite or non-positive input";
    return q;
  }
  if (hEmit < 0.0) {
    q.status = Status::BelowObservation;
    q.note = "emission below the observation radius";
    return q;
  }
  RouteScan const s = scanRoute(ri, f, q.rObs, q.rE);
  q.routeOk = s.ok;
  q.gMin = s.gMin;
  q.rAtMin = s.gMinAt;
  q.gPrimeMin = s.gPrimeMin;
  q.gPrimeMinAt = s.gPrimeMinAt;
  q.monotone = s.ok && (s.gPrimeMin > 0.0);
  q.gObs = (1.0 + f * ri.N(q.rObs)) * q.rObs;
  double const gE = (1.0 + f * ri.N(q.rE)) * q.rE;
  q.launchable = (p <= gE);  // PROTOCOL 2: independent launchability check
  if (!(s.gMin > 0.0)) {
    q.status = Status::InvalidInput;
    q.note = "route minimum not established";
    return q;
  }
  double const margin = p / s.gMin - 1.0;
  if (margin > 0.0) {
    q.reach = Reach::Unreachable;
    q.status = Status::CriticalOrBeyond;
    q.note = "p above the route minimum: turns before the observation level";
    return q;
  }
  if (margin == 0.0) {
    q.reach = Reach::Equality;
    q.status = Status::CriticalEquality;  // explicit, NOT accepted
    q.note = "p equals the route minimum exactly";
    return q;
  }
  q.reach = Reach::Reachable;
  if (!q.monotone) {
    q.status = Status::UnsupportedNonMonotone;
    q.note = "reachable but g is not strictly increasing on the route";
    return q;
  }
  q.status = Status::Eligible;
  return q;
}

static Query classifyMu(RefractiveIndexTable const& ri, double zObs, double hEmit, double muObs,
                        double f) {
  if (!std::isfinite(muObs) || muObs < 0.0 || muObs > 1.0 || !(f > 0.0) || !std::isfinite(f)) {
    Query q;
    q.status = Status::InvalidInput;
    q.note = "mu_obs outside [0,1] or f <= 0";
    return q;
  }
  double const rObs = R_E + zObs;
  double const gObs = (1.0 + f * ri.N(rObs)) * rObs;
  double const p = gObs * std::sqrt(std::max(0.0, 1.0 - muObs * muObs));
  Query q = classifyByP(ri, zObs, hEmit, p, f);
  q.muObs = muObs;
  return q;
}

// emission-direction form: cEOut = dot(direction, outward radial) at emission
static Query classifyEmission(RefractiveIndexTable const& ri, double zObs, double hEmit,
                              double cEOut, double f) {
  Query q;
  if (!std::isfinite(cEOut) || std::fabs(cEOut) > 1.0 || !(f > 0.0) || !(hEmit >= 0.0)) {
    q.status = Status::InvalidInput;
    q.note = "non-finite direction or h < 0 or f <= 0";
    return q;
  }
  double const rE = R_E + zObs + hEmit;
  if (cEOut >= 0.0) {
    q.status = Status::Upgoing;
    q.zObs = zObs;
    q.rObs = R_E + zObs;
    q.rE = rE;
    q.f = f;
    q.note = "direction has a non-negative outward radial component";
    return q;
  }
  double const muE = -cEOut;
  double const gE = (1.0 + f * ri.N(rE)) * rE;
  double const p = gE * std::sqrt(std::max(0.0, 1.0 - muE * muE));
  q = classifyByP(ri, zObs, hEmit, p, f);
  q.launchable = true;
  q.muObs = std::sqrt(std::max(0.0, 1.0 - (p / q.gObs) * (p / q.gObs)));
  return q;
}

// radius where g(r) = p on [rA, rB] (NaN when not bracketed) -- used to show
// that a beyond-critical ray turns ABOVE the observation level
static double findTurning(RefractiveIndexTable const& ri, double f, double p, double rA,
                          double rB) {
  auto gAt = [&](double r) { return (1.0 + f * ri.N(r)) * r; };
  double fa = gAt(rA) - p;
  double fb = gAt(rB) - p;
  if ((fa > 0.0) == (fb > 0.0)) return std::nan("");
  for (int it = 0; it < 200; ++it) {
    double const m = 0.5 * (rA + rB);
    double const fm = gAt(m) - p;
    if ((fm > 0.0) == (fa > 0.0)) {
      rA = m;
      fa = fm;
    } else {
      rB = m;
    }
  }
  return 0.5 * (rA + rB);
}


// ============================================================================
// Provisional, UNCERTIFIED optical-depth integrator (section 5 only).
//   tau = int_{rObs}^{rE} alpha(z(r)) dr / cos(theta(r))
//   cos(theta(r)) = sqrt(1 - (p/g(r))^2),  g = (1 + f N) r,  z = r - R_E
// Panels split at every kernel boundary; nSub sub-panels of 8-point GL each.
// ============================================================================
static double const kXG[8] = {-0.9602898564975363, -0.7966664774136267,
                              -0.5255324099163290, -0.1834346424956498,
                              0.1834346424956498,  0.5255324099163290,
                              0.7966664774136267,  0.9602898564975363};
static double const kWG[8] = {0.1012285362903763, 0.2223810344533745,
                              0.3137066458778873, 0.3626837833783620,
                              0.3626837833783620, 0.3137066458778873,
                              0.2223810344533745, 0.1012285362903763};

static double tauPanel(std::vector<double> const& s, Profile const& prof,
                       RefractiveIndexTable const& ri, double f, double p, double rA,
                       double rB, int nSub) {
  double acc = 0.0;
  double const h = (rB - rA) / double(nSub);
  for (int k = 0; k < nSub; ++k) {
    double const a = rA + double(k) * h;
    double const b = a + h;
    double const c = 0.5 * (a + b);
    double const hh = 0.5 * (b - a);
    for (int i = 0; i < 8; ++i) {
      double const r = c + hh * kXG[i];
      double const g = (1.0 + f * ri.N(r)) * r;
      double const q = p / g;
      double const cosT = std::sqrt(std::max(1e-300, 1.0 - q * q));
      double const al = alphaKernelZ(s, prof.hgt, r - R_E);
      acc += kWG[i] * hh * al / cosT;
    }
  }
  return acc;
}

static double tauProvisional(std::vector<double> const& s, Profile const& prof,
                             RefractiveIndexTable const& ri, double f, double p,
                             double rObs, double rE, int nSub, int* nPanels = nullptr) {
  std::vector<double> zb;
  zb.push_back(0.0);
  for (std::size_t j = 0; j < prof.hgt.size(); ++j) {
    zb.push_back(prof.hgt[j]);
    zb.push_back(prof.hgt[j] - kHalf);
    zb.push_back(prof.hgt[j] + kHalf);
  }
  zb.push_back(H_TOP);
  zb.push_back(H_TOP + kHalf);
  double const zLo = rObs - R_E;
  double const zE = rE - R_E;
  zb.push_back(zE);
  std::sort(zb.begin(), zb.end());
  zb.erase(std::unique(zb.begin(), zb.end()), zb.end());
  double acc = 0.0;
  int np = 0;
  for (std::size_t k = 0; k + 1 < zb.size(); ++k) {
    double const a = R_E + std::max(zLo, zb[k]);
    double const b = R_E + std::min(zE, zb[k + 1]);
    if (!(b > a)) continue;
    acc += tauPanel(s, prof, ri, f, p, a, b, nSub);
    ++np;
  }
  if (nPanels) *nPanels = np;
  return acc;
}

// vacuum distance integral, graded panels (section 3 cross-check)
static double intLquad(double p, double rA, double rB) {
  int const n = 400;
  double acc = 0.0;
  for (int k = 0; k < n; ++k) {
    double const t0 = double(k) / double(n), t1 = double(k + 1) / double(n);
    double const a = rA + (rB - rA) * t0 * t0;
    double const b = rA + (rB - rA) * t1 * t1;
    double const c = 0.5 * (a + b), hh = 0.5 * (b - a);
    for (int i = 0; i < 8; ++i) {
      double const r = c + hh * kXG[i];
      double const w = std::sqrt(std::max(1e-300, (r - p) * (r + p)));
      acc += kWG[i] * hh * r / w;
    }
  }
  return acc;
}


// ============================================================================
int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <real profile .dat>\n", argv[0]);
    return 2;
  }
  Profile prof;
  if (!readProfile(argv[1], prof)) {
    std::printf("STOP: cannot read the profile\n");
    return 1;
  }
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[1]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  RefractiveIndexTable const ri = makeRI(120.0);
  RefractiveIndexTable const vac = makeVacuumRI(120.0);
  RefractiveIndexTable const bad = makeNonMonotoneRI(120.0);

  std::printf("=== Checkpoint 2f-b (T2): domain contract and altitude coordinates ===\n");
  std::printf("profile : %s\n", argv[1]);
  std::printf("SEMANTICS: stable production kernel (T1). The T1 FLOATING-POINT equivalence gate\n");
  std::printf("  remains FAILED (211/10621107 queries at h = H+25-1e-6); this run does not waive\n");
  std::printf("  it and does not assert source-rounding equivalence.\n");
  std::printf("PROBES: development only; no frozen held-out set is used.\n\n");

  int nGate = 0, nFail = 0;
  auto gate = [&](char const* what, double residual, double tol) {
    ++nGate;
    bool const ok = (residual <= tol);
    if (!ok) ++nFail;
    std::printf("  %-52s residual=%.6e  tol=%.1e  %s\n", what, residual, tol,
                ok ? "PASS" : "FAIL");
  };

  // ======================= 1. coordinates and index semantics ==============
  std::printf("--- 1. coordinates and refractive-index semantics ---\n");
  std::printf("  R_E = %.17g m\n", R_E);
  std::printf("  z = r - R_E (profile altitude) ; h = r - r_obs ; z_obs = r_obs - R_E\n");
  std::printf("  kernel altitude  z = h + z_obs   (NOT h when z_obs != 0)\n");
  std::printf("  RI table: [%.17g, %.17g] m, %zu nodes, dz = %.17g m\n", ri.zMin(), ri.zMax(),
              ri.size(), (ri.zMax() - ri.zMin()) / double(ri.size() - 1));
  std::printf("  %-22s %-20s %-20s %-20s\n", "z [m]", "N(z)", "dN/dz [1/m]", "logSlope [1/m]");
  double const zs[4] = {0.0, 1000.0, 20000.0, 119000.0};
  for (int i = 0; i < 4; ++i) {
    double const z = zs[i];
    std::printf("  %-22.17g %-20.17g %-20.17g %-20.17g\n", z, ri.N(R_E + z), ri.dNdz(R_E + z),
                ri.logSlope(ri.cellIndex(R_E + z)));
  }
  // above-range rule: N = N(zMax) * exp(b_last (z - zMax)), floor at 1e-10 -> 0
  double const zTop = ri.zMax();
  double const bLast = ri.logSlope(ri.size() - 2);
  double const nTop = ri.N(zTop);
  double const nAbove = ri.N(zTop + 1000.0);
  double const nRaw = nTop * std::exp(bLast * 1000.0);
  std::printf("  above-range rule: N(zMax) = %.17g ; 1e-10 floor ; N(zMax+1000) = %.17g ;\n",
              nTop, nAbove);
  std::printf("                    unfloored exponential would give %.17g\n", nRaw);
  gate("1 above-range: exactly 0 once below the 1e-10 floor",
       (nTop < 1.0e-10) ? std::fabs(nAbove) : 0.0, 0.0);
  std::printf("\n");

  // ======================= 3. vacuum sanity checks =========================
  // REVISED (checkpoint 2f-b2): the invariant cosine is evaluated in FACTORED
  // form sqrt((r-p)(r+p))/r instead of sqrt(1-(p/r)^2).  Only that formula
  // changed; probes and acceptance tolerances are unchanged.  The unfactored
  // value is kept as an INFORMATIONAL historical comparison, not a gate.
  std::printf("--- 3. vacuum (n=1) checks [REVISED: FACTORED INVARIANT] ---\n");
  {
    using mp = boost::multiprecision::cpp_dec_float_100;
    double const zObs = 0.0, hE = 1.0, muObs = 1.0e-6;  // previous worst row
    double const rObs = R_E + zObs, rE = rObs + hE;
    double const p = rObs * std::sqrt(std::max(0.0, 1.0 - muObs * muObs));
    double const pr = p / rE;
    double const muE = std::sqrt(std::max(0.0, 1.0 - pr * pr));
    double const sObs = std::sqrt(std::max(0.0, (rObs - p) * (rObs + p)));
    double const L = hE * (rE + rObs) / (rE * muE + sObs);
    double const rs = std::hypot(L - rE * muE, p);
    double const line = (rE * muE - L) / rs;
    double const fact = std::sqrt((rs - p) * (rs + p)) / rs;
    mp const P(p), RR(rs);
    mp const ref = sqrt((RR - P) * (RR + P)) / RR;
    double const dFactRef = std::fabs((mp(fact) - ref).convert_to<double>());
    double const dFactLine = std::fabs(fact - line);
    std::printf("  pre-flight, previous worst row (z_obs=0, h=1 m, mu_obs=1e-6):\n");
    std::printf("    mu_fact = %.17g ; mu_line = %.17g ; 100-digit stored-input ref = %.25e\n",
                fact, line, ref.convert_to<double>());
    std::printf("    |mu_fact - ref| = %.3e ; |mu_fact - mu_line| = %.3e\n", dFactRef, dFactLine);
    gate("3 pre-flight |mu_fact - 100-digit ref|", dFactRef, 1.0e-12);
    gate("3 pre-flight |mu_fact - mu_line|", dFactLine, 1.0e-12);
  }
  std::printf("  %-8s %-9s %-10s %-16s %-13s %-14s %-12s %-12s %-12s\n", "z_obs", "h", "mu_obs",
              "p [m]", "mu_e", "L [m]", "|r(L)-r_obs|", "resid_rev", "resid_hist");
  double const zobsV[2] = {0.0, 2000.0};
  double const hV[4] = {1.0, 1000.0, 50000.0, 120000.0};
  double const muV[5] = {1.0, 0.5, 0.1, 2.0e-4, 1.0e-6};
  double worstMu = 0.0, worstR = 0.0, worstL = 0.0, worstLq = 0.0;
  double worstMuWell = 0.0, worstRatio = 0.0, worstMuOld = 0.0;
  long nInvalid = 0;
  for (int iz = 0; iz < 2; ++iz)
    for (int ih = 0; ih < 4; ++ih)
      for (int im = 0; im < 5; ++im) {
        double const zObs = zobsV[iz], hE = hV[ih], muObs = muV[im];
        double const rObs = R_E + zObs, rE = rObs + hE;
        double const p = rObs * std::sqrt(std::max(0.0, 1.0 - muObs * muObs));
        double const pr = p / rE;
        double const muE = std::sqrt(std::max(0.0, 1.0 - pr * pr));
        double const sObs = std::sqrt(std::max(0.0, (rObs - p) * (rObs + p)));
        double const sE = std::sqrt(std::max(0.0, (rE - p) * (rE + p)));
        double const L = hE * (rE + rObs) / (rE * muE + sObs);  // rationalised
        double const Lalt = sE - sObs;
        double dmuNew = 0.0, dmuOld = 0.0;
        for (int k = 0; k <= 4; ++k) {
          double const s = L * double(k) / 4.0;
          // stable: r(s)^2 = (s - rE muE)^2 + p^2  (no cancellation)
          double const rs = std::hypot(s - rE * muE, p);
          double const muLine = (rE * muE - s) / rs;
          // domain guard: r > 0 and 0 <= p <= r ; NO unconditional clamping
          if (!(rs > 0.0) || !(p >= 0.0) || !(p <= rs)) {
            ++nInvalid;
            continue;
          }
          double const muFact = std::sqrt((rs - p) * (rs + p)) / rs;  // REVISED
          double const q = p / rs;
          double const muOld = std::sqrt(std::max(0.0, 1.0 - q * q));  // historical
          dmuNew = std::max(dmuNew, std::fabs(muLine - muFact));
          dmuOld = std::max(dmuOld, std::fabs(muLine - muOld));
        }
        double const rAtL = std::hypot(L - rE * muE, p);
        double const dR = std::fabs(rAtL - rObs);
        bool const quadOk = (muObs >= 1.0e-2);
        double const lq = quadOk ? std::fabs(intLquad(p, rObs, rE) - L) : 0.0;
        worstMu = std::max(worstMu, dmuNew);
        worstMuOld = std::max(worstMuOld, dmuOld);
        if (muObs >= 1.0e-3) worstMuWell = std::max(worstMuWell, dmuNew);
        // sensitivity scale (NOT a bound): recovering cos(theta_obs) from p near
        // r_obs has ABSOLUTE sensitivity ~ ulp(p)/(g_obs mu_obs) = O(eps/mu_obs)
        {
          double const scaleAbs = (muObs > 0.0) ? (double(2.220446049250313e-16) / muObs) : 0.0;
          if (scaleAbs > 0.0) worstRatio = std::max(worstRatio, dmuNew / scaleAbs);
        }
        worstR = std::max(worstR, dR);
        worstL = std::max(worstL, std::fabs(Lalt - L));
        if (quadOk) worstLq = std::max(worstLq, lq);
        std::printf("  %-8.17g %-9.17g %-10.17g %-16.17g %-13.17g %-14.17g %-12.3e %-12.3e "
                    "%-12.3e\n",
                    zObs, hE, muObs, p, muE, L, dR, dmuNew, dmuOld);
      }
  std::printf("\n");
  gate("3 max |mu_line - mu_fact| REVISED (all cases)", worstMu, 1.0e-12);
  gate("3 max |mu_line - mu_fact| REVISED (mu_obs >= 1e-3)", worstMuWell, 1.0e-12);
  gate("3 cosine-domain validity failures (r>0, 0<=p<=r)", double(nInvalid), 0.0);
  gate("3 residual <= eps/mu_obs sensitivity scale (ratio)", worstRatio, 1.0);
  gate("3 max |r(L) - r_obs| [m]", worstR, 1.0e-6);
  gate("3 max |L_alt - L| [m]", worstL, 1.0e-6);
  gate("3 max |L_quad - L| (mu_obs >= 0.01) [m]", worstLq, 1.0e-6);
  std::printf("  informational, NOT a gate: historical unfactored max |mu_line - "
              "sqrt(1-(p/r)^2)| = %.3e\n",
              worstMuOld);
  std::printf("\n");


  // ======================= 4a. development grid ============================
  const double kZobs[2] = {0.0, 2000.0};
  const double kH[7] = {0.0, 1.0, 100.0, 1000.0, 1.0e4, 5.0e4, 1.2e5};
  const double kMu[6] = {2.0e-4, 1.0e-3, 0.01, 0.1, 0.5, 1.0};
  const double kF[3] = {0.96, 1.0, 1.04};
  std::printf("--- 4a. development grid: 2 z_obs x 7 h x 6 mu_obs x 3 f = 252 ---\n");
  long counts[7] = {0, 0, 0, 0, 0, 0, 0};
  int nTot = 0;
  double rgMin[2][3], rgAt[2][3], rgpMin[2][3];
  for (int iz = 0; iz < 2; ++iz)
    for (int iff = 0; iff < 3; ++iff) {
      rgMin[iz][iff] = 1e300;
      rgAt[iz][iff] = 0.0;
      rgpMin[iz][iff] = 1e300;
    }
  for (int iz = 0; iz < 2; ++iz)
    for (int iff = 0; iff < 3; ++iff)
      for (int ih = 0; ih < 7; ++ih)
        for (int im = 0; im < 6; ++im) {
          Query const q = classifyMu(ri, kZobs[iz], kH[ih], kMu[im], kF[iff]);
          ++counts[int(q.status)];
          ++nTot;
          if (q.gMin > 0.0 && q.gMin < rgMin[iz][iff]) {
            rgMin[iz][iff] = q.gMin;
            rgAt[iz][iff] = q.rAtMin;
          }
          if (q.gPrimeMin < rgpMin[iz][iff]) rgpMin[iz][iff] = q.gPrimeMin;
          if (q.status != Status::Eligible) {
            std::printf("    non-Eligible: z_obs=%.17g h=%.17g mu=%.17g f=%.17g -> %s "
                        "(reach=%s, p/g_min-1=%.3e)\n",
                        kZobs[iz], kH[ih], kMu[im], kF[iff], statusName(q.status),
                        reachName(q.reach), q.p / q.gMin - 1.0);
          }
        }
  std::printf("  counts:");
  for (int k = 0; k < 7; ++k) std::printf(" %s=%ld", statusName(Status(k)), counts[k]);
  std::printf("  (total %d)\n\n", nTot);

  std::printf("  route minimum of g(r)=n(r)r per (z_obs, f), over the 42 routes each:\n");
  std::printf("  %-8s %-8s %-22s %-22s %-10s %-20s\n", "z_obs", "f", "min_route g [m]",
              "r_at_min [m]", "interior?", "min g' [1/m]");
  for (int iz = 0; iz < 2; ++iz)
    for (int iff = 0; iff < 3; ++iff) {
      bool const interior = (rgAt[iz][iff] > R_E + kZobs[iz] + 1e-6);
      std::printf("  %-8.17g %-8.17g %-22.17g %-22.17g %-10s %-20.17g\n", kZobs[iz], kF[iff],
                  rgMin[iz][iff], rgAt[iz][iff], interior ? "YES" : "no", rgpMin[iz][iff]);
    }
  gate("4a all 252 grid probes classified", double(252 - nTot), 0.0);
  gate("4a grid: no InvalidInput", double(counts[int(Status::InvalidInput)]), 0.0);
  gate("4a grid: no CriticalEquality case arises", double(counts[int(Status::CriticalEquality)]), 0.0);
  std::printf("\n");

  // ======================= 4b. boundary sequences ==========================
  std::printf("--- 4b. boundary sequences (z_obs = 0, mu_obs = 0.5, f = 1) ---\n");
  std::printf("  ground sequence:\n");
  std::printf("  %-14s %-22s %-12s %-22s %-14s\n", "h [m]", "status", "reach", "p [m]", "p/g_min-1");
  double const hG[12] = {0.0,   1e-6, 1e-3,  1e-2,   0.1,    1.0,
                         10.0,  100.0, 1000.0, -1e-3, -1.0,  -100.0};
  int badNeg = 0;
  for (int i = 0; i < 12; ++i) {
    Query const q = classifyMu(ri, 0.0, hG[i], 0.5, 1.0);
    std::printf("  %-14.17g %-22s %-12s %-22.17g %-14s\n", hG[i], statusName(q.status),
                reachName(q.reach), q.p,
                (q.gMin > 0.0) ? std::to_string(q.p / q.gMin - 1.0).c_str() : "n/a");
    if (hG[i] < 0.0 && q.status != Status::BelowObservation) ++badNeg;
  }
  gate("4b negative heights -> BelowObservation", double(badNeg), 0.0);

  std::printf("  table-ceiling sequence (z_e = H + delta, h = z_e):\n");
  std::printf("  %-14s %-22s %-12s %-22s\n", "delta [m]", "status", "reach", "z_e [m]");
  double const dCeil[8] = {-100.0, -1.0, -1e-3, -1e-2, 0.0, 1e-3, 1.0, 100.0};
  int notElig = 0;
  for (int i = 0; i < 8; ++i) {
    double const zE = H_TOP + dCeil[i];
    Query const q = classifyMu(ri, 0.0, zE, 0.5, 1.0);
    std::printf("  %-14.17g %-22s %-12s %-22.17g\n", dCeil[i], statusName(q.status),
                reachName(q.reach), q.zE);
    if (q.status != Status::Eligible) ++notElig;
  }
  gate("4b ceiling sequence all Eligible", double(notElig), 0.0);
  std::printf("\n");


  // ======================= 4c. equality and beyond-critical ================
  std::printf("--- 4c. equality and beyond-critical: p against the ROUTE minimum ---\n");
  Query const ref = classifyMu(ri, 0.0, 1000.0, 0.5, 1.0);
  std::printf("  real profile: route min g = %.17g m at r = %.17g m (z = %.17g m), "
              "g(r_obs) = %.17g m\n",
              ref.gMin, ref.rAtMin, ref.rAtMin - R_E, ref.gObs);
  std::printf("  %-12s %-22s %-8s %-22s %-12s %-20s %-12s\n", "p/g_min", "p [m]", "launch",
              "status", "reach", "r_turn [m]", "cos(turn)");
  double const ratios[5] = {1.0, 1.0000005, 1.00001, 1.00005, 1.0001};
  int eqReal = 0, beyondReal = 0;
  for (int i = 0; i < 5; ++i) {
    double const p = ratios[i] * ref.gMin;
    Query const q = classifyByP(ri, 0.0, 1000.0, p, 1.0);
    double rTurn = std::nan(""), cosTurn = std::nan("");
    if (q.reach == Reach::Unreachable && q.launchable) {
      // the ray arrives from r_e: the turning point is the OUTER crossing g(r) = p,
      // i.e. on [r_at_min, r_e], not on the inner branch
      rTurn = findTurning(ri, 1.0, p, q.rAtMin, q.rE);
      if (std::isfinite(rTurn)) {
        double const g = (1.0 + ri.N(rTurn)) * rTurn;
        cosTurn = std::sqrt(std::max(0.0, 1.0 - (p / g) * (p / g)));
      }
    }
    if (q.status == Status::CriticalEquality) ++eqReal;
    if (q.status == Status::CriticalOrBeyond) ++beyondReal;
    std::printf("  %-12.17g %-22.17g %-8s %-22s %-12s %-20.17g %-12.3e\n", ratios[i], p,
                q.launchable ? "yes" : "NO", statusName(q.status), reachName(q.reach), rTurn,
                cosTurn);
  }
  gate("4c real profile: ratio 1 -> CriticalEquality", double(1 - eqReal), 0.0);
  gate("4c real profile: ratios >1 -> CriticalOrBeyond", double(4 - beyondReal), 0.0);

  std::printf("  synthetic non-monotone profile (N = 1e-2 exp(-h/7500)), z_obs=0, h=120000, "
              "f=1:\n");
  Query const badRef = classifyMu(bad, 0.0, 120000.0, 0.5, 1.0);
  std::printf("    route min g = %.17g m at r = %.17g m (z = %.17g m) : %s ; g(r_obs) = "
              "%.17g m\n",
              badRef.gMin, badRef.rAtMin, badRef.rAtMin - R_E,
              (badRef.rAtMin > R_E + 1e-6) ? "INTERIOR minimum" : "at r_obs", badRef.gObs);
  std::printf("    a reachable ray on this route: mu_obs=0.5 -> %s (reach=%s)\n",
              statusName(badRef.status), reachName(badRef.reach));
  std::printf("  %-12s %-22s %-8s %-22s %-12s %-20s %-12s\n", "p/g_min", "p [m]", "launch",
              "status", "reach", "r_turn [m]", "cos(turn)");
  double const ratiosB[5] = {1.0, 1.0005, 1.002, 1.005, 1.01};
  int eqB = 0, beyondB = 0;
  for (int i = 0; i < 5; ++i) {
    double const p = ratiosB[i] * badRef.gMin;
    Query const q = classifyByP(bad, 0.0, 120000.0, p, 1.0);
    double rTurn = std::nan(""), cosTurn = std::nan("");
    if (q.reach == Reach::Unreachable && q.launchable) {
      rTurn = findTurning(bad, 1.0, p, q.rAtMin, q.rE);
      if (std::isfinite(rTurn)) {
        double const g = (1.0 + bad.N(rTurn)) * rTurn;
        cosTurn = std::sqrt(std::max(0.0, 1.0 - (p / g) * (p / g)));
      }
    }
    if (q.status == Status::CriticalEquality) ++eqB;
    if (q.status == Status::CriticalOrBeyond) ++beyondB;
    std::printf("  %-12.17g %-22.17g %-8s %-22s %-12s %-20.17g %-12.3e\n", ratiosB[i], p,
                q.launchable ? "yes" : "NO", statusName(q.status), reachName(q.reach), rTurn,
                cosTurn);
  }
  gate("4c synthetic: ratio 1 -> CriticalEquality", double(1 - eqB), 0.0);
  gate("4c synthetic: launchable ratios >1 -> CriticalOrBeyond", double(4 - beyondB), 0.0);
  std::printf("  note: an interior minimum of g is NOT a turning point; the turning point of a\n");
  std::printf("        ray with p > g_min is where g(r) = p, and cos(turn) above is ~0 there.\n\n");


  // ======================= 4d. lambda: evaluation domain only =============
  std::printf("--- 4d. evaluation-domain eligibility (wavelength) vs physical reachability ---\n");
  std::printf("  %-12s %-12s %-12s %-22s %-12s\n", "lambda [nm]", "in [300,700]?",
              "clamped [nm]", "physical status", "reach");
  double const lamSeq[8] = {295.0, 299.999, 300.0, 300.001, 699.999, 700.0, 700.001, 705.0};
  int lamPhysDiffer = 0;
  Status const phys0 = classifyMu(ri, 0.0, 1000.0, 0.5, 1.0).status;
  for (int i = 0; i < 8; ++i) {
    double const lam = lamSeq[i];
    bool const inDom = (lam >= 300.0 && lam <= 700.0);
    double const cl = std::max(prof.lam.front(), std::min(prof.lam.back(), lam));
    Query const q = classifyMu(ri, 0.0, 1000.0, 0.5, 1.0);
    if (q.status != phys0) ++lamPhysDiffer;
    std::printf("  %-12.17g %-12s %-12.17g %-22s %-12s\n", lam, inDom ? "yes" : "OUT", cl,
                statusName(q.status), reachName(q.reach));
  }
  gate("4d physical status independent of lambda", double(lamPhysDiffer), 0.0);

  // ======================= 4e. synthetic status cases ======================
  std::printf("\n--- 4e. one explicit input per status (synthetic where needed) ---\n");
  struct Case {
    char const* name;
    Query q;
    Status expect;
  };
  Case cases[7];
  cases[0].name = "h = -1 m (below observation)";
  cases[0].q = classifyMu(ri, 0.0, -1.0, 0.5, 1.0);
  cases[0].expect = Status::BelowObservation;
  cases[1].name = "emission direction outward (c_e = +0.5)";
  cases[1].q = classifyEmission(ri, 0.0, 1000.0, 0.5, 1.0);
  cases[1].expect = Status::Upgoing;
  cases[2].name = "p = g_min exactly (equality)";
  {
    Query const r0 = classifyMu(ri, 0.0, 1000.0, 0.5, 1.0);
    cases[2].q = classifyByP(ri, 0.0, 1000.0, r0.gMin, 1.0);
  }
  cases[2].expect = Status::CriticalEquality;
  cases[3].name = "p = 1.00001 g_min (beyond critical)";
  {
    Query const r0 = classifyMu(ri, 0.0, 1000.0, 0.5, 1.0);
    cases[3].q = classifyByP(ri, 0.0, 1000.0, 1.00001 * r0.gMin, 1.0);
  }
  cases[3].expect = Status::CriticalOrBeyond;
  cases[4].name = "non-monotone synthetic profile, reachable ray";
  cases[4].q = classifyMu(bad, 0.0, 120000.0, 0.5, 1.0);
  cases[4].expect = Status::UnsupportedNonMonotone;
  cases[5].name = "mu_obs = 1.5 (invalid input)";
  cases[5].q = classifyMu(ri, 0.0, 1000.0, 1.5, 1.0);
  cases[5].expect = Status::InvalidInput;
  cases[6].name = "h = 0 (zero-length path)";
  cases[6].q = classifyMu(ri, 0.0, 0.0, 0.5, 1.0);
  cases[6].expect = Status::Eligible;
  std::vector<double> s480;
  slopesAt(prof, 480.0, s480);
  int mismatch = 0, tauLeak = 0;
  for (int i = 0; i < 7; ++i) {
    bool const ok = (cases[i].q.status == cases[i].expect);
    if (!ok) ++mismatch;
    double tau = std::nan("");
    if (cases[i].q.status == Status::Eligible) {
      tau = tauProvisional(s480, prof, ri, cases[i].q.f, cases[i].q.p, cases[i].q.rObs,
                           cases[i].q.rE, 8);
    }
    if (cases[i].q.status != Status::Eligible && std::isfinite(tau)) ++tauLeak;
    std::printf("  %-52s %-22s expect %-22s %s", cases[i].name, statusName(cases[i].q.status),
                statusName(cases[i].expect), ok ? "PASS" : "FAIL");
    if (cases[i].q.status == Status::Eligible)
      std::printf("   tau = %.17g", tau);
    std::printf("\n");
  }
  gate("4e synthetic status cases match", double(mismatch), 0.0);
  gate("4e rejected inputs produce no tau", double(tauLeak), 0.0);
  gate("4e h = 0 gives tau == 0 exactly",
       std::fabs(tauProvisional(s480, prof, ri, 1.0, cases[6].q.p, cases[6].q.rObs,
                                cases[6].q.rE, 8)),
       0.0);
  std::printf("\n");


  // ======================= 5. saturation (approved) ========================
  std::printf("--- 5. above-profile saturation: p, f, lambda and z_obs held fixed ---\n");
  std::printf("  gate: |tau(z_e) - tau(H+25)| <= 1e-10 (absolute)\n");
  double const lam5[4] = {180.0, 255.0, 480.0, 700.0};
  double const ze5[4] = {H_TOP + 25.0, H_TOP + 100.0, H_TOP + 1000.0, 120000.0};
  long nSat = 0, nSatSkip = 0;
  double worstSat = 0.0, wLam = 0.0, wZe = 0.0, wMu = 0.0, wF = 0.0, wZobs = 0.0;
  for (int iz = 0; iz < 2; ++iz)
    for (int iff = 0; iff < 3; ++iff)
      for (int im = 0; im < 6; ++im)
        for (int il = 0; il < 4; ++il) {
          std::vector<double> s;
          slopesAt(prof, lam5[il], s);
          Query const q0 = classifyMu(ri, kZobs[iz], ze5[0] - kZobs[iz], kMu[im], kF[iff]);
          if (q0.status != Status::Eligible) {
            ++nSatSkip;
            continue;
          }
          double const tau0 = tauProvisional(s, prof, ri, q0.f, q0.p, q0.rObs, q0.rE, 8);
          for (int ie = 1; ie < 4; ++ie) {
            Query const qe = classifyMu(ri, kZobs[iz], ze5[ie] - kZobs[iz], kMu[im], kF[iff]);
            if (qe.status != Status::Eligible) {
              ++nSatSkip;
              continue;
            }
            double const tauE = tauProvisional(s, prof, ri, qe.f, qe.p, qe.rObs, qe.rE, 8);
            double const res = std::fabs(tauE - tau0);
            ++nSat;
            if (res > worstSat) {
              worstSat = res;
              wLam = lam5[il];
              wZe = ze5[ie];
              wMu = kMu[im];
              wF = kF[iff];
              wZobs = kZobs[iz];
            }
          }
        }
  std::printf("  cases=%ld skipped=%ld\n", nSat, nSatSkip);
  std::printf("  worst |tau(z_e)-tau(H+25)| = %.17g", worstSat);
  if (worstSat > 0.0)
    std::printf(" at z_obs=%.17g f=%.17g mu=%.17g lambda=%.17g z_e=%.17g", wZobs, wF, wMu,
                wLam, wZe);
  else
    std::printf(" (no case exceeded zero; the argmax fields are unset)");
  std::printf("\n");
  gate("5 saturation (absolute)", worstSat, 1.0e-10);
  {
    Query const q = classifyMu(ri, 0.0, 120000.0, 0.5, 1.0);
    std::vector<double> s;
    slopesAt(prof, 480.0, s);
    double const t8 = tauProvisional(s, prof, ri, 1.0, q.p, q.rObs, q.rE, 8);
    double const t16 = tauProvisional(s, prof, ri, 1.0, q.p, q.rObs, q.rE, 16);
    std::printf("  convergence indicator (z_obs=0,f=1,mu=0.5,lam=480,z_e=120000): "
                "|tau(16)-tau(8)| = %.3e (informational)\n",
                std::fabs(t16 - t8));
  }
  std::printf("  stable-kernel vs source rule (informational; T1 FP gate stays FAILED): "
              "max |alpha_kern - alpha_src| over z = 0..50025 at lam=480 = ");
  {
    std::vector<double> s;
    slopesAt(prof, 480.0, s);
    double worst = 0.0;
    for (int i = 0; i <= 50025; ++i) {
      double const z = double(i);
      worst = std::max(worst, std::fabs(alphaKernelZ(s, prof.hgt, z) - alphaSrcZ(abs, 480.0, z)));
    }
    std::printf("%.3e\n\n", worst);
  }

  // ---- fixed-emission-cosine comparison (SEPARATE from the saturation test)
  std::printf("  fixed EMISSION cosine comparison (mu_e = 0.5, z_obs = 0, f = 1, lambda = 480):\n");
  std::printf("  %-14s %-22s %-18s %-12s %-22s %-22s\n", "z_e [m]", "p [m]", "mu_obs",
              "status", "tau", "route min g [m]");
  double const zeC[2] = {H_TOP + 100.0, 120000.0};
  for (int i = 0; i < 2; ++i) {
    Query const q = classifyEmission(ri, 0.0, zeC[i], -0.5, 1.0);
    double tau = std::nan("");
    if (q.status == Status::Eligible) {
      std::vector<double> s;
      slopesAt(prof, 480.0, s);
      tau = tauProvisional(s, prof, ri, q.f, q.p, q.rObs, q.rE, 8);
    }
    std::printf("  %-14.17g %-22.17g %-18.17g %-12s %-22.17g %-22.17g\n", zeC[i], q.p,
                q.muObs, statusName(q.status), tau, q.gMin);
  }
  std::printf("  (descriptive: different invariants p give different optical depths; this is\n");
  std::printf("   not a saturation identity)\n\n");

  // ======================= 6. summary ======================================
  std::printf("--- 6. summary ---\n");
  std::printf("  grid counts:");
  for (int k = 0; k < 7; ++k) std::printf(" %s=%ld", statusName(Status(k)), counts[k]);
  std::printf("\n");
  std::printf("  saturation cases=%ld worst=%.17g\n", nSat, worstSat);
  std::printf("  gates: %d evaluated, %d failed\n", nGate, nFail);
  std::printf("%s\n", (nFail == 0) ? "ALL GATES PASS" : "STOP: a gate FAILED -- reported above");
  std::printf("T1 floating-point equivalence gate: STILL FAILED (not waived by this run).\n");
  std::printf("=== end of T2 (2f-b) ===\n");
  return (nFail == 0) ? 0 : 1;
}

