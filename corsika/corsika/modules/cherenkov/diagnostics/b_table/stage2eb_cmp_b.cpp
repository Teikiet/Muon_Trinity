// stage2eb_cmp_b.cpp -- Checkpoint 2e-b: extracted tracer vs the existing diagnostic.
//
//   ./stage2eb_cmp_b <profile .dat>
//
// Compares, for the three fixed reference rays,
//   tau_ref         = tauIntS(ds = 1 m, mode 1)   native-knot slope   [diagnostic]
//   tau_prod_curved = tauIntS(ds = 400 m, mode 0) production alpha     [diagnostic "db" column]
//   tau_tr          = refdiag::RefTracer::trace() production tracer, lifted verbatim
//
// Inputs are the diagnostic's own, unchanged: same RefractiveIndexTable, same
// Fdisp, same alphaProd (Cherenkov.inl:96-113), same p = sin(theta)*bCrit.
// Plane for the tracer: flat ground, point (0,0,R), normal (0,0,1), emission at
// (0,0,R+h) -- the flat-plane analogue of the diagnostic's march to r = R.
//
// No production edits.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include "ref_tracer_b.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0;

static double Nshape(double h) {
  return 2.20e-4 * std::exp(-h / 7500.0) + 0.55e-4 * std::exp(-h / 1600.0);
}
static double Fdisp(double lam_nm) {
  double const l = lam_nm / 1000.0;
  return (1.0 + 0.0032 / (l * l)) / (1.0 + 0.0032 / 0.16);
}
static RefractiveIndexTable makeRI() {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= H_TOP + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(std::log(Nshape(h)));
  }
  return RefractiveIndexTable(alt, lnN, R);
}
static double bCrit(RefractiveIndexTable const& ri, double f) { return (1.0 + ri.N(R) * f) * R; }

struct Profile {
  std::vector<double> lam, hgt;
  std::vector<std::vector<double> > tau;
};
static bool readProfile(std::string const& fn, Profile& p) {
  std::ifstream f(fn.c_str());
  if (!f) return false;
  int nl = 0, nh = 0;
  f >> nl >> nh;
  if (nl <= 0 || nh <= 0) return false;
  p.lam.resize(nl);
  p.hgt.resize(nh);
  for (int i = 0; i < nl; ++i) f >> p.lam[i];
  for (int i = 0; i < nh; ++i) f >> p.hgt[i];
  p.tau.assign(nl, std::vector<double>(nh, 0.0));
  for (int i = 0; i < nl; ++i)
    for (int j = 0; j < nh; ++j) f >> p.tau[i][j];
  return true;
}
static int nearestLam(Profile const& p, double lam) {
  int best = 0;
  double bd = 1e30;
  for (std::size_t i = 0; i < p.lam.size(); ++i) {
    double const d = std::fabs(p.lam[i] - lam);
    if (d < bd) {
      bd = d;
      best = int(i);
    }
  }
  return best;
}
static double alphaRef(Profile const& p, int il, double h) {
  std::vector<double> const& H = p.hgt;
  std::vector<double> const& T = p.tau[il];
  if (h <= H.front()) return (T[1] - T[0]) / (H[1] - H[0]);
  if (h >= H.back())
    return (T[H.size() - 1] - T[H.size() - 2]) / (H.back() - H[H.size() - 2]);
  std::size_t j = 0;
  while (j + 2 < H.size() && H[j + 1] <= h) ++j;
  return (T[j + 1] - T[j]) / (H[j + 1] - H[j]);
}

// production rule, Cherenkov.inl:96-113, via the production class
static double alphaProd(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                        double hHi) {
  double const hl = std::max(hLo, h - 25.0), hh = std::min(hHi, h + 25.0);
  double const adh = hh - hl;
  if (adh <= 0.0) return 0.0;
  double const t1 = abs.getOpticalDepth(lam, hh);
  double const t0 = abs.getOpticalDepth(lam, hl);
  return std::max(0.0, (t1 - t0) / adh);
}

// diagnostic march: RK4 in path length s from the emission down to r = R.
static double tauIntS(AtmosphericAbsorption const& abs, Profile const& prof,
                      RefractiveIndexTable const& ri, double lam, int ilp, double p, double f,
                      double hTop, double ds, int mode, double hLo, double hHi, long* nSteps) {
  double z = R + hTop;
  double const nTop = 1.0 + ri.N(z) * f;
  double const sinT = std::min(1.0, p / (nTop * z));
  double th = std::acos(-std::sqrt(std::max(0.0, 1.0 - sinT * sinT)));
  double tau = 0.0;
  long n = 0;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const nn = 1.0 + ri.N(rr) * f;
    double const np = f * ri.dNdz(rr);
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c;
    dt = -s * (1.0 / rr + (nn > 1e-12 ? np / nn : 0.0));
  };
  while (z > R + 1e-9 && n < 200000000L) {
    double dr1, dt1, dr2, dt2, dr3, dt3, dr4, dt4;
    der(z, th, dr1, dt1);
    der(z + 0.5 * ds * dr1, th + 0.5 * ds * dt1, dr2, dt2);
    der(z + 0.5 * ds * dr2, th + 0.5 * ds * dt2, dr3, dt3);
    der(z + ds * dr3, th + ds * dt3, dr4, dt4);
    double const zPrev = z;
    z += ds / 6.0 * (dr1 + 2 * dr2 + 2 * dr3 + dr4);
    th += ds / 6.0 * (dt1 + 2 * dt2 + 2 * dt3 + dt4);
    double const h = std::max(0.0, 0.5 * (zPrev + z) - R);
    tau += (mode == 0 ? alphaProd(abs, lam, h, hLo, hHi) : alphaRef(prof, ilp, h)) * ds;
    ++n;
  }
  if (nSteps) *nSteps = n;
  return tau;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <profile .dat>\n", argv[0]);
    return 2;
  }
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[1]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  Profile prof;
  if (!readProfile(argv[1], prof)) {
    std::printf("STOP: cannot re-read the native knots\n");
    return 1;
  }
  auto const hr = abs.getHeightRange();
  RefractiveIndexTable ri = makeRI();

  struct C {
    double lam, mu, h;
  };
  C const cs[3] = {{300.0, 0.6, 46000.0}, {480.0, 0.1, 1000.0}, {480.0, 0.1, 975.0}};
  double const dsTr = 8000.0 / 20.0;  // diagnostic proxy step

  std::printf("=== Checkpoint 2e-b: extracted tracer vs the diagnostic ===\n");
  std::printf("tracer plane: point (0,0,%.1f) normal (0,0,1); emission (0,0,R+h)\n", R);
  std::printf("survivalThreshold_in = 0.0 (exp(-tau) < 0.0 is never true)\n\n");

  for (int q = 0; q < 3; ++q) {
    double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
    int const ilq = nearestLam(prof, lm);
    double const ff = Fdisp(lm);
    double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);

    long n1 = 0, n3 = 0, n2 = 0;
    double const ref = tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &n1);
    double const rule1 =
        tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 0, hr.first, hr.second, &n2);
    double const prodc =
        tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsTr, 0, hr.first, hr.second, &n3);

    refdiag::RefTracer tr;
    tr.earthRadiusM_ = R;
    tr.enableAbsorption_ = true;
    tr.refractiveIndexFunc_ = [&](double h_m, double lam_nm) {
      double const f = Fdisp(lam_nm);
      double const n = 1.0 + ri.N(R + h_m) * f;
      return std::make_pair(n, n);  // equal-f convention, as the diagnostic uses
    };
    tr.extinctionCoeffFunc_ = [&](double h_m, double lam_nm) {
      return alphaProd(abs, lam_nm, h_m, hr.first, hr.second);
    };

    double const sinT = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    refdiag::RefTracer::Result const r =
        tr.trace(0.0, 0.0, R + hh, sinT, 0.0, -mu, 0.0, 0.0, R, 0.0, 0.0, 1.0, lm, 0.0, 0.0);

    double const tau_tr = r.accumulated_tau;
    double const den = std::max(std::max(std::fabs(ref), std::fabs(tau_tr)), 1e-30);
    double const e = std::fabs(tau_tr - ref) / den;

    std::printf("ray %d: lambda=%.0f mu=%.2e h=%.0f\n", q, lm, mu, hh);
    std::printf("   tau_ref=%.10f  tau_prod_curved=%.10f  tau_tracer=%.10f  E_tau=%.3e  branch=%d hit=%d\n",
                ref, prodc, tau_tr, e, r.branch, int(r.hit));
    std::printf("   da=rule1m-ref=%.6e  db=prod_curved-rule1m=%.6e  dtot=prod_curved-ref=%.6e\n",
                rule1 - ref, prodc - rule1, prodc - ref);
    std::printf("   steps: ref=%ld rule1m=%ld prod_curved=%ld tracer=%d\n", n1, n2, n3, r.steps);
    std::printf("   path: diag(prod_curved)=%.1f m  tracer=%.1f m  d=%.3f m\n",
                double(n3) * dsTr, r.geometric_path, r.geometric_path - double(n3) * dsTr);
    std::printf("   tracer: turning=%d tau_main=%.10f near_seg=%.3f near_sub=%d cross_alpha=%.6f "
                "cross_ds=%.3f cross_kext=%.6e\n",
                r.turning, r.tau_main, r.near_seg, r.near_sub, r.cross_alpha, r.cross_ds,
                r.cross_kext);
    if (r.branch == 1) {
      double const unclipped = (1.0 - r.cross_alpha) * r.cross_kext * r.cross_ds;
      std::printf("   crossing: (1-alpha)*k*ds=%.6e   tau_tracer-tau_main=%.6e   k*ds=%.6e\n",
                  unclipped, tau_tr - r.tau_main, r.cross_kext * r.cross_ds);
    } else {
      std::printf("   crossing estimate: NOT APPLICABLE (branch %d)\n", r.branch);
    }
    std::fflush(stdout);
  }
  std::printf("\n=== end of 2e-b ===\n");
  return 0;
}
