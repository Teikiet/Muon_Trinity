// stage9b_abs_b.cpp -- Task E2, checkpointed absorption scan.
//
//   ./stage9b_abs_b refcheck <profile>   Checkpoint 1: do the references agree?
//   ./stage9b_abs_b scan     <profile>   Checkpoints 2-4: the scan, the (a)/(b) split, timing
//
// Production path under test: Cherenkov.inl:96-113, the +/-(dh/2) finite difference of
// AtmosphericAbsorption::getOpticalDepth with dh = 50 m.
// Reference: the native-knot slope of the profile's own cumulative tau(h), integrated
// along the same curved path (two parameterisations: w = sqrt(z-R) and a z-march).
//
// No production edits. New diagnostics only.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0;
static double gTauLookup = 0.0, gTauSlope = 0.0; // accumulated seconds
static long gLookupCalls = 0, gSteps = 0;

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
static double nPhase(RefractiveIndexTable const& ri, double zAbs, double f) {
  return 1.0 + ri.N(zAbs) * f;
}
static double bCrit(RefractiveIndexTable const& ri, double f) { return nPhase(ri, R, f) * R; }

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
// exact slope of the native piecewise-linear tau(h)
static double alphaRef(Profile const& p, int il, double h) {
  std::vector<double> const& H = p.hgt;
  std::vector<double> const& T = p.tau[il];
  if (h <= H.front()) return (T[1] - T[0]) / (H[1] - H[0]);
  if (h >= H.back()) return (T[H.size() - 1] - T[H.size() - 2]) / (H.back() - H[H.size() - 2]);
  std::size_t j = 0;
  while (j + 2 < H.size() && H[j + 1] <= h) ++j;
  return (T[j + 1] - T[j]) / (H[j + 1] - H[j]);
}
// the production rule (Cherenkov.inl:96-113), via the production class
static double alphaProd(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                        double hHi) {
  double const hl = std::max(hLo, h - 25.0), hh = std::min(hHi, h + 25.0);
  double const adh = hh - hl;
  if (adh <= 0.0) return 0.0;
  auto const q0 = std::chrono::steady_clock::now();
  double const t1 = abs.getOpticalDepth(lam, hh);
  double const t0 = abs.getOpticalDepth(lam, hl);
  auto const q1 = std::chrono::steady_clock::now();
  gTauLookup += std::chrono::duration<double>(q1 - q0).count();
  gLookupCalls += 2;
  return std::max(0.0, (t1 - t0) / adh);
}

// tau along the ray.  mode 0 = production alpha, mode 1 = native-knot slope.
// param 0 = w-march with step dw, param 1 = z-march with step dz.
static double tauInt(AtmosphericAbsorption const& abs, Profile const& prof,
                     RefractiveIndexTable const& ri, double lam, int ilp, double p, double f,
                     double hTop, double step, int mode, int param, double hLo, double hHi) {
  double tau = 0.0;
  if (param == 0) {
    double const wTop = std::sqrt(hTop);
    for (double w = 0.0; w < wTop - 1e-12;) {
      double const st = std::min(step, wTop - w);
      double const wn = w + st;
      if (!(wn > w)) break;
      double const wm = w + 0.5 * st;
      double const z = R + wm * wm;
      double const h = wm * wm;
      double const s = p / (nPhase(ri, z, f) * z);
      if (s >= 1.0) break;
      double const ds = (2.0 * wm) * st / std::sqrt(std::max(1e-300, 1.0 - s * s));
      if (mode == 0) {
        tau += alphaProd(abs, lam, h, hLo, hHi) * ds;
      } else {
        auto const q0 = std::chrono::steady_clock::now();
        double const a = alphaRef(prof, ilp, h);
        auto const q1 = std::chrono::steady_clock::now();
        gTauSlope += std::chrono::duration<double>(q1 - q0).count();
        tau += a * ds;
      }
      ++gSteps;
      w = wn;
    }
  } else {
    double z = R;
    double const zTop = R + hTop;
    while (z < zTop - 1e-12) {
      double const st = std::min(step, zTop - z);
      double const zn = z + st;
      if (!(zn > z)) break;
      double const zm = z + 0.5 * st;
      double const h = zm - R;
      double const s = p / (nPhase(ri, zm, f) * zm);
      if (s >= 1.0) break;
      double const ds = st / std::sqrt(std::max(1e-300, 1.0 - s * s));
      if (mode == 0) {
        tau += alphaProd(abs, lam, h, hLo, hHi) * ds;
      } else {
        tau += alphaRef(prof, ilp, h) * ds;
      }
      ++gSteps;
      z = zn;
    }
  }
  return tau;
}

// Path-length march: RK4 in s along the same curved path, same integrand.
//   dz/ds = cos(theta),  dtheta/ds = -sin(theta) * (1/r + n'/n),  tau += alpha(h) ds
// Marching down from the emission to R, as the test's referenceRay does; tau is the same
// in either direction.  mode 0 = production alpha, mode 1 = native-knot slope.
static double tauIntS(AtmosphericAbsorption const& abs, Profile const& prof,
                      RefractiveIndexTable const& ri, double lam, int ilp, double p, double f,
                      double hTop, double ds, int mode, double hLo, double hHi, long* nSteps,
                      std::vector<double>* rec = nullptr) {
  double z = R + hTop;
  double const nTop = 1.0 + ri.N(z) * f;
  double const sinT = std::min(1.0, p / (nTop * z));
  double th = std::acos(-std::sqrt(std::max(0.0, 1.0 - sinT * sinT))); // downgoing
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
    if (rec) {
      rec->push_back(h);
      rec->push_back(ds);
    }
    ++n;
  }
  if (nSteps) *nSteps = n;
  return tau;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::printf("usage: %s <refcheck|scan> <profile .dat>\n", argv[0]);
    return 2;
  }
  std::string const mode = argv[1];
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[2]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  Profile prof;
  if (!readProfile(argv[2], prof)) {
    std::printf("STOP: cannot re-read the profile's native knots\n");
    return 1;
  }
  auto const hr = abs.getHeightRange();
  RefractiveIndexTable ri = makeRI();
  double const lam = 300.0, hE = 46000.0;
  double const f = Fdisp(lam);
  int const ilp = nearestLam(prof, lam);

  if (mode == "refcheck") {
    std::printf("=== Task E2 Checkpoint 1: reference convergence ===\n");
    std::printf("3 rays: lambda=%.0f nm, h=%.0f m, mu in {1e-3, 2e-2, 0.6}\n", lam, hE);
    std::printf("integrand: the native-knot slope of the profile's own cumulative tau(h)\n\n");
    std::printf("   %-9s %-16s %-16s %-16s %-12s %-12s\n", "mu", "w(dw=2e-3)", "w(dw=1e-3)",
                "z(dz=1m)", "|w2-w1|", "|w1-z|");
    double const mus3[3] = {1e-3, 2e-2, 0.6};
    double wWorst = 0.0, wzWorst = 0.0, sWorst = 0.0, sRWorst = 0.0;
    for (int q = 0; q < 3; ++q) {
      double const p = std::sqrt(std::max(0.0, 1.0 - mus3[q] * mus3[q])) * bCrit(ri, f);
      double const w2 =
          tauInt(abs, prof, ri, lam, ilp, p, f, hE, 2.0e-3, 1, 0, hr.first, hr.second);
      double const w1 =
          tauInt(abs, prof, ri, lam, ilp, p, f, hE, 1.0e-3, 1, 0, hr.first, hr.second);
      double const zz = tauInt(abs, prof, ri, lam, ilp, p, f, hE, 1.0, 1, 1, hr.first, hr.second);
      double const d21 = std::fabs(w2 - w1), dz = std::fabs(w1 - zz);
      if (d21 > wWorst) wWorst = d21;
      if (dz > wzWorst) wzWorst = dz;
      std::printf("   %-9.2e %-16.8f %-16.8f %-16.8f %-12.3e %-12.3e\n", mus3[q], w2, w1, zz, d21,
                  dz);
      // Checkpoint 1b: the path-length march, plus the Richardson value of w
      double const rich = 2.0 * w1 - w2;
      auto const r0 = std::chrono::steady_clock::now();
      long n1 = 0, n2 = 0;
      double const s1 = tauIntS(abs, prof, ri, lam, ilp, p, f, hE, 1.0, 1, hr.first, hr.second, &n1);
      auto const r1 = std::chrono::steady_clock::now();
      double const s2 = tauIntS(abs, prof, ri, lam, ilp, p, f, hE, 0.5, 1, hr.first, hr.second, &n2);
      auto const r2 = std::chrono::steady_clock::now();
      double const dss = std::fabs(s1 - s2), dsr = std::fabs(s2 - rich);
      if (dss > sWorst) sWorst = dss;
      if (dsr > sRWorst) sRWorst = dsr;
      std::printf("      s(1)=%.8f s(0.5)=%.8f R=%.8f |s1-s.5|=%.3e |s.5-R|=%.3e steps=%ld/%ld "
                  "t=%.1f/%.1f s\n",
                  s1, s2, rich, dss, dsr, n1, n2,
                  std::chrono::duration<double>(r1 - r0).count(),
                  std::chrono::duration<double>(r2 - r1).count());
    }
    std::printf("\nlargest |w(2e-3) - w(1e-3)| = %.4e\n", wWorst);
    std::printf("largest |w - z|            = %.4e\n", wzWorst);
    std::printf("largest |s(1) - s(0.5)|     = %.4e\n", sWorst);
    std::printf("largest |s(0.5) - R|        = %.4e\n", sRWorst);
    std::printf("gate 3.3e-4 (1/3 of the 1e-3 tolerance): %s\n",
                (sWorst <= 3.3e-4 && sRWorst <= 3.3e-4)
                    ? "PASS -- use s(ds=1 m) as the reference"
                    : "STOP -- the path-length reference is not converged");
    return 0;
  }

  // ---- 2d Step 1: the production step sequence on the three sanity rays --------
  if (mode == "seq2d") {
    struct C {
      double lam, mu, h;
    };
    C const cs[3] = {{300.0, 0.6, 46000.0}, {480.0, 0.1, 1000.0}, {480.0, 0.1, 975.0}};
    double const dsP = 400.0; // proxy: min(H_scale/20, 1000) with H_scale ~ 8 km
    std::printf("=== 2d Step 1: production step sequence (proxy ds = %.0f m) ===\n\n", dsP);
    for (int q = 0; q < 3; ++q) {
      double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
      int const ilq = nearestLam(prof, lm);
      double const ff = Fdisp(lm);
      double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
      long n1 = 0;
      double const ref =
          tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &n1);
      std::vector<double> rc;
      long np = 0;
      tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsP, 1, hr.first, hr.second, &np, &rc);
      long const ns = long(rc.size() / 2);
      double sumds = 0.0, v1 = 0.0, v2 = 0.0;
      for (long k = 0; k < ns; ++k) {
        double const h = rc[2 * k], ds = rc[2 * k + 1];
        sumds += ds;
        v1 += alphaProd(abs, lm, h, hr.first, hr.second) * ds;
        v2 += alphaRef(prof, ilq, h) * ds;
      }
      long n3 = 0;
      double const v3 =
          tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsP / 100.0, 1, hr.first, hr.second, &n3);
      std::printf("lambda=%.0f mu=%.2e h=%.0f\n", lm, mu, hh);
      std::printf("   steps=%ld  sum(ds)=%.6f  1m path=%.6f (n1=%ld)  ratio=%.9f  "
                  "path/steps=%.6f  mean ds=%.6f\n",
                  np, sumds, double(n1) * 1.0, n1, sumds / (double(n1) * 1.0), sumds / double(np),
                  sumds / double(np));
      std::printf("   ds    first3 = %.3f %.3f %.3f   last3 = %.3f %.3f %.3f\n", rc[1], rc[3],
                  rc[5], rc[2 * ns - 1], rc[2 * ns - 3], rc[2 * ns - 5]);
      std::printf("   h_mid first3 = %.4f %.4f %.4f   last3 = %.4f %.4f %.4f\n", rc[0], rc[2],
                  rc[4], rc[2 * ns - 2], rc[2 * ns - 4], rc[2 * ns - 6]);
      std::printf("   tau_ref=%.8f   (i) kext_prod=%.8f   (ii) alphaRef=%.8f   (iii) exact=%.8f\n",
                  ref, v1, v2, v3);
      std::fflush(stdout);
    }
    std::printf("\n=== end of 2d Step 1 ===\n");
    return 0;
  }

  // ---- Step 1 of the 2d brief: the reference cache -----------------------------
  if (mode == "cache" || mode == "cache_check") {
    char const* fn = "ref_cache.bin";
    double const hC[7] = {100.0, 990.0, 1000.0, 5000.0, 20000.0, 46000.0, 50000.0};
    double const mC[5] = {2.0e-4, 1.0e-3, 2.0e-2, 1.0e-1, 0.6};
    std::vector<double> lC;
    for (std::size_t i = 0; i < prof.lam.size(); ++i)
      if (prof.lam[i] >= 300.0 - 1e-9 && prof.lam[i] <= 700.0 + 1e-9 && (i % 4) == 0)
        lC.push_back(prof.lam[i]);
    if (mode == "cache") {
      FILE* f = std::fopen(fn, "wb");
      if (!f) {
        std::printf("STOP: cannot open %s for writing\n", fn);
        return 1;
      }
      long nrec = 0;
      for (std::size_t i = 0; i < lC.size(); ++i)
        for (int a = 0; a < 5; ++a)
          for (int b = 0; b < 7; ++b) {
            double const lm = lC[i], mu = mC[a], hh = hC[b];
            int const ilq = nearestLam(prof, lm);
            double const ff = Fdisp(lm);
            double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
            std::vector<double> rc;
            long n = 0;
            double const tr =
                tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &n, &rc);
            double path = 0.0;
            for (std::size_t k = 1; k < rc.size(); k += 2) path += rc[k];
            long const ns = long(rc.size() / 2);
            double const hdr[6] = {lm, mu, hh, tr, path, (ns > 0) ? 1.0 : 0.0};
            std::fwrite(hdr, sizeof(double), 6, f);
            std::fwrite(&ns, sizeof(long), 1, f);
            if (ns) std::fwrite(rc.data(), sizeof(double), rc.size(), f);
            ++nrec;
          }
      std::fclose(f);
      std::ifstream in(fn, std::ios::binary | std::ios::ate);
      long long sz = -1;
      if (in) sz = (long long)in.tellg();
      std::printf("cache: records=%ld  file size=%lld bytes\n", nrec, sz);
      return 0;
    }
    FILE* f = std::fopen(fn, "rb");
    if (!f) {
      std::printf("STOP: cannot open %s for reading\n", fn);
      return 1;
    }
    long nrec = 0, nbad = 0;
    double worst = 0.0;
    for (;;) {
      double hdr[6];
      long ns = 0;
      if (std::fread(hdr, sizeof(double), 6, f) != 6) break;
      if (std::fread(&ns, sizeof(long), 1, f) != 1) break;
      std::vector<double> rc(2 * std::size_t(ns));
      if (ns && std::fread(rc.data(), sizeof(double), rc.size(), f) != rc.size()) break;
      int const ilq = nearestLam(prof, hdr[0]);
      double t = 0.0;
      for (long k = 0; k < ns; ++k) t += alphaRef(prof, ilq, rc[2 * k]) * rc[2 * k + 1];
      double const rel =
          (hdr[3] != 0.0) ? std::fabs(t - hdr[3]) / std::fabs(hdr[3]) : std::fabs(t);
      if (rel > worst) worst = rel;
      if (!(rel <= 1e-12)) ++nbad;
      ++nrec;
    }
    std::fclose(f);
    std::printf("cache_check: records=%ld  bad=%ld  worst relative=%.3e  gate 1e-12: %s\n", nrec,
                nbad, worst, (nbad == 0) ? "PASS" : "FAIL");
    return 0;
  }

  // ---- Steps 2-3 of the 2c brief: the five columns on the sanity rays ----------
  if (mode == "chk2c") {
    struct C {
      double lam, mu, h;
    };
    C const cs[3] = {{300.0, 0.6, 46000.0}, {480.0, 0.1, 1000.0}, {480.0, 0.1, 975.0}};
    // The tracer's own step (Cherenkov.inl:608-615): dx = min(H_scale/20, 1000) m, or
    // H_scale/100 for near-radial rays, with H_scale = computeLocalScaleHeight(h_cur).
    // Proxy here: H_scale ~ 8 km, so the path step is 400 m.
    double const dsTr = 8000.0 / 20.0;
    std::printf("=== E2 2c Steps 2-3: the five columns on the sanity rays ===\n");
    std::printf("tracer path-step proxy ds = %.1f m  (Cherenkov.inl:608-615, H_scale ~ 8 km)\n\n", dsTr);
    for (int q = 0; q < 3; ++q) {
      double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
      int const ilq = nearestLam(prof, lm);
      double const ff = Fdisp(lm);
      double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
      long n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0;
      double const ref = tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &n1);
      double const rule1 =
          tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 0, hr.first, hr.second, &n2);
      double const prodc =
          tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsTr, 0, hr.first, hr.second, &n3);
      double const tabc =
          tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsTr, 1, hr.first, hr.second, &n4);
      double const tab1 =
          tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &n5);
      double const da = rule1 - ref, db = prodc - rule1, dtot = prodc - ref;
      double const dtab = tabc - ref, dtab1m = tab1 - ref;
      double const path = double(n3) * dsTr;
      std::printf("lambda=%.0f mu=%.2e h=%.0f\n", lm, mu, hh);
      std::printf("   tau_ref=%.8f  da=%.3e  db=%.3e  dtot=%.3e  dtab=%.3e  dtab1m=%.3e\n", ref, da, db,
                  dtot, dtab, dtab1m);
      std::printf("   steps ref=%ld rule1m=%ld prod_curved=%ld tab_curved=%ld tab_1m=%ld   "
                  "tracer path=%.0f m  mean step=%.1f m\n",
                  n1, n2, n3, n4, n5, path, path / double(std::max(1L, n3)));
      std::fflush(stdout);
    }
    std::printf("\n=== end of Steps 2-3 ===\n");
    return 0;
  }

  // ---- Step 1: knot-offset test on the 5 worst rays of cp2_rule_on_1m.out ----
  if (mode == "knot") {
    struct R {
      double lam, mu;
    };
    R const rays[5] = {{480.0, 1.0e-1},
                       {680.0, 2.0e-2},
                       {460.0, 1.0e-1},
                       {500.0, 1.0e-1},
                       {660.0, 2.0e-2}};
    double const hs2[7] = {950.0, 975.0, 1000.0, 1025.0, 1050.0, 1100.0, 1500.0};
    std::printf("=== Step 1: knot-offset test (production rule on 1 m reference steps) ===\n");
    std::printf("%-7s %-9s %-8s %-17s %-13s %-13s\n", "lambda", "mu", "h", "tau_ref", "dtau",
                "dtau/max(1,tau)");
    for (int q = 0; q < 5; ++q) {
      double const lm = rays[q].lam, mu = rays[q].mu;
      int const ilq = nearestLam(prof, lm);
      double const ff = Fdisp(lm);
      double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
      for (int b = 0; b < 7; ++b) {
        double const hh = hs2[b];
        long nR = 0, nP = 0;
        double const tr =
            tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &nR);
        double const tp =
            tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 0, hr.first, hr.second, &nP);
        double const d = tp - tr;
        std::printf("%-7.0f %-9.2e %-8.0f %-17.8f %+-13.3e %+-13.3e\n", lm, mu, hh, tr, d,
                    d / std::max(1.0, tr));
        std::fflush(stdout);
      }
    }
    std::printf("\nnative-knot slope (alphaRef) by layer:\n");
    for (int k = 0; k < 2; ++k) {
      double const lm = (k == 0) ? 480.0 : 680.0;
      int const ilq = nearestLam(prof, lm);
      std::printf("  lambda=%.0f:  h=500:%.6e  h=950:%.6e  h=1050:%.6e  h=1500:%.6e  h=5000:%.6e\n",
                  lm, alphaRef(prof, ilq, 500.0), alphaRef(prof, ilq, 950.0),
                  alphaRef(prof, ilq, 1050.0), alphaRef(prof, ilq, 1500.0),
                  alphaRef(prof, ilq, 5000.0));
      std::fflush(stdout);
    }
    std::printf("\n=== end of Step 1 ===\n");
    return 0;
  }

  // ---- scan -----------------------------------------------------------------
  auto const t0 = std::chrono::steady_clock::now();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  auto const t1 = std::chrono::steady_clock::now();
  double const buildS = std::chrono::duration<double>(t1 - t0).count();
  std::printf("=== Task E2 Checkpoint 2: absorption scan (production path) ===\n");
  std::printf("geometry table build %.1f s (valid=%d)\n", buildS, int(tab.valid()));
  if (!tab.valid()) return 1;
  std::vector<double> lams;
  for (std::size_t i = 0; i < prof.lam.size(); ++i)
    if (prof.lam[i] >= 300.0 - 1e-9 && prof.lam[i] <= 700.0 + 1e-9 && (i % 4) == 0)
      lams.push_back(prof.lam[i]);
  double const mus[5] = {2.0e-4, 1.0e-3, 2.0e-2, 1.0e-1, 0.6};
  double const hs[6] = {100.0, 1000.0, 5000.0, 20000.0, 46000.0, 50000.0};
  std::printf("lambda nodes used: %zu (every 4th file knot from 300 to 700 nm)\n\n", lams.size());

  int tot[5][6] = {{0}}, fail[5][6] = {{0}}, excl[5][6] = {{0}};
  double worst[5][6];
  for (int a = 0; a < 5; ++a)
    for (int b = 0; b < 6; ++b) worst[a][b] = 0.0;
  double gW = 0.0, gLam = 0, gMu = 0, gH = 0;
  double rec[5] = {0.0}, recA[5] = {0.0}, recB[5] = {0.0}, recL[5] = {0.0};
  double recLam[5] = {0.0}, recMu[5] = {0.0}, recH[5] = {0.0};
  int nHit = 0, nExcl = 0;
  double tProd = 0.0, tRef = 0.0;

  for (std::size_t il = 0; il < lams.size(); ++il) {
    double const lm = lams[il];
    int const ilq = nearestLam(prof, lm);
    double const ff = Fdisp(lm);
    for (int a = 0; a < 5; ++a) {
      double const mu = mus[a];
      double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
      for (int b = 0; b < 6; ++b) {
        double const hh = hs[b];
        double const re = R + hh;
        double const sinT = p / (nPhase(ri, re, ff) * re);
        if (sinT >= 1.0) {
          ++excl[a][b];
          ++nExcl;
          continue;
        }
        double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
        CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, ff, ff);
        if (r.status != CurvedRayTable::Status::Hit) {
          ++excl[a][b]; // mu = 2e-4 grazing rejections live in their own column
          ++nExcl;
          continue;
        }
        ++nHit;
        long nsP = 0, nsR = 0;
        auto const q0 = std::chrono::steady_clock::now();
        double const tp =
            tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 0, hr.first, hr.second, &nsP);
        auto const q1 = std::chrono::steady_clock::now();
        double const tr =
            tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 1.0, 1, hr.first, hr.second, &nsR);
        auto const q2 = std::chrono::steady_clock::now();
        tProd += std::chrono::duration<double>(q1 - q0).count();
        tRef += std::chrono::duration<double>(q2 - q1).count();
        // RELATIVE tolerance: a ray passes if |dtau| <= 1e-3 * max(1, tau_ref).  An absolute
        // limit is meaningless at tau ~ 75, where the transmission is e^-75.
        double const d = std::fabs(tp - tr);
        double const lim = 1e-3 * std::max(1.0, tr);
        double const rel = d / std::max(1.0, tr);
        double const dT = std::exp(-tr) * d;
        ++tot[a][b];
        if (d > lim) ++fail[a][b];
        if (rel > worst[a][b]) worst[a][b] = rel;
        if (rel > gW) {
          gW = rel;
          gLam = lm;
          gMu = mu;
          gH = hh;
        }
        if (rel > rec[4]) {
          int q = 4;
          while (q > 0 && rel > rec[q - 1]) {
            rec[q] = rec[q - 1];
            recA[q] = recA[q - 1];
            recB[q] = recB[q - 1];
            recL[q] = recL[q - 1];
            recLam[q] = recLam[q - 1];
            recMu[q] = recMu[q - 1];
            recH[q] = recH[q - 1];
            --q;
          }
          rec[q] = rel;
          recA[q] = rel; // (a): same path, same steps, only the alpha rule differs
          long nsB = 0;
          double const tr2 =
              tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, 2.0, 1, hr.first, hr.second, &nsB);
          recB[q] = std::fabs(tr2 - tr) / std::max(1.0, tr); // (b): step size only
          recL[q] = double(nsR) * 1.0;                       // path length of the reference
          recLam[q] = lm;
          recMu[q] = mu;
          recH[q] = hh;
        }
      }
    }
  }

  char const* ml[5] = {"2e-4", "1e-3", "2e-2", "1e-1", "0.6"};
  std::printf("fails / total (worst relative |dtau|) per mu x h; E = excluded (grazing rejection)\n");
  std::printf("%-8s", "mu \\ h");
  for (int b = 0; b < 6; ++b) std::printf(" h=%-13.0f", hs[b]);
  std::printf("\n");
  for (int a = 0; a < 5; ++a) {
    std::printf("%-8s", ml[a]);
    for (int b = 0; b < 6; ++b) {
      char buf[32];
      if (tot[a][b] == 0) std::snprintf(buf, sizeof(buf), "E(%d)", excl[a][b]);
      else std::snprintf(buf, sizeof(buf), "%d/%-2d %.1e", fail[a][b], tot[a][b], worst[a][b]);
      std::printf(" %-15s", buf);
    }
    std::printf("\n");
  }
  std::printf("\nhits %d, excluded %d\n", nHit, nExcl);
  std::printf("worst relative |dtau|/max(1,tau_ref) = %.4e  at lambda=%.0f mu=%.2e h=%.0f\n", gW, gLam,
              gMu, gH);
  std::printf("(relative limit 1e-3; the absolute |dtau| and |dT| of that ray are in the split below)\n");
  std::printf("\nCheckpoint 3: the 5 worst rays split into (a) alpha-rule and (b) step-size:\n");
  std::printf("   %-10s %-10s %-8s %-11s %-11s %-11s %-12s\n", "|dtau|", "lambda", "mu", "h", "(a)",
              "(b)", "path L[m]");
  for (int q = 0; q < 5; ++q)
    std::printf("   %-10.3e %-10.0f %-8.2e %-11.0f %-11.3e %-11.3e %-12.0f\n", rec[q], recLam[q],
                recMu[q], recH[q], recA[q], recB[q], recL[q]);
  std::printf("\nCheckpoint 4: timing (per integration step, z-march at dr=1 m)\n");
  std::printf("   production rule %.4f us/step (%ld calls, %.1f calls/step)\n",
              gSteps ? gTauLookup / double(gSteps) * 1e6 : 0.0, gLookupCalls,
              gSteps ? double(gLookupCalls) / double(gSteps) : 0.0);
  std::printf("   direct slope    %.4f us/step\n", gSteps ? gTauSlope / double(gSteps) * 1e6 : 0.0);
  std::printf("   wall: production %.1f s, reference %.1f s, geometry build %.1f s\n", tProd, tRef,
              buildS);
  std::printf("\n=== end of Task E2 ===\n");
  return 0;
}
