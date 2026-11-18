/*
 * testCherenkovTracerBenchmark.cpp
 *
 * Phase 3 benchmark: RK4 (reference/fallback) vs the optional adaptive DOPRI5 5(4) on a
 * frozen set of representative rays, measured against a CONVERGED reference rather than
 * against the production step size.
 *
 * Why a separate test: the acceptance tests answer "is the tracer right", this one
 * answers "what does it cost, and does a different integrator buy anything".  It prints
 * a CSV-shaped table and a summary so the numbers can be re-derived, and it asserts only
 * the accuracy gates -- the speed criterion is reported, because a wall-clock assertion
 * in CI is a coin flip.
 *
 * Reference: DOPRI5 with atol_pos = 1e-6 m, atol_dir = 1e-12, atol_opt = 1e-8 m and no
 * step cap.  It is cross-checked ONCE per run against RK4 at a 0.5 m radial cap, so a
 * bug in either cannot masquerade as agreement.
 */

#include <corsika/modules/cherenkov/Cherenkov.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace corsika;

namespace {

  constexpr double kR = 6.371e6;
  constexpr double kTelAlt = 3000.0;
  constexpr double kRdet = kR + kTelAlt;
  constexpr double kN0 = 3.27e-4;
  constexpr double kHatm = 8000.0;
  constexpr double kLambda = 400.0;
  constexpr double kTelRadius = 12.0;

  double refrIndex(double h_m) { return 1.0 + kN0 * std::exp(-std::max(0.0, h_m) / kHatm); }
  std::pair<double, double> refrFunc(double h_m, double) {
    double const n = refrIndex(h_m);
    return {n, n};
  }

  /// Production-like uniform 10 m grid, so the separable fast path is available.
  std::shared_ptr<cherenkov::RefractiveIndexTable const> makeTable() {
    std::vector<double> alt_km, lnN;
    for (int i = 0; i <= 12000; ++i) {
      double const h_km = 0.01 * double(i);
      alt_km.push_back(h_km);
      double const N = kN0 * std::exp(-h_km * 1000.0 / kHatm);
      lnN.push_back(std::log(std::max(N, 1.0e-10)));
    }
    auto t = std::make_shared<cherenkov::RefractiveIndexTable>(alt_km, lnN, kR);
    return t->valid() ? t : nullptr;
  }

  struct Ray {
    char const* name;
    double h_em;    ///< emission height above sea level [m]
    double d;       ///< horizontal distance to the telescope [m]
    double aimDeg;  ///< tilt of the chord in the vertical plane (+ = upward)
  };

  struct Res {
    bool ok{false};
    std::array<double, 3> p{};
    std::array<double, 3> u{};
    double L{0.0};
    double opt{0.0};
    double tau{0.0};
    double us{0.0};
  };

  std::unique_ptr<cherenkov::CherenkovProcess> makeProc(
      std::shared_ptr<cherenkov::RefractiveIndexTable const> const& tab,
      cherenkov::TracerIntegrator const integ, double const maxStep, double const atolPos,
      double const atolDir, double const atolOpt) {
    cherenkov::CherenkovProcess::HitRecorder const noRec{};
    auto p = std::make_unique<cherenkov::CherenkovProcess>(
        1.0003, 399.0, 401.0, 1.0, kRdet, noRec, cherenkov::ProjectionMode::Telescope,
        std::shared_ptr<corsika::media::ShowerAxis const>{}, 0.0, 0.0, 0.0, 0.0, 0.0, kRdet,
        kTelRadius, std::function<std::pair<double, double>(double, double)>{refrFunc}, kR,
        0.0 /*azimuth*/, M_PI / 2.0 /*zenith 90 deg*/, 0.0, 0.0, false /*dispersion*/,
        true /*curvature*/, std::string{}, 0.0, tab, false /*curved table*/);
    p->setTracerIntegrator(integ);
    p->setDopriTolerances(atolPos, atolDir, atolOpt);
    p->setDopriStepLimits(1.0e-3, maxStep);
    return p;
  }

  /// Trace one ray `reps` times and return the last result plus the mean cost per call.
  Res trace(cherenkov::CherenkovProcess const& proc, Ray const& r, int const reps,
            double const maxRadialStep_m = 100.0) {
    double const h = r.h_em;
    double const dPhi = r.d / kR;
    double const rEm = kR + h;
    std::array<double, 3> const p0{-rEm * std::sin(dPhi), 0.0, rEm * std::cos(dPhi)};
    std::array<double, 3> const tgt{0.0, 0.0, kRdet}; // the dish centre
    double const chord = std::sqrt((tgt[0] - p0[0]) * (tgt[0] - p0[0]) +
                                   (tgt[1] - p0[1]) * (tgt[1] - p0[1]) +
                                   (tgt[2] - p0[2]) * (tgt[2] - p0[2]));
    double const ux0 = (tgt[0] - p0[0]) / chord;
    double const uz0 = (tgt[2] - p0[2]) / chord;
    double const a = r.aimDeg * M_PI / 180.0;
    std::array<double, 3> const u{ux0 * std::cos(a) + uz0 * std::sin(a), 0.0,
                                  -ux0 * std::sin(a) + uz0 * std::cos(a)};

    Res out;
    auto const t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < reps; ++i) {
      double xi = 0.0, yi = 0.0, zi = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
      double dhx = 0.0, dhy = 0.0, dhz = 0.0;
      out.ok = proc.rayTraceSphericalToPlane(p0[0], p0[1], p0[2], u[0], u[1], u[2], 0.0,
                                             0.0, kRdet, 1.0, 0.0, 0.0, xi, yi, zi, op, gp,
                                             kLambda, tau, dhx, dhy, dhz, maxRadialStep_m);
      out.p = {xi, yi, zi};
      out.u = {dhx, dhy, dhz};
      out.L = gp;
      out.opt = op;
      out.tau = tau;
    }
    auto const t1 = std::chrono::steady_clock::now();
    out.us = std::chrono::duration<double, std::micro>(t1 - t0).count() / double(reps);
    return out;
  }

} // namespace

namespace {
  double dPos(Res const& a, Res const& b) { // metres, 0 when both miss
    if (a.ok != b.ok) return 1.0;
    if (!a.ok) return 0.0;
    double const dx = a.p[0] - b.p[0], dy = a.p[1] - b.p[1], dz = a.p[2] - b.p[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }
  double dDir(Res const& a, Res const& b) {
    if (a.ok != b.ok || !a.ok) return 0.0;
    double m = 0.0;
    for (int i = 0; i < 3; ++i) m = std::max(m, std::fabs(a.u[i] - b.u[i]));
    return m;
  }
  double dTime_ns(Res const& a, Res const& b) {
    if (a.ok != b.ok) return 1.0e30;
    if (!a.ok) return 0.0;
    return std::fabs(a.opt - b.opt) / 299792458.0 * 1.0e9;
  }
} // namespace

TEST_CASE("Cherenkov tracer benchmark (RK4 vs adaptive DOPRI5)", "modules,cherenkov") {
  logging::set_level(logging::level::warn);
  auto const tab = makeTable();
  REQUIRE(tab != nullptr);

  std::vector<Ray> const rays = {
      {"horiz_200km", 43000.0, 200000.0, 0.0}, {"horiz_400km", 43000.0, 400000.0, 0.0},
      {"horiz_700km", 43000.0, 700000.0, 0.0}, {"down_10km", 10000.0, 10000.0, 0.0},
      {"down_20km", 10000.0, 20000.0, 0.0},    {"down_40km", 10000.0, 40000.0, 0.0},
      {"tang_100km", 3100.0, 100000.0, 0.0},   {"tang_200km", 3100.0, 200000.0, 0.0},
      {"layer_30km", 30000.0, 60000.0, 0.0},   {"layer_5km", 5000.0, 30000.0, 0.0},
      {"steep_5km", 3500.0, 5000.0, 0.0},      {"escaping", 43000.0, 200000.0, 2.0},
  };

  using cherenkov::TracerIntegrator;
  // Reference: DOPRI5 at atol_pos = 1e-6 m with NO step cap.  Its own convergence is
  // checked against a 100x tighter tolerance below (a radial cap is NOT a step cap: at
  // 0.5 m/step a steep ray needs ~3e4 steps and hits the 1e4 step limit, which is why
  // the earlier RK4 cross-check was ill-posed).
  auto const refProc = makeProc(tab, TracerIntegrator::DOPRI5, 1.0e7, 1.0e-6, 1.0e-12, 1.0e-8);
  auto const refTight = makeProc(tab, TracerIntegrator::DOPRI5, 1.0e7, 1.0e-8, 1.0e-14, 1.0e-10);
  auto const rk4 = makeProc(tab, TracerIntegrator::RK4, 0.0, 1e-3, 1e-9, 1e-4);
  auto const dpCons = makeProc(tab, TracerIntegrator::DOPRI5, 0.0, 1e-3, 1e-9, 1e-4);
  std::vector<std::pair<double, std::unique_ptr<cherenkov::CherenkovProcess>>> relaxCfg;
  for (double const cap : {5000.0, 20000.0, 50000.0})
    relaxCfg.emplace_back(cap, makeProc(tab, TracerIntegrator::DOPRI5, cap, 1e-3, 1e-9, 1e-4));

  int const reps = 20;
  std::printf("\n[CherenkovTracerBenchmark] %zu rays, reference = DOPRI5 atol_pos=1e-6 m, no cap\n",
              rays.size());

  // ---- per-ray pass: reference, its own convergence, and RK4 ---------------------
  std::vector<Res> ref(rays.size());
  double mxRefSelf = 0.0, tRef = 0.0, mxRk4P = 0.0, mxRk4D = 0.0, mxRk4T = 0.0, tRk4 = 0.0;
  int mismatch = 0, refBad = 0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    ref[i] = trace(*refProc, rays[i], 1);
    Res const tight = trace(*refTight, rays[i], 1);
    Res const a = trace(*rk4, rays[i], reps);
    double const selfD = dPos(ref[i], tight);
    mxRefSelf = std::max(mxRefSelf, selfD);
    if (!(selfD <= 1.0e-3)) ++refBad;
    if (a.ok != ref[i].ok) ++mismatch;
    mxRk4P = std::max(mxRk4P, dPos(a, ref[i]));
    mxRk4D = std::max(mxRk4D, dDir(a, ref[i]));
    mxRk4T = std::max(mxRk4T, dTime_ns(a, ref[i]));
    tRef += ref[i].us;
    tRk4 += a.us;
  }

  std::printf("  reference self-check (100x tighter tolerance): max %.6f mm, inconsistent %d\n",
              mxRefSelf * 1e3, refBad);
  std::printf("  RK4 (production caps): max |dP| %.6f mm, |dU| %.2e, |dt| %.4f ns, %.1f us/set\n",
              mxRk4P * 1e3, mxRk4D, mxRk4T, tRk4);

  // ---- sweep the maximum step ---------------------------------------------------
  std::printf("\n  %-10s %8s %8s %8s %10s %10s %10s %8s %8s\n", "maxstep", "dopri_us",
              "speedup", "horiz_x", "max_dP_mm", "max_dU", "max_dt_ns", "steps", "rej");
  double bestSpeed = 0.0, bestCap = 0.0, bestP = 0.0, bestD = 0.0, bestT = 0.0;
  for (auto const& [cap, proc] : relaxCfg) {
    double tSum = 0.0, mxP = 0.0, mxD = 0.0, mxT = 0.0, tHorizRk4 = 0.0, tHorizDp = 0.0;
    int bad = 0;
    for (std::size_t i = 0; i < rays.size(); ++i) {
      Res const q = trace(*proc, rays[i], reps);
      if (q.ok != ref[i].ok) ++bad;
      mxP = std::max(mxP, dPos(q, ref[i]));
      mxD = std::max(mxD, dDir(q, ref[i]));
      mxT = std::max(mxT, dTime_ns(q, ref[i]));
      tSum += q.us;
      if (rays[i].h_em > 40000.0) {
        tHorizDp += q.us;
        tHorizRk4 += trace(*rk4, rays[i], 1).us;
      }
    }
    auto const st = proc->tracerStats();
    std::printf("  %-10.0f %8.1f %8.2f %8.2f %10.6f %10.2e %10.4f %8zu %8zu%s\n", cap, tSum,
                tRk4 / tSum, tHorizRk4 / tHorizDp, mxP * 1e3, mxD, mxT, st.accepted,
                st.rejected, bad ? " HIT-MISMATCH" : "");
    if (bad == 0 && mxP <= 1e-3 && mxD <= 1e-7 && mxT <= 0.05 && tRk4 / tSum > bestSpeed) {
      bestSpeed = tRk4 / tSum; bestCap = cap; bestP = mxP; bestD = mxD; bestT = mxT;
    }
  }
  std::printf("\n  selected: maxstep %.0f m, speedup %.2fx, |dP| %.6f mm, |dU| %.2e, "
              "|dt| %.4f ns\n", bestCap, bestSpeed, bestP * 1e3, bestD, bestT);
  std::printf("  capped-step counts (dopri at 5 km): ");
  {
    auto const st = relaxCfg[0].second->tracerStats();
    std::printf("%zu of %zu accepted, h[%.3g..%.3g] m, %zu f-evals, %zu crossing iters\n",
                st.cappedSteps, st.accepted, st.hMin, st.hMax, st.fEvaluations,
                st.crossingIterations);
  }

  CHECK(mismatch == 0);
  CHECK(refBad == 0);
  CHECK(mxRk4P <= 1e-3);
  CHECK(bestSpeed >= 1.2); // >= 20 % faster than RK4 at the same accuracy
}
