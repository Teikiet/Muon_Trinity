/*
 * testCherenkovImpactController.cpp
 *
 * Phase 4, step 3: the impact-aware accuracy controller.  Answers three questions with
 * measurements rather than assertions of intent:
 *
 *   1. Does it hold the HIT error under 1 mm, including for shallow rays?  (vs an
 *      independent, self-checked reference)
 *   2. What does it cost (time, steps, evaluations)?
 *   3. Is the extra crossing correction (land the crossing step on the plane) worth it
 *      once the controller is in place?
 *
 * Why a new test and not an extension of the benchmark: the benchmark runs one plane
 * (telescope-normal) and measures |dP| against a reference with the SAME interpolated
 * crossing.  Here the plane normal is swept, because sin(theta_arr) is what amplifies a
 * position error into the hit, and the reference lands its crossing so it cannot hide the
 * defect being measured.
 *
 * Reference: atol_pos 1e-8 m, atol_dir 1e-14, atol_opt 1e-10, no step cap, crossing
 * re-landed.  Its own convergence is checked against a 100x tighter setting and asserted,
 * which is the Phase-1 lesson: never gate against an unconverged reference.
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
  const double kPi = 3.14159265358979323846;

  std::pair<double, double> refrFunc(double h_m, double) {
    double const n = 1.0 + kN0 * std::exp(-std::max(0.0, h_m) / kHatm);
    return {n, n};
  }

  std::shared_ptr<cherenkov::RefractiveIndexTable const> makeTable() {
    std::vector<double> alt_km, lnN;
    for (int i = 0; i <= 12000; ++i) {
      double const h_km = 0.01 * double(i);
      alt_km.push_back(h_km);
      lnN.push_back(std::log(std::max(kN0 * std::exp(-h_km * 1000.0 / kHatm), 1.0e-10)));
    }
    auto t = std::make_shared<cherenkov::RefractiveIndexTable>(alt_km, lnN, kR);
    return t->valid() ? t : nullptr;
  }

  struct Ray {
    char const* name;
    double h_em, d, aimDeg;
  };

  /// One ray + one target plane + one configuration.
  struct Hit {
    bool ok = false;
    double x = 0, y = 0, z = 0;
    double len = 0;
  };

  void traceOne(cherenkov::CherenkovProcess const& proc, Ray const& r,
                std::array<double, 3> const& nrm, double const maxRadialStep_m, Hit& out) {
    double const dPhi = r.d / kR, rEm = kR + r.h_em;
    std::array<double, 3> const p0{-rEm * std::sin(dPhi), 0.0, rEm * std::cos(dPhi)};
    std::array<double, 3> const tgt{0.0, 0.0, kRdet};
    double const chord = std::sqrt((tgt[0] - p0[0]) * (tgt[0] - p0[0]) +
                                   (tgt[2] - p0[2]) * (tgt[2] - p0[2]));
    double const ux0 = (tgt[0] - p0[0]) / chord, uz0 = (tgt[2] - p0[2]) / chord;
    double const a = r.aimDeg * kPi / 180.0;
    std::array<double, 3> const u{ux0 * std::cos(a) + uz0 * std::sin(a), 0.0,
                                  -ux0 * std::sin(a) + uz0 * std::cos(a)};
    double xi = 0, yi = 0, zi = 0, op = 0, gp = 0, tau = 0, dhx = 0, dhy = 0, dhz = 0;
    out.ok = proc.rayTracePlaneDopri(p0[0], p0[1], p0[2], u[0], u[1], u[2], 0.0, 0.0,
                                     kRdet, nrm[0], nrm[1], nrm[2], xi, yi, zi, op, gp,
                                     kLambda, tau, dhx, dhy, dhz, maxRadialStep_m);
    out.x = xi;
    out.y = yi;
    out.z = zi;
    out.len = gp;
  }

  double dPos(Hit const& a, Hit const& b) {
    if (a.ok != b.ok) return 1.0e30;
    if (!a.ok) return 0.0;
    double const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  std::unique_ptr<cherenkov::CherenkovProcess> makeProc(
      std::shared_ptr<cherenkov::RefractiveIndexTable const> const& tab, double const atolPos,
      double const atolDir, double const atolOpt, double const maxStep,
      double const impactBudget, bool const relanding) {
    cherenkov::CherenkovProcess::HitRecorder const noRec{};
    auto p = std::make_unique<cherenkov::CherenkovProcess>(
        1.0003, 399.0, 401.0, 1.0, kRdet, noRec, cherenkov::ProjectionMode::Telescope,
        std::shared_ptr<corsika::media::ShowerAxis const>{}, 0.0, 0.0, 0.0, 0.0, 0.0, kRdet,
        12.0, std::function<std::pair<double, double>(double, double)>{refrFunc}, kR, 0.0,
        kPi / 2.0, 0.0, 0.0, false, true, std::string{}, 0.0, tab, false);
    p->setTracerIntegrator(cherenkov::TracerIntegrator::DOPRI5);
    p->setDopriTolerances(atolPos, atolDir, atolOpt);
    p->setDopriStepLimits(1.0e-3, maxStep);
    p->setImpactBudget(impactBudget);
    p->setCrossingRelanding(relanding);
    return p;
  }
} // namespace


TEST_CASE("Cherenkov impact-aware controller: hit error under 1 mm", "modules,cherenkov") {
  logging::set_level(logging::level::warn);
  auto const tab = makeTable();
  REQUIRE(tab != nullptr);

  // A stress set: long horizontals, mid paths, near-tangent rays at the dish altitude and
  // a steep one.  Nothing here is tuned to the controller.
  std::vector<Ray> const rays = {
      {"horiz_200km", 43000.0, 200000.0, 0.0}, {"horiz_400km", 43000.0, 400000.0, 0.0},
      {"down_20km", 10000.0, 20000.0, 0.0},    {"down_40km", 10000.0, 40000.0, 0.0},
      {"tang_100km", 3100.0, 100000.0, 0.0},   {"tang_200km", 3100.0, 200000.0, 0.0},
      {"layer_30km", 30000.0, 60000.0, 0.0},   {"layer_5km", 5000.0, 30000.0, 0.0},
      {"steep_5km", 3500.0, 5000.0, 0.0},
  };
  // Plane orientations: the normal is tilted off the aim direction, so sin(theta_arr) runs
  // from ~0 (head-on) toward ~1 (the ray runs in the plane = the shallow case).
  std::vector<double> const thetas = {0.0, 15.0, 45.0, 75.0, 80.0, 85.0};
  std::vector<std::array<double, 3>> normals;
  for (double const th : thetas) {
    double const a = th * kPi / 180.0;
    normals.push_back({std::cos(a), 0.0, std::sin(a)});
  }
  std::size_t const nTh = normals.size(), nR = rays.size(), nCase = nR * nTh;

  // Reference: re-landed crossing (so it cannot hide the interpolant defect) at 1e-8, and a
  // 100x tighter copy for the self-check.
  auto const refP = makeProc(tab, 1.0e-8, 1.0e-14, 1.0e-10, 1.0e7, 0.0, true);
  auto const refT = makeProc(tab, 1.0e-10, 1.0e-16, 1.0e-12, 1.0e7, 0.0, true);
  // Production-like configurations, differing ONLY in the two new switches.
  auto const cA = makeProc(tab, 1.0e-3, 1.0e-9, 1.0e-4, 20000.0, 0.0, false);    // shipped
  auto const cB = makeProc(tab, 1.0e-3, 1.0e-9, 1.0e-4, 20000.0, 5.0e-4, false); // +controller
  auto const cC = makeProc(tab, 1.0e-3, 1.0e-9, 1.0e-4, 20000.0, 0.0, true);     // +landing
  auto const cD = makeProc(tab, 1.0e-3, 1.0e-9, 1.0e-4, 20000.0, 5.0e-4, true);  // +both

  std::vector<double> eA(nCase), eB(nCase), eC(nCase), eD(nCase);
  double refSelf = 0.0;
  int missA = 0, missB = 0, missC = 0, missD = 0;

  std::printf("\n[CherenkovImpactController] %zu rays x %zu plane orientations\n", nR, nTh);
  std::printf("  %-13s %6s %11s %11s %11s %11s\n", "ray", "theta", "shipped", "impact",
              "landing", "both");
  for (std::size_t i = 0; i < nR; ++i) {
    for (std::size_t j = 0; j < nTh; ++j) {
      std::size_t const c = i * nTh + j;
      Hit ref, tight, a, b, cc, d;
      traceOne(*refP, rays[i], normals[j], 0.0, ref);
      traceOne(*refT, rays[i], normals[j], 0.0, tight);
      refSelf = std::max(refSelf, dPos(ref, tight));
      traceOne(*cA, rays[i], normals[j], 0.0, a);
      traceOne(*cB, rays[i], normals[j], 0.0, b);
      traceOne(*cC, rays[i], normals[j], 0.0, cc);
      traceOne(*cD, rays[i], normals[j], 0.0, d);
      eA[c] = dPos(a, ref);
      eB[c] = dPos(b, ref);
      eC[c] = dPos(cc, ref);
      eD[c] = dPos(d, ref);
      if (eA[c] > 1e20) ++missA;
      if (eB[c] > 1e20) ++missB;
      if (eC[c] > 1e20) ++missC;
      if (eD[c] > 1e20) ++missD;
      std::printf("  %-13s %6.0f %11.6f %11.6f %11.6f %11.6f\n", rays[i].name, thetas[j],
                  eA[c] * 1e3, eB[c] * 1e3, eC[c] * 1e3, eD[c] * 1e3);
    }
  }
  auto const mx = [](std::vector<double> const& v) {
    double m = 0.0;
    for (double const x : v)
      if (x < 1e20) m = std::max(m, x);
    return m;
  };
  std::printf("  reference self-check (1e-10 vs 1e-8): max %.6f mm\n", refSelf * 1e3);
  std::printf("\n  %-8s %12s %12s %12s %12s\n", "theta", "shipped", "impact", "landing",
              "both");
  for (std::size_t j = 0; j < nTh; ++j) {
    double a = 0, b = 0, cc = 0, d = 0;
    for (std::size_t i = 0; i < nR; ++i) {
      std::size_t const k = i * nTh + j;
      if (eA[k] < 1e20) a = std::max(a, eA[k]);
      if (eB[k] < 1e20) b = std::max(b, eB[k]);
      if (eC[k] < 1e20) cc = std::max(cc, eC[k]);
      if (eD[k] < 1e20) d = std::max(d, eD[k]);
    }
    std::printf("  %-8.0f %12.6f %12.6f %12.6f %12.6f\n", thetas[j], a * 1e3, b * 1e3,
                cc * 1e3, d * 1e3);
  }
  std::printf("  max over all    shipped %.6f | impact %.6f | landing %.6f | both %.6f mm\n",
              mx(eA) * 1e3, mx(eB) * 1e3, mx(eC) * 1e3, mx(eD) * 1e3);
  std::printf("  miss mismatches %d / %d / %d / %d (shipped/impact/landing/both)\n", missA,
              missB, missC, missD);


  // ---- cost: per-case wall clock (best of 3) and the step statistics -------------
  auto const timeCfg = [&](cherenkov::CherenkovProcess const& p, int const reps) {
    double best = 1.0e30;
    for (int rep = 0; rep < reps; ++rep) {
      auto const t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; i < nR; ++i)
        for (std::size_t j = 0; j < nTh; ++j) {
          Hit h;
          traceOne(p, rays[i], normals[j], 0.0, h);
        }
      auto const t1 = std::chrono::steady_clock::now();
      double const us =
          std::chrono::duration<double, std::micro>(t1 - t0).count() / double(nCase);
      best = std::min(best, us);
    }
    return best;
  };
  cA->resetTracerStats();
  cB->resetTracerStats();
  cC->resetTracerStats();
  cD->resetTracerStats();
  double const tA = timeCfg(*cA, 3);
  double const tB = timeCfg(*cB, 3);
  double const tC = timeCfg(*cC, 3);
  double const tD = timeCfg(*cD, 3);
  auto const stA = cA->tracerStats();
  auto const stB = cB->tracerStats();
  auto const stC = cC->tracerStats();
  auto const stD = cD->tracerStats();
  std::printf("\n  cost per case   shipped %.1f us | impact %.1f us (%.2fx) | landing %.1f us "
              "(%.2fx) | both %.1f us (%.2fx)\n",
              tA, tB, tB / tA, tC, tC / tA, tD, tD / tA);
  std::printf("  steps per case  shipped %.1f | impact %.1f | landing %.1f | both %.1f\n",
              double(stA.accepted) / double(nCase), double(stB.accepted) / double(nCase),
              double(stC.accepted) / double(nCase), double(stD.accepted) / double(nCase));
  std::printf("  impact control  accepted %zu, shortened %zu, uncertified %zu | both: accepted "
              "%zu, shortened %zu, uncertified %zu\n",
              stB.impactAccepted, stB.impactRejections, cB->impactUncertified(),
              stD.impactAccepted, stD.impactRejections, cD->impactUncertified());

  // ---- the claims ---------------------------------------------------------------
  double shallowA = 0.0, shallowB = 0.0, shallowC = 0.0, shallowD = 0.0;
  for (std::size_t i = 0; i < nR; ++i) {
    std::size_t const k = i * nTh + (nTh - 1); // 75 deg = the shallowest plane
    if (eA[k] < 1e20) shallowA = std::max(shallowA, eA[k]);
    if (eB[k] < 1e20) shallowB = std::max(shallowB, eB[k]);
    if (eC[k] < 1e20) shallowC = std::max(shallowC, eC[k]);
    if (eD[k] < 1e20) shallowD = std::max(shallowD, eD[k]);
  }
  std::printf("  shallow rows (75 deg): shipped %.6f | impact %.6f | landing %.6f | both "
              "%.6f mm\n",
              shallowA * 1e3, shallowB * 1e3, shallowC * 1e3, shallowD * 1e3);

  CHECK(refSelf <= 5.0e-5);          // the reference is converged to better than 0.05 mm
  CHECK(missB == 0);                 // the controller never changes hit/miss
  CHECK(missD == 0);
  CHECK(mx(eB) <= 1.0e-3);           // THE claim: 1 mm with the controller alone
  CHECK(mx(eD) <= 1.0e-3);
  CHECK(shallowB <= 1.0e-3);         // ... including the shallowest plane
  CHECK(shallowD <= 1.0e-3);
  CHECK(mx(eB) <= mx(eA) + 1.0e-9);  // never worse than the shipped rule
  REQUIRE(stB.impactAccepted > 0);   // the controller actually engaged
  REQUIRE(stB.impactRejections > 0); // ... and it actually shortened steps
}
