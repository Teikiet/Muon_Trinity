/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * \file testCherenkovPlaneTracerAcceptance.cpp
 *
 * Step 2l (2g-h2l): frozen acceptance set for the Telescope-mode curved tracer,
 * CherenkovProcess::rayTraceSphericalToPlane.
 *
 * Until this set existed the plane tracer had no test at all: no test in the tree
 * constructs a CherenkovProcess in ProjectionMode::Telescope.
 *
 * GEOMETRY (frozen): telescope at R + 3000 m, pointing exactly at zenith 90 deg,
 *   azimuth 0 -> the telescope plane is VERTICAL, normal = +x, in-plane axes
 *   X = +y (horizontal transverse) and Y = +z (vertical).  A horizon photon
 *   crosses the plane face-on (no tangency), but any refraction-induced ALTITUDE
 *   shift lands directly in the dish coordinate yg -- so impact position, not
 *   only transmission, is the IACT criterion.
 *
 * TRUTH (frozen, independent of the tracer): RK4 in the spherically symmetric
 *   medium, ds = 1 m, n(h) = 1 + 3.27e-4 exp(-h/8 km), stop at the sign change
 *   of the signed distance to the telescope plane (that distance is linear along
 *   a step, so the interpolation is exact), tau exact in altitude against the
 *   MODTRAN table.  Validated here by the analytic n = 1 probe and by a ds/2
 *   convergence check.
 *
 * METRIC: eps = |T_prod - T_ref| / T_ref, <= 1 % pass, 1-5 % flag, > 5 % fail;
 *   |dT| <= 1e-3 when T_ref < 0.01.  Geometry is gated separately (relative path
 *   length, impact position, arrival direction) because it feeds timing; since the
 *   Phase-1 fixes (stage-local kappa, Hermite crossing) the geometry gates are
 *   |dL|/L <= 1e-9 and impact <= 1 mm, against measured 1.7e-12 and 6.2e-07 m.
 *   CHER_TRUTH_STEP_M=<m> refines the truth so its own convergence can be shown;
 *   the default (1 m) keeps every archived number reproducible.
 */

#include <corsika/modules/cherenkov/Cherenkov.hpp>
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <catch2/catch_all.hpp>

#include <PhysicalUnitsCatch2.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace corsika;

namespace {

  constexpr double kEarthRadius_m = 6.371e6;
  constexpr double kTelescopeAltitude_m = 3000.0;
  constexpr double kRdet = kEarthRadius_m + kTelescopeAltitude_m;
  // Truth step.  Settable from the environment so the FROZEN truth's own convergence
  // can be demonstrated (CHER_TRUTH_STEP_M=0.25): an acceptance metric must not be
  // reference-limited.  The default keeps every archived number reproducible.
  double const kStepRef_m = std::getenv("CHER_TRUTH_STEP_M")
                                ? std::atof(std::getenv("CHER_TRUTH_STEP_M"))
                                : 1.0;
  constexpr double kRefrN0 = 3.27e-4;
  constexpr double kRefrH_m = 8000.0;
  constexpr double kLambda_nm = 400.0;
  constexpr double kTelescopeRadius_m = 12.0;   // primary dish radius
  constexpr double kTelescopeRadiusWide_m = 30.0;
  constexpr double kEpsPass = 0.01;
  constexpr double kEpsAccept = 0.05;
  constexpr double kFovHalf_deg = 5.0;

  constexpr double kEmissions_km[] = {0.5, 2.0, 2.9, 3.0, 3.1, 5.0, 10.0, 20.0};
  constexpr double kDistances_km[] = {1.0, 10.0, 50.0, 100.0, 200.0, 400.0};

  double refrIndex(double h_m) { return 1.0 + kRefrN0 * std::exp(-std::max(0.0, h_m) / kRefrH_m); }
  double refrDnDr(double h_m) {
    return -kRefrN0 / kRefrH_m * std::exp(-std::max(0.0, h_m) / kRefrH_m);
  }
  std::pair<double, double> refrFunc(double h_m, double) {
    double const n = refrIndex(h_m);
    return {n, n};
  }

  // ------------------------------------------------------------------ table
  struct Table {
    std::size_t nLambda{0}, nHeight{0};
    std::vector<double> lambda_nm, height_m;
    std::vector<std::vector<double>> tau;
    std::string path;
    bool loaded() const { return nLambda > 0 && nHeight > 0; }
    std::size_t rowOf(double lambdaNm) const {
      std::size_t best = 0;
      double bestD = 1.0e30;
      for (std::size_t i = 0; i < lambda_nm.size(); ++i) {
        double const d = std::fabs(lambda_nm[i] - lambdaNm);
        if (d < bestD) { bestD = d; best = i; }
      }
      return best;
    }
    std::size_t slab(double h_m) const {
      std::size_t ih = 0;
      while (ih + 1 < height_m.size() && height_m[ih + 1] <= h_m) ++ih;
      return ih;
    }
    double tauAt(double h_m, std::size_t row) const {
      std::size_t const ih = slab(h_m);
      if (ih + 1 >= height_m.size()) return tau[row][ih];
      double const dh = height_m[ih + 1] - height_m[ih];
      double const f = (h_m - height_m[ih]) / dh;
      return tau[row][ih] * (1.0 - f) + tau[row][ih + 1] * f;
    }
    double alpha(double h_m, std::size_t row) const {
      std::size_t const ih = slab(h_m);
      if (ih + 1 >= height_m.size()) return 0.0;
      return (tau[row][ih + 1] - tau[row][ih]) / (height_m[ih + 1] - height_m[ih]);
    }
  };

  Table loadTable() {
    std::vector<std::string> const candidates = {
        "/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/"
        "abstable_MODTRAN_bottomup.dat",
        "modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat",
        "../modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat",
        "../../modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat"};
    for (auto const& path : candidates) {
      std::ifstream in(path);
      if (!in.is_open()) continue;
      std::size_t nl = 0, nh = 0;
      if (!(in >> nl >> nh) || nl < 45 || nh == 0) continue;
      Table t;
      t.nLambda = nl;
      t.nHeight = nh;
      t.lambda_nm.resize(nl);
      t.height_m.resize(nh);
      bool ok = true;
      for (std::size_t i = 0; i < nl && ok; ++i) ok = static_cast<bool>(in >> t.lambda_nm[i]);
      for (std::size_t i = 0; i < nh && ok; ++i) ok = static_cast<bool>(in >> t.height_m[i]);
      t.tau.assign(nl, std::vector<double>(nh, 0.0));
      for (std::size_t i = 0; i < nl && ok; ++i)
        for (std::size_t j = 0; j < nh && ok; ++j) ok = static_cast<bool>(in >> t.tau[i][j]);
      if (!ok) continue;
      t.path = path;
      return t;
    }
    return Table{};
  }

}  // namespace

namespace {

  // ------------------------------------------------------- frozen truth
  struct Truth {
    bool ok{false};
    char const* why{""};
    double L{0.0};
    double optical{0.0};   // int n_group ds (n_group == n_phase here)
    double tau{0.0};
    std::array<double, 3> p{0.0, 0.0, 0.0};
    std::array<double, 3> u{0.0, 0.0, 1.0};
    double hMin{0.0};
  };

  /// RK4, ds = kStepRef_m, stop at the sign change of n.(p - c).
  Truth traceTruthPlane(Table const& tab, std::size_t const row, std::array<double, 3> const& p0,
                        std::array<double, 3> const& u0, std::array<double, 3> const& c,
                        std::array<double, 3> const& nrm) {
    Truth out;
    out.p = p0;
    out.u = u0;
    out.hMin = std::sqrt(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]) - kEarthRadius_m;
    auto dist = [&](std::array<double, 3> const& q) {
      return (q[0] - c[0]) * nrm[0] + (q[1] - c[1]) * nrm[1] + (q[2] - c[2]) * nrm[2];
    };
    auto deriv = [](std::array<double, 3> const& q, std::array<double, 3> const& v,
                    std::array<double, 3>& dq, std::array<double, 3>& dv) {
      double const r = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
      double const k = refrDnDr(r - kEarthRadius_m) / refrIndex(r - kEarthRadius_m);
      double const rx = q[0] / r, ry = q[1] / r, rz = q[2] / r;
      double const ur = v[0] * rx + v[1] * ry + v[2] * rz;
      dq = v;
      dv = {k * (rx - ur * v[0]), k * (ry - ur * v[1]), k * (rz - ur * v[2])};
    };
    std::array<double, 3> q = p0, v = u0;
    double s = 0.0;
    double dPrev = dist(q);
    double hPrev = out.hMin;
    for (std::size_t it = 0; it < 3000000ULL; ++it) {
      std::array<double, 3> a1, b1, a2, b2, a3, b3, a4, b4, tq, tv;
      deriv(q, v, a1, b1);
      for (int i = 0; i < 3; ++i) { tq[i] = q[i] + 0.5 * kStepRef_m * a1[i]; tv[i] = v[i] + 0.5 * kStepRef_m * b1[i]; }
      deriv(tq, tv, a2, b2);
      for (int i = 0; i < 3; ++i) { tq[i] = q[i] + 0.5 * kStepRef_m * a2[i]; tv[i] = v[i] + 0.5 * kStepRef_m * b2[i]; }
      deriv(tq, tv, a3, b3);
      for (int i = 0; i < 3; ++i) { tq[i] = q[i] + kStepRef_m * a3[i]; tv[i] = v[i] + kStepRef_m * b3[i]; }
      deriv(tq, tv, a4, b4);
      for (int i = 0; i < 3; ++i) {
        q[i] += kStepRef_m / 6.0 * (a1[i] + 2.0 * a2[i] + 2.0 * a3[i] + a4[i]);
        v[i] += kStepRef_m / 6.0 * (b1[i] + 2.0 * b2[i] + 2.0 * b3[i] + b4[i]);
      }
      double const vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      for (int i = 0; i < 3; ++i) v[i] /= vn;
      double const r = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
      double const h = r - kEarthRadius_m;
      out.hMin = std::min(out.hMin, h);
      if (h < -100.0) { out.why = "below sea level"; return out; }
      double const dh = h - hPrev;
      double const dt = (std::fabs(dh) > 1.0e-6)
                            ? std::fabs(tab.tauAt(h, row) - tab.tauAt(hPrev, row)) *
                                  kStepRef_m / std::fabs(dh)
                            : tab.alpha(0.5 * (h + hPrev), row) * kStepRef_m;
      out.tau += dt;
      out.optical += refrIndex(0.5 * (h + hPrev)) * kStepRef_m;
      s += kStepRef_m;
      double const dCur = dist(q);
      if (dPrev * dCur <= 0.0) {
        double const alpha = std::fabs(dPrev) / std::max(1e-300, std::fabs(dPrev) + std::fabs(dCur));
        out.ok = true;
        out.L = s - kStepRef_m + alpha * kStepRef_m;
        out.tau -= (1.0 - alpha) * dt;
        out.optical -= (1.0 - alpha) * refrIndex(0.5 * (h + hPrev)) * kStepRef_m;
        for (int i = 0; i < 3; ++i) out.p[i] = q[i] - (1.0 - alpha) * kStepRef_m * v[i];
        out.u = v;
        return out;
      }
      dPrev = dCur;
      hPrev = h;
      if (h > 100000.0) { out.why = "escaped"; return out; }
    }
    out.why = "step limit";
    return out;
  }

}  // namespace

TEST_CASE("Cherenkov plane tracer acceptance (telescope mode, 90 deg pointing)",
          "modules,cherenkov") {
  logging::set_level(logging::level::warn);
  Table const tab = loadTable();
  REQUIRE(tab.loaded());
  std::size_t const row = tab.rowOf(kLambda_nm);

  std::printf("\n[CherenkovPlaneTracerAcceptance] telescope at %.0f m, pointing zenith 90 deg, "
              "azimuth 0, dish r = %.0f m / %.0f m, lambda = %.0f nm\n",
              kTelescopeAltitude_m, kTelescopeRadius_m, kTelescopeRadiusWide_m, kLambda_nm);
  std::printf("  truth: RK4 ds = 1 m, plane stop, tau exact in altitude (no window)\n");
  std::printf("  telescope basis at zenith 90 deg: Z = +x (normal), X = +y, Y = +z (vertical)\n");

  // Telescope plane: point = telescope centre, normal = +x.  In-plane axes +y, +z.
  std::array<double, 3> const planePt{0.0, 0.0, kRdet};
  std::array<double, 3> const planeN{1.0, 0.0, 0.0};
  auto planeCoords = [&](std::array<double, 3> const& q, double& xg, double& yg) {
    xg = q[1] - planePt[1];
    yg = q[2] - planePt[2];
  };

  cherenkov::CherenkovProcess::HitRecorder const noRec{};
  // (a) absorption ON: needed because this tracer only accumulates tau when
  //     absorption is enabled; it then samples survival internally from a
  //     time-seeded thread-local RNG, so tau is taken from the first call that
  //     SURVIVES (the geometry and tau of a given ray are deterministic).
  cherenkov::CherenkovProcess const proc(
      1.0003, 399.0, 401.0, 1.0, kRdet, noRec, cherenkov::ProjectionMode::Telescope,
      std::shared_ptr<corsika::media::ShowerAxis const>{}, 0.0, 0.0, 0.0,
      planePt[0], planePt[1], planePt[2], kTelescopeRadius_m,
      std::function<std::pair<double, double>(double, double)>{refrFunc}, kEarthRadius_m,
      0.0 /*azimuth*/, M_PI / 2.0 /*zenith 90 deg*/, 0.0, 0.0, false /*dispersion*/,
      true /*curvature*/, tab.path, 0.0,
      std::shared_ptr<cherenkov::RefractiveIndexTable const>{}, false);
  // (b) absorption OFF: deterministic geometry probe for the large sweep sets.
  cherenkov::CherenkovProcess const geom(
      1.0003, 399.0, 401.0, 1.0, kRdet, noRec, cherenkov::ProjectionMode::Telescope,
      std::shared_ptr<corsika::media::ShowerAxis const>{}, 0.0, 0.0, 0.0,
      planePt[0], planePt[1], planePt[2], kTelescopeRadius_m,
      std::function<std::pair<double, double>(double, double)>{refrFunc}, kEarthRadius_m,
      0.0, M_PI / 2.0, 0.0, 0.0, false, true, std::string{}, 0.0,
      std::shared_ptr<cherenkov::RefractiveIndexTable const>{}, false);

  struct Row {
    double h{0.0}, d{0.0};
    double Lref{0.0}, Ltrc{0.0};
    double tauRef{0.0}, tauTrc{0.0};
    double eps{0.0};
    double xgRef{0.0}, ygRef{0.0}, xgTrc{0.0}, ygTrc{0.0};
    double offRef{0.0}, offTrc{0.0};
    double thetaArr{0.0};
    int calls{0};
    bool tauScored{false};
    bool ok{false};
    bool inDishRef{false}, inDishTrc{false};
    char const* status{""};
  };
  std::vector<Row> rows;

  double usSum = 0.0;
  int callCount = 0;

  for (double const hKm : kEmissions_km) {
    for (double const dKm : kDistances_km) {
      Row r;
      r.h = hKm * 1000.0;
      r.d = dKm * 1000.0;
      double const dPhi = r.d / kEarthRadius_m;
      double const rEm = kEarthRadius_m + r.h;
      std::array<double, 3> const p0{-rEm * std::sin(dPhi), 0.0, rEm * std::cos(dPhi)};
      double const chord = std::sqrt((planePt[0] - p0[0]) * (planePt[0] - p0[0]) +
                                     (planePt[1] - p0[1]) * (planePt[1] - p0[1]) +
                                     (planePt[2] - p0[2]) * (planePt[2] - p0[2]));
      std::array<double, 3> const u{(planePt[0] - p0[0]) / chord,
                                    (planePt[1] - p0[1]) / chord,
                                    (planePt[2] - p0[2]) / chord};
      Truth const tr = traceTruthPlane(tab, row, p0, u, planePt, planeN);
      if (!tr.ok) {
        r.status = tr.why;
        rows.push_back(r);
        continue;
      }
      r.Lref = tr.L;
      r.tauRef = tr.tau;
      planeCoords(tr.p, r.xgRef, r.ygRef);
      r.offRef = std::sqrt(r.xgRef * r.xgRef + r.ygRef * r.ygRef);
      r.thetaArr = std::acos(std::max(-1.0, std::min(1.0, -tr.u[2]))) * 180.0 / M_PI;
      r.inDishRef = r.offRef <= kTelescopeRadius_m;

      // Geometry + impact from the ABSORPTION-OFF instance: deterministic, one call.
      {
        double xi = 0.0, yi = 0.0, zi = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
        double dhx = 0.0, dhy = 0.0, dhz = 0.0;  // arrival direction output
        auto const t0 = std::chrono::steady_clock::now();
        bool const ok = geom.rayTraceSphericalToPlane(
            p0[0], p0[1], p0[2], u[0], u[1], u[2], planePt[0], planePt[1], planePt[2],
            planeN[0], planeN[1], planeN[2], xi, yi, zi, op, gp, kLambda_nm, tau,
            dhx, dhy, dhz);
        auto const t1 = std::chrono::steady_clock::now();
        usSum += std::chrono::duration<double, std::micro>(t1 - t0).count();
        ++callCount;
        ++r.calls;
        if (ok) {
          r.ok = true;
          r.Ltrc = gp;
          std::array<double, 3> const qi{xi, yi, zi};
          planeCoords(qi, r.xgTrc, r.ygTrc);
          r.offTrc = std::sqrt(r.xgTrc * r.xgTrc + r.ygTrc * r.ygTrc);
          r.inDishTrc = r.offTrc <= kTelescopeRadius_m;
        }
      }
      if (!r.ok) {
        r.status = "tracer: no path";
        rows.push_back(r);
        continue;
      }
      // tau: this tracer only accumulates it when absorption is enabled, and it then
      // samples survival internally from a time-seeded thread-local RNG.  Retry until
      // a call SURVIVES (the geometry and tau of a given ray are deterministic), and
      // score it only where the truth is thick enough for an absolute |dT| gate.
      double const Tref = std::exp(-r.tauRef);
      if (Tref >= kEpsPass) {
        int const maxCalls = std::min(400, std::max(1, static_cast<int>(2.0 / Tref)));
        for (int k = 0; k < maxCalls; ++k) {
          double xi = 0.0, yi = 0.0, zi = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
          double dhx = 0.0, dhy = 0.0, dhz = 0.0;  // arrival direction output
          auto const t0 = std::chrono::steady_clock::now();
          bool const ok = proc.rayTraceSphericalToPlane(
              p0[0], p0[1], p0[2], u[0], u[1], u[2], planePt[0], planePt[1], planePt[2],
              planeN[0], planeN[1], planeN[2], xi, yi, zi, op, gp, kLambda_nm, tau,
              dhx, dhy, dhz);
          auto const t1 = std::chrono::steady_clock::now();
          usSum += std::chrono::duration<double, std::micro>(t1 - t0).count();
          ++callCount;
          ++r.calls;
          if (ok) {
            r.tauTrc = tau;
            r.tauScored = true;
            break;
          }
        }
      }
      if (r.tauScored) {
        double const Ttrc = std::exp(-r.tauTrc);
        r.eps = (Tref >= kEpsPass) ? std::fabs(Ttrc - Tref) / Tref : std::fabs(Ttrc - Tref);
        r.status = "ok";
      } else {
        r.eps = std::numeric_limits<double>::quiet_NaN();
        r.status = "geometry ok, tau not scored (T_ref < 1e-2)";
      }
      rows.push_back(r);
    }
  }

  // ------------------------------------------------------------- report
  std::printf("\n  set A: forward/aim construction, %zu pairs (absorption on, tau scored)\n",
              rows.size());
  std::printf("  %7s %7s %8s %10s %10s %10s %10s %9s %9s %9s %6s\n", "h_em", "d_km",
              "arr_deg", "L_ref", "L_tracer", "dL/L", "y_off_ref", "y_off_trc", "eps", "tau_ref", "");
  int noPath = 0, fovTotal = 0, fovPass = 0, fovFlag = 0, fovFail = 0, dishBad = 0, guard = 0;
  double worstEps = 0.0, worstEpsH = 0.0, worstEpsD = 0.0;
  double maxDL = 0.0, maxImpact = 0.0, worstTauRel = 0.0;
  double worstImpactH = 0.0, worstImpactD = 0.0;
  for (auto const& r : rows) {
    if (!r.ok) { if (r.Lref > 0.0) ++noPath; std::printf("  %7.0f %7.0f %8.3f %10s %10s %10s %10s %9s %9s %9s %s\n", r.h, r.d / 1000.0, r.thetaArr, "-", "-", "-", "-", "-", "-", "-", r.status); continue; }
    double const dl = (r.Lref > 0.0) ? (r.Ltrc - r.Lref) / r.Lref : 0.0;
    maxDL = std::max(maxDL, std::fabs(dl));
    if (!std::isfinite(r.eps)) {
      std::printf("  %7.0f %7.0f %8.3f %10.2f %10.2f %10.2e %10.2f %9.2f %9s %9.4f %s\n", r.h,
                  r.d / 1000.0, r.thetaArr, r.Lref, r.Ltrc, dl, r.ygRef, r.ygTrc, "n/a", r.tauRef,
                  r.status);
      bool const inFov = std::fabs(r.thetaArr - 90.0) <= kFovHalf_deg;
      if (inFov) ++fovTotal;
      if (r.inDishRef != r.inDishTrc) {
        if (std::fabs(r.offRef - kTelescopeRadius_m) <= 0.5 ||
            std::fabs(r.offRef - kTelescopeRadiusWide_m) <= 0.5) ++guard;
        else ++dishBad;
      }
      continue;
    }
    double const dimp = std::sqrt((r.xgTrc - r.xgRef) * (r.xgTrc - r.xgRef) +
                                  (r.ygTrc - r.ygRef) * (r.ygTrc - r.ygRef));
    if (dimp > maxImpact) { maxImpact = dimp; worstImpactH = r.h; worstImpactD = r.d; }
    if (r.tauScored && r.tauRef > 1.0e-9)
      worstTauRel = std::max(worstTauRel, std::fabs(r.tauTrc - r.tauRef) / r.tauRef);
    std::printf("  %7.0f %7.0f %8.3f %10.2f %10.2f %10.2e %10.2f %9.2f %9.2e %9.4f %s\n",
                r.h, r.d / 1000.0, r.thetaArr, r.Lref, r.Ltrc, dl, r.ygRef, r.ygTrc, r.eps,
                r.tauRef, r.inDishRef != r.inDishTrc ? "DISH-MISMATCH" : "");
    bool const inFov = std::fabs(r.thetaArr - 90.0) <= kFovHalf_deg;
    if (inFov) {
      ++fovTotal;
      if (r.eps <= kEpsPass) ++fovPass;
      else if (r.eps <= kEpsAccept) ++fovFlag;
      else ++fovFail;
    }
    if (r.eps > worstEps) { worstEps = r.eps; worstEpsH = r.h; worstEpsD = r.d / 1000.0; }
    if (r.inDishRef != r.inDishTrc) {
      if (std::fabs(r.offRef - kTelescopeRadius_m) <= 0.5 ||
          std::fabs(r.offRef - kTelescopeRadiusWide_m) <= 0.5) ++guard;
      else ++dishBad;
    }
  }
  std::printf("    no path: %d | in-FoV %d: pass %d flagged %d fail %d | worst eps %.3e at %.0f m / "
              "%.0f km\n",
              noPath, fovTotal, fovPass, fovFlag, fovFail, worstEps, worstEpsH, worstEpsD);
  std::printf("    max |dL|/L = %.3e   max impact error = %.3e m (at %.0f m / %.0f km)   "
              "max |d tau|/tau = %.3e\n",
              maxDL, maxImpact, worstImpactH, worstImpactD / 1000.0, worstTauRel);
  std::printf("    dish verdict mismatches (r = %.0f m): %d (guard band 0.5 m: %d)\n",
              kTelescopeRadius_m, dishBad, guard);
  std::printf("    cost: %.1f us per call over %d calls\n",
              callCount > 0 ? usSum / callCount : 0.0, callCount);

  // ------------------------------------ n = 1 probe: analytic plane crossing
  std::printf("\n  diagnostic: n = 1 straight-ray probe (analytic plane crossing)\n");
  double worstOne = 0.0;
  {
    auto const refrOne = [](double, double) { return std::pair<double, double>{1.0, 1.0}; };
    cherenkov::CherenkovProcess const one(
        1.0, 399.0, 401.0, 1.0, kRdet, noRec, cherenkov::ProjectionMode::Telescope,
        std::shared_ptr<corsika::media::ShowerAxis const>{}, 0.0, 0.0, 0.0, planePt[0],
        planePt[1], planePt[2], kTelescopeRadius_m,
        std::function<std::pair<double, double>(double, double)>{refrOne}, kEarthRadius_m, 0.0,
        M_PI / 2.0, 0.0, 0.0, false, true, std::string{}, 0.0,
        std::shared_ptr<cherenkov::RefractiveIndexTable const>{}, false);
    std::printf("  %7s %7s %12s %12s %10s %14s\n", "h_em", "d_km", "t_analytic", "L_tracer",
                "dL/L", "|P_end-P_ref|");
    for (double const hKm : {2.0, 3.1, 5.0, 10.0, 20.0}) {
      for (double const dKm : {10.0, 100.0, 200.0}) {
        double const h = hKm * 1000.0, d = dKm * 1000.0;
        double const dPhi = d / kEarthRadius_m;
        double const rEm = kEarthRadius_m + h;
        std::array<double, 3> const p0{-rEm * std::sin(dPhi), 0.0, rEm * std::cos(dPhi)};
        double const chord = std::sqrt((planePt[0] - p0[0]) * (planePt[0] - p0[0]) +
                                       (planePt[2] - p0[2]) * (planePt[2] - p0[2]));
        std::array<double, 3> const u{(planePt[0] - p0[0]) / chord, 0.0,
                                      (planePt[2] - p0[2]) / chord};
        double const den = u[0] * planeN[0] + u[1] * planeN[1] + u[2] * planeN[2];
        double const num = (planePt[0] - p0[0]) * planeN[0] + (planePt[1] - p0[1]) * planeN[1] +
                           (planePt[2] - p0[2]) * planeN[2];
        double const tAna = num / den;
        double xi = 0.0, yi = 0.0, zi = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
        double dhx = 0.0, dhy = 0.0, dhz = 0.0;  // arrival direction output
        bool const ok = one.rayTraceSphericalToPlane(p0[0], p0[1], p0[2], u[0], u[1], u[2],
                                                     planePt[0], planePt[1], planePt[2], planeN[0],
                                                     planeN[1], planeN[2], xi, yi, zi, op, gp,
                                                     kLambda_nm, tau, dhx, dhy, dhz);
        double const miss = std::sqrt((xi - (p0[0] + tAna * u[0])) * (xi - (p0[0] + tAna * u[0])) +
                                      (yi - (p0[1] + tAna * u[1])) * (yi - (p0[1] + tAna * u[1])) +
                                      (zi - (p0[2] + tAna * u[2])) * (zi - (p0[2] + tAna * u[2])));
        double const dl = (gp - tAna) / tAna;
        worstOne = std::max(worstOne, std::fabs(dl));
        std::printf("  %7.0f %7.0f %12.3f %12.3f %10.2e %14.2e %s\n", h, dKm, tAna, gp, dl, miss,
                    ok ? "ok" : "no path");
      }
    }
    std::printf("  worst |dL|/L with n = 1: %.3e (analytic answer is 0)\n", worstOne);
  }
  CHECK(worstOne <= 1.0e-9);

  // ------------------- set B: continuous sweep of the aim through 90 deg ----
  // A spurious reflection / branch switch hides between grid points; sweep the
  // aim in 0.5 deg steps across the pointing direction and watch for a jump.
  std::printf("\n  set B: aim sweep through the pointing direction (geometry only, absorption "
              "off)\n");
  std::printf("  %7s %8s %10s %10s %10s %10s\n", "h_em", "aim_deg", "L_ref", "L_tracer",
              "dL/L", "dL jump");
  double sweepMaxDL = 0.0, sweepMaxJump = 0.0;
  for (double const hKm : {0.5, 3.0, 10.0}) {
    double prevL = 0.0, prevLtrc = 0.0;
    for (int k = -10; k <= 10; ++k) {
      double const aimDeg = 0.5 * k;
      double const aim = aimDeg * M_PI / 180.0;
      double const h = hKm * 1000.0, d = 50000.0;
      double const dPhi = d / kEarthRadius_m;
      double const rEm = kEarthRadius_m + h;
      std::array<double, 3> const p0{-rEm * std::sin(dPhi), 0.0, rEm * std::cos(dPhi)};
      double const chord = std::sqrt((planePt[0] - p0[0]) * (planePt[0] - p0[0]) +
                                     (planePt[2] - p0[2]) * (planePt[2] - p0[2]));
      double const ux0 = (planePt[0] - p0[0]) / chord;
      double const uz0 = (planePt[2] - p0[2]) / chord;
      std::array<double, 3> const u{ux0 * std::cos(aim) + uz0 * std::sin(aim), 0.0,
                                    -ux0 * std::sin(aim) + uz0 * std::cos(aim)};
      Truth const tr = traceTruthPlane(tab, row, p0, u, planePt, planeN);
      double xi = 0.0, yi = 0.0, zi = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
      double dhx = 0.0, dhy = 0.0, dhz = 0.0;  // arrival direction output
      bool const ok = geom.rayTraceSphericalToPlane(p0[0], p0[1], p0[2], u[0], u[1], u[2],
                                                    planePt[0], planePt[1], planePt[2], planeN[0],
                                                    planeN[1], planeN[2], xi, yi, zi, op, gp,
                                                    kLambda_nm, tau, dhx, dhy, dhz);
      if (!tr.ok || !ok) {
        std::printf("  %7.0f %8.2f %10s %10s %10s %10s  (truth %s / tracer %s)\n", hKm * 1000.0,
                    aimDeg, "-", "-", "-", "-", tr.ok ? "ok" : tr.why, ok ? "ok" : "no path");
        prevL = 0.0;
        prevLtrc = 0.0;
        continue;
      }
      double const dl = (tr.L > 0.0) ? (gp - tr.L) / tr.L : 0.0;
      sweepMaxDL = std::max(sweepMaxDL, std::fabs(dl));
      double jump = 0.0;
      if (prevL > 0.0) {
        jump = std::fabs((gp - prevLtrc) - (tr.L - prevL)) / tr.L;
        sweepMaxJump = std::max(sweepMaxJump, jump);
      }
      std::printf("  %7.0f %8.2f %10.2f %10.2f %10.2e %10.2e\n", hKm * 1000.0, aimDeg, tr.L, gp,
                  dl, jump);
      prevL = tr.L;
      prevLtrc = gp;
    }
  }
  std::printf("    sweep: max |dL|/L = %.3e   max step-to-step jump = %.3e\n", sweepMaxDL,
              sweepMaxJump);

  // ------------------- set C: off-plane azimuths (exercise the 3-D frame) ---
  std::printf("\n  set C: emission azimuth 37 deg and 137 deg off the telescope plane\n");
  std::printf("  %7s %7s %8s %10s %10s %10s %9s\n", "h_em", "d_km", "az_deg", "L_ref", "L_tracer",
              "dL/L", "x_off");
  double azMaxDL = 0.0, azMaxImpact = 0.0;
  int azBad = 0;
  for (double const azDeg : {37.0, 137.0}) {
    double const az = azDeg * M_PI / 180.0;
    for (double const hKm : {0.5, 3.1, 20.0}) {
      for (double const dKm : {10.0, 50.0, 200.0}) {
        double const h = hKm * 1000.0, d = dKm * 1000.0;
        double const dPhi = d / kEarthRadius_m;
        double const rEm = kEarthRadius_m + h;
        std::array<double, 3> const p0{-rEm * std::sin(dPhi) * std::cos(az),
                                       -rEm * std::sin(dPhi) * std::sin(az),
                                       rEm * std::cos(dPhi)};
        double const chord = std::sqrt((planePt[0] - p0[0]) * (planePt[0] - p0[0]) +
                                       (planePt[1] - p0[1]) * (planePt[1] - p0[1]) +
                                       (planePt[2] - p0[2]) * (planePt[2] - p0[2]));
        std::array<double, 3> const u{(planePt[0] - p0[0]) / chord, (planePt[1] - p0[1]) / chord,
                                      (planePt[2] - p0[2]) / chord};
        Truth const tr = traceTruthPlane(tab, row, p0, u, planePt, planeN);
        double xi = 0.0, yi = 0.0, zi = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
        double dhx = 0.0, dhy = 0.0, dhz = 0.0;  // arrival direction output
        bool const ok = geom.rayTraceSphericalToPlane(p0[0], p0[1], p0[2], u[0], u[1], u[2],
                                                      planePt[0], planePt[1], planePt[2],
                                                      planeN[0], planeN[1], planeN[2], xi, yi, zi,
                                                      op, gp, kLambda_nm, tau, dhx, dhy, dhz);
        if (!tr.ok || !ok) {
          ++azBad;
          std::printf("  %7.0f %7.0f %8.0f %10s %10s %10s %9s  (truth %s / tracer %s)\n", h, dKm,
                      azDeg, "-", "-", "-", "-", tr.ok ? "ok" : tr.why, ok ? "ok" : "no path");
          continue;
        }
        double const dl = (tr.L > 0.0) ? (gp - tr.L) / tr.L : 0.0;
        double const dimp = std::sqrt((xi - tr.p[0]) * (xi - tr.p[0]) +
                                      (yi - tr.p[1]) * (yi - tr.p[1]) +
                                      (zi - tr.p[2]) * (zi - tr.p[2]));
        azMaxDL = std::max(azMaxDL, std::fabs(dl));
        azMaxImpact = std::max(azMaxImpact, dimp);
        std::printf("  %7.0f %7.0f %8.0f %10.2f %10.2f %10.2e %9.2f\n", h, dKm, azDeg, tr.L, gp, dl,
                    dimp);
      }
    }
  }
  std::printf("    off-plane: max |dL|/L = %.3e   max impact error = %.2f m   failures = %d\n",
              azMaxDL, azMaxImpact, azBad);

  // ------------------------------------------------------------ assertions
  CHECK(noPath == 0);
  CHECK(fovFail == 0);
  CHECK(dishBad == 0);
  // Phase 1 gates.  Measured after the stage-local-kappa RK4 + Hermite closure:
  //   max |dL|/L            = 1.701e-12
  //   max impact error      = 6.184e-07 m
  //   max step-to-step jump = 3.231e-11
  // Before, the same set gave 4.117e-07 / 3.695e-01 m / 4.819e-08, i.e. the impact
  // error was 370x above the 1 mm target because the bending integrand froze n'/n at
  // the step start (first order in ds/H) and the plane crossing was interpolated
  // linearly along the chord.  The gates below are set 1e3-1e6 above the measurement
  // so that they test the physics, not the platform.
  CHECK(maxDL <= 1.0e-9);
  CHECK(maxImpact <= 1.0e-3);
  CHECK(sweepMaxJump <= 1.0e-9);
  CHECK(azBad == 0);
}

// ============================================================================
// Fast path (2g-h2m): the curved-ray TABLE in Telescope mode, with the optical
// depth tabulated inside the same lookup.  Only the cells the table can serve
// (emission above the telescope, downgoing arrival) are exercised; the rest is
// the documented fallback (EmissionBelowObs / Upgoing).
// ============================================================================
TEST_CASE("Cherenkov curved-ray table in telescope mode (fast path)", "modules,cherenkov") {
  logging::set_level(logging::level::warn);
  Table const tab = loadTable();
  REQUIRE(tab.loaded());
  std::size_t const row = tab.rowOf(kLambda_nm);

  // Refractive-index table from the same analytic profile, 0..130 km on 100 m.
  std::vector<double> altKm, lnN;
  for (double h = 0.0; h <= 130000.0; h += 100.0) {
    altKm.push_back(h / 1000.0);
    lnN.push_back(std::log(kRefrN0) - h / kRefrH_m);
  }
  auto const ri = std::make_shared<cherenkov::RefractiveIndexTable>(altKm, lnN, kEarthRadius_m);
  REQUIRE(ri->valid());

  std::array<double, 3> const planePt{0.0, 0.0, kRdet};
  std::array<double, 3> const planeN{1.0, 0.0, 0.0};
  cherenkov::CherenkovProcess::HitRecorder const noRec{};
  cherenkov::CherenkovProcess const proc(
      1.0003, 399.0, 401.0, 1.0, kRdet, noRec, cherenkov::ProjectionMode::Telescope,
      std::shared_ptr<corsika::media::ShowerAxis const>{}, 0.0, 0.0, 0.0, planePt[0], planePt[1],
      planePt[2], kTelescopeRadius_m,
      std::function<std::pair<double, double>(double, double)>{refrFunc}, kEarthRadius_m, 0.0,
      M_PI / 2.0, 0.0, 0.0, false /*dispersion*/, true /*curvature*/, tab.path, 0.0, ri,
      true /*enableCurvedTable*/);

  std::printf("\n[fast path] curved-ray table in Telescope mode, 90 deg pointing\n");
  std::printf("  %7s %7s %10s %10s %9s %9s %10s %10s %4s %s\n", "h_em", "d_km", "L_table",
              "L_truth", "dL/L", "impact_m", "tau_tab", "tau_ref", "tau_", "path");

  int hits = 0, fallbacks = 0;
  int hitsBelow = 0, fallbacksBelow = 0, contFail = 0, legacyDrop = 0;
  double maxDL = 0.0, maxImpact = 0.0, maxTauRel = 0.0, worstEps = 0.0;
  double maxDLBelow = 0.0, maxImpactBelow = 0.0, maxTauRelBelow = 0.0, worstEpsBelow = 0.0;
  double usLookup = 0.0, usBelow = 0.0;
  std::vector<double> usAll;
  std::size_t const nCells = sizeof(kDistances_km) / sizeof(double);
  // 0.5 / 2.0 / 2.9 km are BELOW the 3000 m telescope: those photons rise to it and
  // are served by the inward block.  3.1 / 5 / 10 / 20 km are the downward cells.
  for (double const hKm : {0.5, 2.0, 2.9, 3.1, 5.0, 10.0, 20.0}) {
    bool const below = hKm * 1000.0 < kRdet - kEarthRadius_m;
    for (std::size_t kd = 0; kd < nCells; ++kd) {
      double const dKm = kDistances_km[kd];
      double const h = hKm * 1000.0, d = dKm * 1000.0;
      double const dPhi = d / kEarthRadius_m;
      double const rEm = kEarthRadius_m + h;
      std::array<double, 3> const p0{-rEm * std::sin(dPhi), 0.0, rEm * std::cos(dPhi)};
      double const chord = std::sqrt((planePt[0] - p0[0]) * (planePt[0] - p0[0]) +
                                     (planePt[2] - p0[2]) * (planePt[2] - p0[2]));
      std::array<double, 3> const u{(planePt[0] - p0[0]) / chord, 0.0,
                                    (planePt[2] - p0[2]) / chord};
      Truth const tr = traceTruthPlane(tab, row, p0, u, planePt, planeN);
      if (!tr.ok) continue;
      double xg = 0.0, yg = 0.0, zg = 0.0, op = 0.0, gp = 0.0, tau = 0.0;
      bool usedTable = false;
      // Wave 2c attribution (Step A): which source answered this row's tau.  The two
      // counters live on the process, so a before/after delta names the source per row
      // without any new API and without trusting the printed aggregate.
      auto const stBefore = proc.curvedTableStats();
      auto const t0 = std::chrono::steady_clock::now();
      bool const hit = proc.traceCurvedGround(p0[0], p0[1], p0[2], u[0], u[1], u[2], kLambda_nm, xg,
                                              yg, zg, op, gp, tau, usedTable);
      auto const t1 = std::chrono::steady_clock::now();
      auto const stAfter = proc.curvedTableStats();
      char const* tauSrc = (stAfter.absTableTau > stBefore.absTableTau)     ? "abs"
                           : (stAfter.absFallbackTau > stBefore.absFallbackTau) ? "chrd"
                                                                              : "-";
      double const us = std::chrono::duration<double, std::micro>(t1 - t0).count();
      if (!hit) {
        // Two very different outcomes.  usedTable=true means the TABLE found the
        // sphere crossing and the plane continuation failed: that is a tracer defect
        // and must be zero.  usedTable=false means the lookup rejected the photon
        // (documented: descending-emission branch or inside the apparent horizon) and
        // the LEGACY plane tracer took over, which samples survival itself -- with
        // tau = 4..22 on these cells it mostly absorbs the photon, which is a
        // legitimate outcome, not a missing path.
        if (usedTable) {
          ++contFail;
          std::printf("  %7.0f %7.0f %10s %10.2f %9s %9s %10s %10.4f %4s  TABLE HIT, NO PLANE%s\n",
                      h, dKm, "-", tr.L, "-", "-", "-", tr.tau, tauSrc,
                      below ? " (BELOW)" : "");
        } else {
          ++legacyDrop;
          std::printf("  %7.0f %7.0f %10s %10.2f %9s %9s %10s %10.4f %4s  legacy drop (table "
                      "rejected)%s\n",
                      h, dKm, "-", tr.L, "-", "-", "-", tr.tau, tauSrc,
                      below ? " (BELOW)" : "");
        }
        continue;
      }
      (usedTable ? hits : fallbacks)++;
      if (below) (usedTable ? hitsBelow : fallbacksBelow)++;
      if (usedTable) usLookup += us;
      if (usedTable && below) usBelow += us;
      if (usedTable) usAll.push_back(us);
      double const dl = (tr.L > 0.0) ? (gp - tr.L) / tr.L : 0.0;
      double const dimp = std::sqrt((xg - tr.p[0]) * (xg - tr.p[0]) +
                                    (yg - tr.p[1]) * (yg - tr.p[1]) +
                                    (zg - tr.p[2]) * (zg - tr.p[2]));
      maxDL = std::max(maxDL, std::fabs(dl));
      maxImpact = std::max(maxImpact, dimp);
      double const Tref = std::exp(-tr.tau);
      double const Ttab = std::exp(-tau);
      double const eps = (Tref >= kEpsPass) ? std::fabs(Ttab - Tref) / Tref : std::fabs(Ttab - Tref);
      if (tr.tau > 1.0e-9) maxTauRel = std::max(maxTauRel, std::fabs(tau - tr.tau) / tr.tau);
      worstEps = std::max(worstEps, eps);
      if (below && usedTable) {
        maxDLBelow = std::max(maxDLBelow, std::fabs(dl));
        maxImpactBelow = std::max(maxImpactBelow, dimp);
        if (tr.tau > 1.0e-9)
          maxTauRelBelow = std::max(maxTauRelBelow, std::fabs(tau - tr.tau) / tr.tau);
        worstEpsBelow = std::max(worstEpsBelow, eps);
      }
      std::printf("  %7.0f %7.0f %10.1f %10.1f %9.1e %9.2f %10.4f %10.4f %4s %s%s  %8.2f us\n", h,
                  dKm, gp, tr.L, dl, dimp, tau, tr.tau, tauSrc, usedTable ? "table" : "fallback",
                  below ? " (BELOW)" : "", us);
    }
  }
  {
    auto const st = proc.curvedTableStats();
    std::printf("    table: available=%d attempted=%d hits=%zu bypassed=%zu lookups=%zu "
                "legacyOK=%zu legacyFail=%zu reason='%s'\n",
                st.tableAvailable ? 1 : 0, st.tableAttempted ? 1 : 0, st.hits, st.bypassed,
                st.attempted, st.legacyFallbackSucceeded, st.legacyFallbackFailed,
                st.unavailableReason.c_str());
    std::printf("    gate: passed=%zu rejected=%zu   table hits above=%zu below=%zu\n",
                st.gatePassed, st.gateRejected, st.tableHitAbove, st.tableHitBelow);
    std::printf("    tau source: absorption table=%zu   exact chord (last resort)=%zu   "
                "missing (no abs table)=%zu\n",
                st.absTableTau, st.absFallbackTau, st.tauMissing);
    // Step A: attribute every ray the absorption table refused.  The buffer holds the
    // last 32 rejections in call order, so these ARE the acceptance set's rejects.
    auto const rej = proc.absorptionRejections();
    std::printf("    absorption-table rejections recorded: %zu\n", rej.size());
    for (std::size_t i = 0; i < rej.size() && i < 32; ++i) {
      std::printf("    abs-reject[%zu]: hA=%7.0f hB=%7.0f xi=%.7f mu=%.3e lam=%.0f  %s\n", i,
                  rej[i].hA, rej[i].hB, rej[i].xi,
                  std::sqrt(std::max(0.0, 1.0 - rej[i].xi * rej[i].xi)), rej[i].lambdaNm,
                  cherenkov::CurvedAbsorptionTable::statusName(rej[i].status));
    }
    for (int k = 0; k < 9; ++k) {
      if (st.rejected[k] > 0) {
        std::printf("    rejected[%d] = %zu\n", k, st.rejected[k]);
      }
    }
  }
  std::printf("    table hits %d, fallbacks %d\n", hits, fallbacks);
  std::printf("    below-observation rows (inward block): table %d, fallback %d, continuation "
              "failures %d, legacy drops %d\n",
              hitsBelow, fallbacksBelow, contFail, legacyDrop);
  std::printf("    max |dL|/L = %.3e   max impact error = %.2f m   max |dtau|/tau = %.3e   "
              "worst eps = %.3e\n",
              maxDL, maxImpact, maxTauRel, worstEps);
  std::printf("    below: max |dL|/L = %.3e   max impact = %.2f m   max |dtau|/tau = %.3e   "
              "worst eps = %.3e\n",
              maxDLBelow, maxImpactBelow, maxTauRelBelow, worstEpsBelow);
  // The first lookup carries the one-time table build (~40 s); report the median of
  // the whole route (table lookup + the short plane continuation) instead of a mean,
  // which a single outlier would dominate.
  std::sort(usAll.begin(), usAll.end());
  double const usMed = usAll.empty() ? 0.0 : usAll[usAll.size() / 2];
  double const usP90 = usAll.empty() ? 0.0 : usAll[(usAll.size() * 9) / 10];
  std::printf("    table-route cost: median %.0f ns, p90 %.0f ns, max %.0f ns over %zu lookups "
              "(mean %.0f) -- the first one builds the table (~40 s)\n",
              usMed * 1.0e3, usP90 * 1.0e3, usAll.empty() ? 0.0 : usAll.back() * 1.0e3,
              usAll.size(), usAll.empty() ? 0.0 : 1.0e3 * usLookup / double(usAll.size()));
  CHECK(hits > 0);
  CHECK(hitsBelow > 0);
  CHECK(contFail == 0);
  CHECK(maxDL <= 1.0e-3);
  CHECK(maxImpact <= 1.0);
  CHECK(maxTauRel <= 5.0e-2);
  CHECK(maxDLBelow <= 1.0e-3);
  CHECK(maxImpactBelow <= 1.0);
  CHECK(maxTauRelBelow <= 5.0e-2);
}
