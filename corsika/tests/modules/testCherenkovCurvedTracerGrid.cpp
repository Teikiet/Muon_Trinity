/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * \file testCherenkovCurvedTracerGrid.cpp
 *
 * Step 2 (2g-h2h): the existing curved tracer with no table, checked against an
 * independent refracted-path reference on an endpoint-defined grid.
 *
 * ACCURACY RULE (agreed): relative transmission error
 *     eps = |T_prod - T_ref| / T_ref      (|dT| when T_ref < 0.01)
 *     eps <= 1 %  pass ; 1-5 % acceptable (flagged) ; > 5 % fail.
 *
 * PRODUCTION UNDER TEST (read-only):
 *   CherenkovProcess::rayTraceSpherical  (Cherenkov.inl:896-1234) -- the legacy
 *   spherical tracer, used by traceCurvedGround whenever the curved table is
 *   not eligible (absorption on, or enableCurvedTable = false).  It is called
 *   with curvature ON and enableCurvedTable = false, exactly as
 *   Cherenkov.inl:1629-1650 does on the legacy route, and it returns
 *   accumulated_tau without sampling survival (Cherenkov.inl:2029-2031).
 *   Its extinction comes from extinctionCoeffFunc_ = a +-25 m window average of
 *   the cumulative tau (defined Cherenkov.inl:96-113, used at :1121, :1177, :1204).
 *
 * REFERENCE (this file only): the step-1 spherical integrator plus refraction.
 *   A spherically symmetric index n(r) = 1 + N0 exp(-(r-R)/H) (N0 = 3.27e-4,
 *   H = 8 km -- the profile already used by testCherenkovCurvedRayTable), ray
 *   equation  du/ds = (n'/n)(r_hat - (u.r_hat)u)  integrated with RK4 at
 *   ds = 1 m, tau = int alpha ds with alpha = dtau/dh from the same MODTRAN
 *   table (400 nm row).  Halving the step is reported for every cell.
 *
 * GRID: emission {0.5,2,5,10,20} km x detector {0.5,2,5,10,20,35} km x ground
 * distance {0,1,5,20,50,100,200,400} km.  Cells are defined by their endpoints,
 * so every cell is reachable and same-altitude horizontal paths are included.
 */

#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/Cherenkov.hpp>

#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <catch2/catch_all.hpp>

#include <PhysicalUnitsCatch2.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace corsika;

namespace {

  constexpr double kEarthRadius_m = 6.371e6;
  constexpr double kLambda_nm = 400.0;
  constexpr double kStepRef_m = 1.0;   // reference RK4 step
  constexpr double kRefrN0 = 3.27e-4;  // refractivity at sea level
  constexpr double kRefrH_m = 8000.0;  // scale height
  constexpr double kTauRefDet_m = 0.216105;
  constexpr double kEpsPass = 0.01;
  constexpr double kEpsAccept = 0.05;

  constexpr double kEmissions_km[] = {0.5, 2.0, 5.0, 10.0, 20.0};
  constexpr double kDetectors_km[] = {0.5, 2.0, 5.0, 10.0, 20.0, 35.0};
  constexpr double kDistances_km[] = {0.0, 1.0, 5.0, 20.0, 50.0, 100.0, 200.0, 400.0};

  // ------------------------------------------------------------------ table
  struct Table {
    std::size_t nLambda{0};
    std::size_t nHeight{0};
    std::vector<double> height_m;
    std::vector<std::vector<double>> tau; // [ilambda][iheight], cumulative from 0 m
    std::string path;
    bool loaded() const { return nLambda > 0 && nHeight > 0; }

    std::size_t slab(double const h_m) const {
      std::size_t ih = 0;
      while (ih + 1 < height_m.size() && height_m[ih + 1] <= h_m) ++ih;
      return ih;
    }
    /// Piecewise-constant exact extinction of the 400 nm row.
    double alpha(double const h_m) const {
      std::size_t const ih = slab(h_m);
      if (ih + 1 >= height_m.size()) return 0.0;
      return (tau[44][ih + 1] - tau[44][ih]) / (height_m[ih + 1] - height_m[ih]);
    }
    double tauAt(double const h_m) const {
      std::size_t const ih = slab(h_m);
      if (ih + 1 >= height_m.size()) return tau[44][ih];
      double const dh = height_m[ih + 1] - height_m[ih];
      double const f = (h_m - height_m[ih]) / dh;
      return tau[44][ih] * (1.0 - f) + tau[44][ih + 1] * f;
    }
    /// The production window: (tau(h+dh/2) - tau(h-dh/2)) / dh with dh = 50 m,
    /// clamped to the table bounds -- identical to Cherenkov.inl:96-113.
    double alphaWindow(double const h_m, double const dh = 50.0) const {
      double const hLo = std::max(height_m.front(), h_m - 0.5 * dh);
      double const hHi = std::min(height_m.back(), h_m + 0.5 * dh);
      double const actual = hHi - hLo;
      if (actual <= 0.0) return 0.0;
      return std::max(0.0, (tauAt(hHi) - tauAt(hLo)) / actual);
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
      std::size_t nl = 0;
      std::size_t nh = 0;
      if (!(in >> nl >> nh) || nl < 45 || nh == 0) continue;
      Table t;
      t.nLambda = nl;
      t.nHeight = nh;
      std::vector<double> lam(nl);
      t.height_m.resize(nh);
      bool ok = true;
      for (std::size_t i = 0; i < nl && ok; ++i) ok = static_cast<bool>(in >> lam[i]);
      for (std::size_t i = 0; i < nh && ok; ++i) ok = static_cast<bool>(in >> t.height_m[i]);
      t.tau.assign(nl, std::vector<double>(nh, 0.0));
      for (std::size_t i = 0; i < nl && ok; ++i) {
        for (std::size_t j = 0; j < nh && ok; ++j) ok = static_cast<bool>(in >> t.tau[i][j]);
      }
      if (!ok) continue;
      t.path = path;
      return t;
    }
    return Table{};
  }

  // ------------------------------------------------------- refraction model
  double refrIndex(double const h_m) {
    return 1.0 + kRefrN0 * std::exp(-std::max(0.0, h_m) / kRefrH_m);
  }
  double refrDnDr(double const h_m) {
    return -kRefrN0 / kRefrH_m * std::exp(-std::max(0.0, h_m) / kRefrH_m);
  }
  std::pair<double, double> refrFunc(double const h_m, double) {
    double const n = refrIndex(h_m);
    return {n, n};
  }


  // ------------------------------------------------ reference ray integrator
  enum class Stop { AtRadius, AtArc };

  struct RefOut {
    bool ok{false};
    char const* why{""};
    double tau{0.0};        // exact alpha
    double tauWindow{0.0};  // production +-25 m window
    double L_m{0.0};        // path length
    double arc_m{0.0};      // ground arc distance
    double hMin_m{0.0};
    double hEnd_m{0.0};
    double invariantDrift{0.0};
    std::size_t steps{0};
  };

  /// RK4 for the ray equation in a spherically symmetric medium, plane coords
  /// (xp, zp) in metres about the Earth centre, unit direction (ux, uz).
  RefOut traceRefracted(Table const& tab, double const xStart, double const zStart,
                        double const uxStart, double const uzStart, Stop const stop,
                        double const stopValue, double const ds) {
    RefOut out;
    double xp = xStart, zp = zStart, ux = uxStart, uz = uzStart;
    double r = std::hypot(xp, zp);
    out.hMin_m = r - kEarthRadius_m;
    double arc = 0.0;
    double s = 0.0;
    bool const descending = (uz * (zp / r) + ux * (xp / r)) < 0.0;
    double const p0 = refrIndex(out.hMin_m) * r * std::abs(ux * (zp / r) - uz * (xp / r));

    auto deriv = [&](double const x, double const z, double const vx, double const vz,
                     double& dx, double& dz, double& dvx, double& dvz) {
      double const rr = std::hypot(x, z);
      double const h = rr - kEarthRadius_m;
      double const k = refrDnDr(h) / refrIndex(h);
      double const rx = x / rr, rz = z / rr;
      double const ur = vx * rx + vz * rz;
      dx = vx;
      dz = vz;
      dvx = k * (rx - ur * vx);
      dvz = k * (rz - ur * vz);
    };

    for (std::size_t iter = 0; iter < 2000000ULL; ++iter) {
      if (iter + 1 == 2000000ULL) {
        out.why = "path longer than 2000 km";
        break;
      }
      double const h = r - kEarthRadius_m;
      if (h < -100.0) {
        out.why = "ray below sea level";
        break;
      }
      if (stop == Stop::AtArc && out.arc_m >= stopValue) {
        out.ok = true;
        break;
      }
      double dax, daz, davx, davz;
      deriv(xp, zp, ux, uz, dax, daz, davx, davz);
      double dbx, dbz, dbvx, dbvz;
      deriv(xp + 0.5 * ds * dax, zp + 0.5 * ds * daz, ux + 0.5 * ds * davx, uz + 0.5 * ds * davz,
            dbx, dbz, dbvx, dbvz);
      double dcx, dcz, dcvx, dcvz;
      deriv(xp + 0.5 * ds * dbx, zp + 0.5 * ds * dbz, ux + 0.5 * ds * dbvx, uz + 0.5 * ds * dbvz,
            dcx, dcz, dcvx, dcvz);
      double ddx, ddz, ddvx, ddvz;
      deriv(xp + ds * dcx, zp + ds * dcz, ux + ds * dcvx, uz + ds * dcvz, ddx, ddz, ddvx, ddvz);

      double const hPrev = h;
      xp += ds / 6.0 * (dax + 2.0 * dbx + 2.0 * dcx + ddx);
      zp += ds / 6.0 * (daz + 2.0 * dbz + 2.0 * dcz + ddz);
      ux += ds / 6.0 * (davx + 2.0 * dbvx + 2.0 * dcvx + ddvx);
      uz += ds / 6.0 * (davz + 2.0 * dbvz + 2.0 * dcvz + ddvz);
      double const norm = std::hypot(ux, uz);
      ux /= norm;
      uz /= norm;
      r = std::hypot(xp, zp);
      double const hNew = r - kEarthRadius_m;
      out.hMin_m = std::min(out.hMin_m, hNew);

      // tau increment, exact in altitude: along one step dh/ds is constant, so
      // int alpha ds = |tau(h_new) - tau(h_prev)| * ds / |dh|, which uses the
      // table's cumulative tau directly and removes the staircase quadrature error.
      double const hMid = 0.5 * (hPrev + hNew);
      double const dhStep = hNew - hPrev;
      double tauStep = 0.0;
      double tauWinStep = 0.0;
      if (std::fabs(dhStep) > 1.0e-6) {
        tauStep = std::fabs(tab.tauAt(hNew) - tab.tauAt(hPrev)) * ds / std::fabs(dhStep);
      } else {
        tauStep = tab.alpha(hMid) * ds;
      }
      tauWinStep = 0.5 * (tab.alphaWindow(hPrev) + tab.alphaWindow(hNew)) * ds;
      out.tau += tauStep;
      out.tauWindow += tauWinStep;

      // ground arc
      double const cross = std::abs(ux * (zp / r) - uz * (xp / r));
      arc += cross / r * ds * kEarthRadius_m;
      out.arc_m = arc;
      s += ds;
      out.steps = iter + 1;
      out.L_m = s;
      out.hEnd_m = hNew;

      // radius crossing (interpolated) when stopping on radius
      if (stop == Stop::AtRadius && s > 0.0) {
        double const ur = ux * (xp / r) + uz * (zp / r);
        if ((descending && ur >= 1.0e-12) || (!descending && ur <= -1.0e-12)) {
          out.why = "ray turned before reaching the target altitude";
          break;
        }
        bool const crossed =
            (descending && hNew <= stopValue) || (!descending && hNew >= stopValue);
        bool const wasOnSide = descending ? (hPrev >= stopValue) : (hPrev <= stopValue);
        if (crossed && wasOnSide) {
          double const frac = (hPrev - stopValue) / (hPrev - hNew);
          double const keep = std::max(0.0, std::min(1.0, frac));
          // replace the full step by the part that stops at the target altitude
          out.tau -= tauStep;
          out.tauWindow -= tauWinStep;
          if (std::fabs(stopValue - hPrev) > 1.0e-6) {
            out.tau += std::fabs(tab.tauAt(stopValue) - tab.tauAt(hPrev)) * (ds * keep) /
                       std::fabs(stopValue - hPrev);
          } else {
            out.tau += tab.alpha(hPrev) * ds * keep;
          }
          out.tauWindow +=
              0.5 * (tab.alphaWindow(hPrev) + tab.alphaWindow(stopValue)) * ds * keep;
          out.L_m = s - ds * (1.0 - keep);
          out.hEnd_m = stopValue;
          out.ok = true;
          break;
        }
      }
    }
    if (out.ok) {
      double const pEnd = refrIndex(out.hEnd_m) * (kEarthRadius_m + out.hEnd_m) *
                          std::abs(ux * (zp / r) - uz * (xp / r));
      out.invariantDrift = (p0 > 0.0) ? std::abs(pEnd - p0) / p0 : 0.0;
    }
    return out;
  }


  // -------------------------------------------------------------- tracer call
  struct TracerOut {
    bool ok{false};
    double tau{0.0};
    double geomPath_m{0.0};
  };

  TracerOut runTracer(cherenkov::CherenkovProcess const& proc, double const h1_m,
                      double const ux, double const uz) {
    double xg = 0.0, yg = 0.0, zg = 0.0, optical = 0.0, geom = 0.0, tau = 0.0;
    TracerOut out;
    out.ok = proc.rayTraceSpherical(0.0, 0.0, kEarthRadius_m + h1_m, ux, 0.0, uz, xg, yg, zg,
                                    optical, geom, kLambda_nm, tau);
    out.tau = tau;
    out.geomPath_m = geom;
    return out;
  }

  /// Chord direction and lengths for the endpoint-defined cell.
  struct CellGeom {
    double ux{0.0}, uz{0.0};   // unit direction at emission (chord to the endpoint)
    double chord_m{0.0};
    double thetaDeg{0.0};      // zenith angle at emission
  };

  CellGeom cellGeometry(double const h1_m, double const h2_m, double const d_m) {
    double const dPhi = d_m / kEarthRadius_m;
    double const r1 = kEarthRadius_m + h1_m;
    double const r2 = kEarthRadius_m + h2_m;
    double const x1 = 0.0, z1 = r1;
    double const x2 = r2 * std::sin(dPhi), z2 = r2 * std::cos(dPhi);
    CellGeom g;
    g.chord_m = std::hypot(x2 - x1, z2 - z1);
    g.ux = (x2 - x1) / g.chord_m;
    g.uz = (z2 - z1) / g.chord_m;
    g.thetaDeg =
        std::acos(std::max(-1.0, std::min(1.0, -g.uz))) * 180.0 / M_PI; // CORSIKA zenith
    return g;
  }

  struct Worst {
    double eps{0.0};
    double h1_m{0.0}, h2_m{0.0}, d_m{0.0};
    char const* note{""};
  };


} // namespace

TEST_CASE("Cherenkov curved tracer vs refracted reference", "modules,cherenkov") {
  logging::set_level(logging::level::warn);

  Table const tab = loadTable();
  REQUIRE(tab.loaded());
  CHECK(tab.tauAt(2000.0) == Approx(kTauRefDet_m).margin(1.0e-6));

  std::printf("\n[CherenkovTracer] table %s\n", tab.path.c_str());
  std::printf("  n(h) = 1 + %.3e exp(-h/%.0f m); reference RK4 ds = %.1f m; lambda = %.0f nm\n",
              kRefrN0, kRefrH_m, kStepRef_m, kLambda_nm);
  std::printf("  production window +-25 m: Cherenkov.inl:96-113 (extinctionCoeffFunc_), "
              "used at :1121, :1177, :1204\n");
  std::printf("  grid: emission %zu x detector %zu x distance %zu = %zu cells\n",
              sizeof(kEmissions_km) / sizeof(double), sizeof(kDetectors_km) / sizeof(double),
              sizeof(kDistances_km) / sizeof(double),
              sizeof(kEmissions_km) / sizeof(double) * sizeof(kDetectors_km) / sizeof(double) *
                  sizeof(kDistances_km) / sizeof(double));

  int cells = 0, skipped = 0, notExpressible = 0, noPath = 0, nonFinite = 0;
  int pass = 0, accept = 0, fail = 0;
  double maxWindowBias = 0.0, maxWindowBiasAt[3] = {0, 0, 0};
  double maxHalving = 0.0, maxInvariant = 0.0;
  double halvingAt[3] = {0.0, 0.0, 0.0};
  Worst worst;
  Worst worstNonHorizontal;

  std::printf("  %6s %6s %7s %7s %9s %9s %11s %11s %9s %s\n", "h_em", "h_det", "d_km",
              "theta", "L_ref_m", "L_trc_m", "T_tracer", "T_ref", "eps", "class");

  for (double const hDet_km : kDetectors_km) {
    double const hDet = hDet_km * 1000.0;
    // one process per detector altitude: this is the plane the tracer stops on
    cherenkov::CherenkovProcess::HitRecorder const noRecorder{};
    cherenkov::CherenkovProcess proc(
        1.0003, 399.0, 401.0, 1.0, kEarthRadius_m + hDet, noRecorder,
        cherenkov::ProjectionMode::Ground, std::shared_ptr<media::ShowerAxis const>{}, 0.0, 0.0,
        0.0, 0.0, 0.0, kEarthRadius_m + hDet, 0.0,
        std::function<std::pair<double, double>(double, double)>{refrFunc}, kEarthRadius_m, -1.0,
        -1.0, 0.0, 0.0, false /*dispersion*/, true /*enableCurvature*/, tab.path, 0.0,
        std::shared_ptr<cherenkov::RefractiveIndexTable const>{}, false /*curvedTable*/);
    (void)proc;

    for (double const hEm_km : kEmissions_km) {
      double const hEm = hEm_km * 1000.0;
      for (double const d_km : kDistances_km) {
        double const d = d_km * 1000.0;
        if (hEm == hDet && d == 0.0) {
          ++skipped;
          continue;
        }
        ++cells;
        bool const horizontal = (hEm == hDet);
        CellGeom const g = cellGeometry(hEm, hDet, d);
        RefOut ref;
        TracerOut trc;
        double const ux = horizontal ? 1.0 : g.ux;
        double const uz = horizontal ? 0.0 : g.uz;
        if (horizontal) {
          ref = traceRefracted(tab, 0.0, kEarthRadius_m + hEm, 1.0, 0.0, Stop::AtArc, d,
                               kStepRef_m);
        } else {
          ref = traceRefracted(tab, 0.0, kEarthRadius_m + hEm, g.ux, g.uz, Stop::AtRadius, hDet,
                               kStepRef_m);
        }
        trc = runTracer(proc, hEm, ux, uz);
        if (!std::isfinite(trc.tau) || !std::isfinite(ref.tau)) {
          ++nonFinite;
        }
        // reference step-halving check
        RefOut const refHalf =
            horizontal ? traceRefracted(tab, 0.0, kEarthRadius_m + hEm, 1.0, 0.0, Stop::AtArc, d,
                                        0.5 * kStepRef_m)
                       : traceRefracted(tab, 0.0, kEarthRadius_m + hEm, g.ux, g.uz,
                                        Stop::AtRadius, hDet, 0.5 * kStepRef_m);
        if (ref.ok && refHalf.ok && ref.tau > 1.0e-9) {
          double const hv = std::fabs(ref.tau - refHalf.tau) / ref.tau;
          if (hv > maxHalving) {
            maxHalving = hv;
            halvingAt[0] = hEm;
            halvingAt[1] = hDet;
            halvingAt[2] = d;
          }
        }
        maxInvariant = std::max(maxInvariant, ref.invariantDrift);
        if (ref.ok && ref.tau > 1.0e-9) {
          double const bias = std::fabs(ref.tauWindow - ref.tau) / ref.tau;
          if (bias > maxWindowBias) {
            maxWindowBias = bias;
            maxWindowBiasAt[0] = hEm;
            maxWindowBiasAt[1] = hDet;
            maxWindowBiasAt[2] = d;
          }
        }
        if (!ref.ok) {
          std::printf(" %6.0f %6.0f %7.0f %7s %9s %9s %11s %11s %9s reference: %s\n", hEm, hDet,
                      d, "-", "-", "-", "-", "-", "-", ref.why);
          ++fail;
          continue;
        }
        double const tRef = std::exp(-ref.tau);
        double const tTrc = trc.ok ? std::exp(-trc.tau) : std::numeric_limits<double>::quiet_NaN();
        char const* cls = "";
        double eps = 0.0;
        if (!trc.ok) {
          ++noPath;
          cls = (horizontal ? "NO-PATH (horizontal)" : "NO-PATH");
          ++fail;
        } else if (horizontal) {
          ++notExpressible;
          eps = (tRef >= kEpsPass) ? std::fabs(tTrc - tRef) / tRef : std::fabs(tTrc - tRef);
          cls = "NOT-EXPRESSIBLE";
          ++fail;
        } else {
          eps = (tRef >= kEpsPass) ? std::fabs(tTrc - tRef) / tRef : std::fabs(tTrc - tRef);
          if (eps <= kEpsPass) {
            ++pass;
            cls = "pass";
          } else if (eps <= kEpsAccept) {
            ++accept;
            cls = "ACCEPT";
          } else {
            ++fail;
            cls = "FAIL";
          }
        }
        if (std::isfinite(eps) && eps > worst.eps) {
          worst.eps = eps;
          worst.h1_m = hEm;
          worst.h2_m = hDet;
          worst.d_m = d;
          worst.note = cls;
        }
        if (!horizontal && trc.ok && std::isfinite(eps) && eps > worstNonHorizontal.eps) {
          worstNonHorizontal.eps = eps;
          worstNonHorizontal.h1_m = hEm;
          worstNonHorizontal.h2_m = hDet;
          worstNonHorizontal.d_m = d;
          worstNonHorizontal.note = cls;
        }
        std::printf(" %6.0f %6.0f %7.0f %7.2f %9.1f %9.1f %11.6f %11.6f %9.2e %s\n", hEm, hDet, d,
                    g.thetaDeg, ref.L_m, trc.ok ? trc.geomPath_m : -1.0, tTrc, tRef, eps, cls);
      }
    }
  }
  (void)notExpressible;
  (void)noPath;

  // ---------------------------------------------------------------- timings
  double tVertical_us = 0.0;
  double tGrazing_us = 0.0;
  double tauVertical = 0.0;
  double tauGrazing = 0.0;
  {
    cherenkov::CherenkovProcess::HitRecorder const noRecorder{};
    cherenkov::CherenkovProcess tproc(
        1.0003, 399.0, 401.0, 1.0, kEarthRadius_m + 2000.0, noRecorder,
        cherenkov::ProjectionMode::Ground, std::shared_ptr<media::ShowerAxis const>{}, 0.0, 0.0,
        0.0, 0.0, 0.0, kEarthRadius_m + 2000.0, 0.0,
        std::function<std::pair<double, double>(double, double)>{refrFunc}, kEarthRadius_m, -1.0,
        -1.0, 0.0, 0.0, false, true, tab.path, 0.0,
        std::shared_ptr<cherenkov::RefractiveIndexTable const>{}, false);
    auto bench = [&](double const h1_m, double const d_m, double& tauOut) {
      CellGeom const g = cellGeometry(h1_m, 2000.0, d_m);
      std::size_t const n = 200;
      auto const t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; i < n; ++i) {
        double xg = 0, yg = 0, zg = 0, opt = 0, geom = 0, tau = 0;
        tproc.rayTraceSpherical(0.0, 0.0, kEarthRadius_m + h1_m, g.ux, 0.0, g.uz, xg, yg, zg,
                                opt, geom, kLambda_nm, tau);
        tauOut += tau;
      }
      auto const t1 = std::chrono::steady_clock::now();
      tauOut /= static_cast<double>(n);
      return std::chrono::duration<double, std::micro>(t1 - t0).count() /
             static_cast<double>(n);
    };
    tVertical_us = bench(10000.0, 0.0, tauVertical);
    tGrazing_us = bench(5000.0, 100000.0, tauGrazing);
  }

  // ---------------------------------------------------------------- verdict
  int const classified = pass + accept + fail;
  std::printf("\n[CherenkovTracer] verdict (rule: eps = |dT|/T_ref, <=1 %% pass, 1-5 %% "
              "acceptable, >5 %% fail)\n");
  std::printf("  cells %d (+%d zero-length skipped): pass %d, acceptable %d, fail %d "
              "(of which horizontal not-expressible %d, no-path %d)\n",
              cells, skipped, pass, accept, fail, notExpressible, noPath);
  std::printf("  pass rate: %.1f %% within 1 %%, %.1f %% within 5 %%\n",
              classified > 0 ? 100.0 * pass / classified : 0.0,
              classified > 0 ? 100.0 * (pass + accept) / classified : 0.0);
  std::printf("  worst eps = %.3e at h_emit = %.0f m, h_det = %.0f m, d = %.0f m  [%s]\n",
              worst.eps, worst.h1_m, worst.h2_m, worst.d_m, worst.note);
  std::printf("  worst non-horizontal eps = %.3e at h_emit = %.0f m, h_det = %.0f m, "
              "d = %.0f m  [%s]\n",
              worstNonHorizontal.eps, worstNonHorizontal.h1_m, worstNonHorizontal.h2_m,
              worstNonHorizontal.d_m, worstNonHorizontal.note);
  std::printf("  +-25 m window bias: max |tau_win - tau_exact|/tau_exact = %.3e at "
              "h_emit = %.0f m, h_det = %.0f m, d = %.0f m\n",
              maxWindowBias, maxWindowBiasAt[0], maxWindowBiasAt[1], maxWindowBiasAt[2]);
  std::printf("  reference convergence: max |dtau(1 m) - dtau(0.5 m)|/dtau = %.3e at h_emit = %.0f m, "
              "h_det = %.0f m, d = %.0f m\n",
              maxHalving, halvingAt[0], halvingAt[1], halvingAt[2]);
  std::printf("  Bouguer invariant drift (n r sin(psi)): max %.3e\n", maxInvariant);
  std::printf("  non-finite cells = %d\n", nonFinite);
  std::printf("  t_tracer: vertical 10->2 km %.2f us/call (tau = %.6f), grazing 5->2 km "
              "d = 100 km %.2f us/call (tau = %.6f)\n",
              tVertical_us, tauVertical, tGrazing_us, tauGrazing);

  CHECK(nonFinite == 0);
  // The straight-chord reference of step 1 converged to 7e-12.  On refracted
  // rays whose endpoint is reached at a shallow angle the crossing interpolation
  // dominates and 1e-6 is not reachable with a 1 m step; 1e-3 is 10x below the
  // 1 % budget and is the level used here.
  CHECK(maxHalving < 1.0e-3);
  CHECK(maxInvariant < 1.0e-6);
}


