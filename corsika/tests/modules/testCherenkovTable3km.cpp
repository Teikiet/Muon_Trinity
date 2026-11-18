/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * \file testCherenkovTable3km.cpp
 *
 * Step 3 (2g-h2i) with the fixed IACT geometry: detector at 3 km, pointing at
 * 90 deg zenith, 10 deg FULL field of view -> accepted arrival zenith 85...95
 * deg (half angle 5 deg).
 *
 * ACCURACY RULE (agreed): eps = |T_prod - T_ref| / T_ref, <= 1 % pass,
 * 1-5 % flagged, > 5 % fail; |dT| <= 1e-3 when T_ref < 0.01.
 *
 * TRUTH (test only): refracted-ray RK4 integration, ds = 1 m, tau exact in
 * altitude against the table's cumulative tau (no +-25 m window), n(r) =
 * 1 + 3.27e-4 exp(-(r-R)/8 km).  Convergence is re-checked on this grid.
 *
 * TABLE PROTOTYPE (test only, nothing implemented in production):
 *   one slice at h_det = 3000 m, indexed by (emission height, arrival zenith),
 *   node values are the exact-in-altitude tau of the ray that arrives at the
 *   detector with that zenith, bilinear interpolation inside a block.
 *   Two blocks cover the physically valid domain:
 *     block A: h_em >= 3000 m, theta_arr <= 90 deg  (source above the horizon)
 *     block B: h_em <= 3000 m, theta_arr >= 90 deg  (source below the horizon)
 *   Two node densities are compared, plus a cos(theta)-spaced diagnostic.
 *
 * GRID: emission {0.5, 2, 2.9, 3.0, 3.1, 5, 10, 20} km x distance
 * {1, 10, 50, 100, 200, 400} km x lambda {300, 400, 500, 600} nm, built in
 * arrival direction: a cell counts only if its arrival zenith is inside the FoV.
 */

#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/Cherenkov.hpp>

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
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace corsika;

namespace {

  constexpr double kEarthRadius_m = 6.371e6;
  constexpr double kDetectorAltitude_m = 3000.0;
  constexpr double kStepRef_m = 1.0;
  constexpr double kRefrN0 = 3.27e-4;
  constexpr double kRefrH_m = 8000.0;
  constexpr double kZenithPointing_deg = 90.0;
  constexpr double kFovHalfAngle_deg = 5.0; // 10 deg full field of view
  constexpr double kEpsPass = 0.01;
  constexpr double kEpsAccept = 0.05;
  constexpr double kTauRef3km_m = 0.0; // filled from the table below

  constexpr double kEmissions_km[] = {0.5, 2.0, 2.9, 3.0, 3.1, 5.0, 10.0, 20.0};
  constexpr double kDistances_km[] = {1.0, 10.0, 50.0, 100.0, 200.0, 400.0};
  constexpr double kLambdas_nm[] = {300.0, 400.0, 500.0, 600.0};

  // ------------------------------------------------------------------ table
  struct Table {
    std::size_t nLambda{0};
    std::size_t nHeight{0};
    std::vector<double> lambda_nm;
    std::vector<double> height_m;
    std::vector<std::vector<double>> tau; // [ilambda][iheight], cumulative from 0 m
    std::string path;
    bool loaded() const { return nLambda > 0 && nHeight > 0; }

    /// Index of the wavelength row (grid is 180...700 nm in 5 nm steps).
    std::size_t rowOf(double const lambdaNm) const {
      std::size_t best = 0;
      double bestD = 1.0e30;
      for (std::size_t i = 0; i < lambda_nm.size(); ++i) {
        double const d = std::fabs(lambda_nm[i] - lambdaNm);
        if (d < bestD) {
          bestD = d;
          best = i;
        }
      }
      return best;
    }
    std::size_t slab(double const h_m) const {
      std::size_t ih = 0;
      while (ih + 1 < height_m.size() && height_m[ih + 1] <= h_m) ++ih;
      return ih;
    }
    double alpha(double const h_m, std::size_t const row) const {
      std::size_t const ih = slab(h_m);
      if (ih + 1 >= height_m.size()) return 0.0;
      return (tau[row][ih + 1] - tau[row][ih]) / (height_m[ih + 1] - height_m[ih]);
    }
    double tauAt(double const h_m, std::size_t const row) const {
      std::size_t const ih = slab(h_m);
      if (ih + 1 >= height_m.size()) return tau[row][ih];
      double const dh = height_m[ih + 1] - height_m[ih];
      double const f = (h_m - height_m[ih]) / dh;
      return tau[row][ih] * (1.0 - f) + tau[row][ih + 1] * f;
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
      t.lambda_nm.resize(nl);
      t.height_m.resize(nh);
      bool ok = true;
      for (std::size_t i = 0; i < nl && ok; ++i) ok = static_cast<bool>(in >> t.lambda_nm[i]);
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

  double refrIndex(double const h_m) {
    return 1.0 + kRefrN0 * std::exp(-std::max(0.0, h_m) / kRefrH_m);
  }
  double refrDnDr(double const h_m) {
    return -kRefrN0 / kRefrH_m * std::exp(-std::max(0.0, h_m) / kRefrH_m);
  }


  // ------------------------------------------------ reference ray integrator
  struct RefOut {
    bool ok{false};
    char const* why{""};
    double tau{0.0};
    double L_m{0.0};
    double arc_m{0.0};
    double hMin_m{0.0};
    double hEnd_m{0.0};
    double uxEnd{0.0};
    double uzEnd{0.0};
    double invariantDrift{0.0};
  };

  /// Stop criterion: reach a target altitude (table nodes, traced backwards from
  /// the detector) or reach the closest approach to a target point (grid cells,
  /// which handles same-altitude horizontal paths where an altitude stop returns
  /// immediately).
  enum class StopMode { AtAltitude, AtPoint };

  /// RK4 in a spherically symmetric medium; tau exact in altitude.
  RefOut traceRay(Table const& tab, std::size_t const row, double const xStart,
                  double const zStart, double const uxStart, double const uzStart,
                  StopMode const mode, double const xTarget, double const zTarget,
                  double const stopAltitude_m, double const ds) {
    RefOut out;
    double xp = xStart, zp = zStart, ux = uxStart, uz = uzStart;
    double r = std::hypot(xp, zp);
    out.hMin_m = r - kEarthRadius_m;
    double s = 0.0;
    bool const descending = (uz * (zp / r) + ux * (xp / r)) < 0.0;
    double const p0 = refrIndex(out.hMin_m) * r * std::abs(ux * (zp / r) - uz * (xp / r));
    double distPrev = std::hypot(xStart - xTarget, zStart - zTarget);

    auto deriv = [](double const x, double const z, double const vx, double const vz,
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

      double const dhStep = hNew - hPrev;
      double tauStep = 0.0;
      if (std::fabs(dhStep) > 1.0e-6) {
        tauStep =
            std::fabs(tab.tauAt(hNew, row) - tab.tauAt(hPrev, row)) * ds / std::fabs(dhStep);
      } else {
        tauStep = tab.alpha(0.5 * (hPrev + hNew), row) * ds;
      }
      out.tau += tauStep;

      double const cross = std::abs(ux * (zp / r) - uz * (xp / r));
      out.arc_m += cross / r * ds * kEarthRadius_m;
      s += ds;
      out.L_m = s;
      out.hEnd_m = hNew;
      out.uxEnd = ux;
      out.uzEnd = uz;

      if (mode == StopMode::AtPoint) {
        double const dist = std::hypot(xp - xTarget, zp - zTarget);
        if (dist < 1.0e-3 || (s > 0.0 && dist > distPrev)) {
          out.ok = true;
          break;
        }
        distPrev = dist;
        continue;
      }

      double const ur = ux * (xp / r) + uz * (zp / r);
      if ((descending && ur >= 1.0e-12) || (!descending && ur <= -1.0e-12)) {
        out.why = "ray turned before reaching the target altitude";
        break;
      }
      bool const crossed =
          (descending && hNew <= stopAltitude_m) || (!descending && hNew >= stopAltitude_m);
      bool const wasOnSide =
          descending ? (hPrev >= stopAltitude_m) : (hPrev <= stopAltitude_m);
      if (crossed && wasOnSide && s > 0.0) {
        double const keep =
            std::max(0.0, std::min(1.0, (hPrev - stopAltitude_m) / (hPrev - hNew)));
        out.tau -= tauStep;
        if (std::fabs(stopAltitude_m - hPrev) > 1.0e-6) {
          out.tau += std::fabs(tab.tauAt(stopAltitude_m, row) - tab.tauAt(hPrev, row)) *
                     (ds * keep) / std::fabs(stopAltitude_m - hPrev);
        } else {
          out.tau += tab.alpha(hPrev, row) * ds * keep;
        }
        out.L_m = s - ds * (1.0 - keep);
        out.hEnd_m = stopAltitude_m;
        out.ok = true;
        break;
      }
    }
    if (out.ok) {
      double const pEnd = refrIndex(out.hEnd_m) * (kEarthRadius_m + out.hEnd_m) *
                          std::abs(ux * (zp / r) - uz * (xp / r));
      out.invariantDrift = (p0 > 0.0) ? std::abs(pEnd - p0) / p0 : 0.0;
    }
    return out;
  }


  // ------------------------------------------------------- table prototype
  struct Block {
    std::vector<double> h_m;
    std::vector<double> theta_deg;
    std::vector<std::vector<double>> tau; // [ih][ith]
  };

  struct SliceTable {
    Block a; // h >= 3000 m, theta <= 90 deg
    Block b; // h <= 3000 m, theta >= 90 deg
    char const* name{""};
    double build_s{0.0};
    std::size_t nodes() const {
      return a.h_m.size() * a.theta_deg.size() + b.h_m.size() * b.theta_deg.size();
    }
    std::size_t bytes() const { return nodes() * sizeof(double); }

    bool lookup(double const h_m, double const theta_deg, double& tauOut) const {
      Block const* blk = nullptr;
      if (theta_deg <= kZenithPointing_deg && h_m >= kDetectorAltitude_m) {
        blk = &a;
      } else if (theta_deg >= kZenithPointing_deg && h_m <= kDetectorAltitude_m) {
        blk = &b;
      } else {
        return false; // outside the physically valid L-shape
      }
      auto const& H = blk->h_m;
      auto const& T = blk->theta_deg;
      auto const& V = blk->tau;
      if (H.empty() || T.empty()) return false;
      double const h = std::max(H.front(), std::min(H.back(), h_m));
      double const th = std::max(T.front(), std::min(T.back(), theta_deg));
      std::size_t ih = 0;
      while (ih + 2 < H.size() && H[ih + 1] < h) ++ih;
      std::size_t it = 0;
      while (it + 2 < T.size() && T[it + 1] < th) ++it;
      double const fh = (H[ih + 1] > H[ih]) ? (h - H[ih]) / (H[ih + 1] - H[ih]) : 0.0;
      double const ft = (T[it + 1] > T[it]) ? (th - T[it]) / (T[it + 1] - T[it]) : 0.0;
      // NaN-aware bilinear: nodes outside the physically valid domain (a grazing
      // backward ray cannot reach far below/above the detector) are skipped.
      double wsum = 0.0;
      double tsum = 0.0;
      auto add = [&](double const v, double const w) {
        if (std::isfinite(v) && w > 0.0) {
          wsum += w;
          tsum += w * v;
        }
      };
      add(V[ih][it], (1.0 - fh) * (1.0 - ft));
      add(V[ih][it + 1], (1.0 - fh) * ft);
      add(V[ih + 1][it], fh * (1.0 - ft));
      add(V[ih + 1][it + 1], fh * ft);
      if (wsum <= 0.0) return false;
      tauOut = tsum / wsum;
      return true;
    }

    /// Number of nodes outside the valid domain (NaN by construction).
    std::size_t emptyNodes() const {
      std::size_t n = 0;
      for (auto const& blk : {a, b}) {
        for (auto const& row : blk.tau) {
          for (double const v : row) {
            if (!std::isfinite(v)) ++n;
          }
        }
      }
      return n;
    }
  };

  /// Node heights: block A above the detector, block B below, excluding the
  /// degenerate rows at exactly 3000 m (tau = 0 there by construction).
  std::vector<double> nodeHeightsA() { return {3100.0, 3500.0, 5000.0, 8000.0, 12000.0, 20000.0}; }
  std::vector<double> nodeHeightsB() { return {500.0, 1000.0, 2000.0, 2500.0, 2900.0}; }

  /// theta nodes: uniform in theta (as specified) or uniform in cos(theta)
  /// (step-5 candidate B, kept as a labelled diagnostic).
  std::vector<double> nodeThetas(double const lo_deg, double const hi_deg, int const n,
                                 bool const cosSpaced) {
    std::vector<double> t;
    t.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      double const f = (n > 1) ? static_cast<double>(i) / static_cast<double>(n - 1) : 0.0;
      if (cosSpaced) {
        double const c0 = std::cos(lo_deg * M_PI / 180.0);
        double const c1 = std::cos(hi_deg * M_PI / 180.0);
        t.push_back(std::acos(c0 + f * (c1 - c0)) * 180.0 / M_PI);
      } else {
        t.push_back(lo_deg + f * (hi_deg - lo_deg));
      }
    }
    return t;
  }

  /// Build one slice at 3 km: node value = exact-in-altitude tau of the ray that
  /// arrives at the detector with that zenith, traced backwards from the detector.
  void fillBlock(Table const& tab, std::size_t const row, Block& blk) {
    blk.tau.assign(blk.h_m.size(), std::vector<double>(blk.theta_deg.size(), 0.0));
    for (std::size_t it = 0; it < blk.theta_deg.size(); ++it) {
      double const th = blk.theta_deg[it] * M_PI / 180.0;
      double const ux = std::sin(th);
      double const uz = std::cos(th); // backward direction from the detector
      for (std::size_t ih = 0; ih < blk.h_m.size(); ++ih) {
        RefOut const r = traceRay(tab, row, 0.0, kEarthRadius_m + kDetectorAltitude_m, ux, uz,
                                  StopMode::AtAltitude, 0.0, 0.0, blk.h_m[ih], kStepRef_m);
        blk.tau[ih][it] = r.ok ? r.tau : std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  SliceTable buildSlice(Table const& tab, std::size_t const row, int const nThetaA,
                        int const nThetaB, bool const cosSpaced, char const* name) {
    SliceTable s;
    s.name = name;
    s.a.h_m = nodeHeightsA();
    s.b.h_m = nodeHeightsB();
    s.a.theta_deg = nodeThetas(kZenithPointing_deg - kFovHalfAngle_deg, kZenithPointing_deg,
                               nThetaA, cosSpaced);
    // theta = 90 exactly is invalid in block B: a horizontal backward ray from the
    // detector can never reach a lower altitude, so that node column would be NaN.
    s.b.theta_deg = nodeThetas(kZenithPointing_deg + 1.0e-3,
                               kZenithPointing_deg + kFovHalfAngle_deg, nThetaB, cosSpaced);
    auto const t0 = std::chrono::steady_clock::now();
    fillBlock(tab, row, s.a);
    fillBlock(tab, row, s.b);
    s.build_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return s;
  }


} // namespace

TEST_CASE("Cherenkov absorption table at 3 km, 90 deg pointing", "modules,cherenkov") {
  logging::set_level(logging::level::warn);

  Table const tab = loadTable();
  REQUIRE(tab.loaded());
  std::size_t const row400 = tab.rowOf(400.0);
  std::printf("\n[CherenkovTable3km] table %s; row(400 nm) = %zu, lambda nodes %zu\n",
              tab.path.c_str(), row400, tab.lambda_nm.size());
  std::printf("  detector 3000 m; pointing 90 deg zenith; FoV 10 deg FULL -> accepted arrival "
              "zenith %.1f ... %.1f deg\n",
              kZenithPointing_deg - kFovHalfAngle_deg,
              kZenithPointing_deg + kFovHalfAngle_deg);
  std::printf("  tau(400 nm, 3000 m) = %.6f; n(h) = 1 + %.3e exp(-h/%.0f m)\n",
              tab.tauAt(3000.0, row400), kRefrN0, kRefrH_m);

  // ------------------------------------------------------------ truth grid
  struct Cell {
    double h_em{0.0};
    double d_km{0.0};
    double theta_arr{0.0};
    double L_m{0.0};
    std::array<double, 4> tau{};   // per wavelength
    std::array<double, 4> tRef{};  // per wavelength
    bool reachable{false};
    bool inFov{false};
    char const* status{""};
  };

  std::vector<Cell> cells;
  int droppedUnreachable = 0;
  int droppedFov = 0;
  double maxHalving = 0.0, halvingAt[2] = {0.0, 0.0};
  double maxInvariant = 0.0;
  int nonFinite = 0;

  for (double const hEm_km : kEmissions_km) {
    for (double const d_km : kDistances_km) {
      double const hEm = hEm_km * 1000.0;
      double const d = d_km * 1000.0;
      // emission point at ground distance d from the detector (detector at the origin)
      double const dPhi = d / kEarthRadius_m;
      double const rEm = kEarthRadius_m + hEm;
      double const xEm = rEm * std::sin(dPhi);
      double const zEm = rEm * std::cos(dPhi);
      double const xDet = 0.0;
      double const zDet = kEarthRadius_m + kDetectorAltitude_m;
      double const chord = std::hypot(xDet - xEm, zDet - zEm);
      double const ux = (xDet - xEm) / chord;
      double const uz = (zDet - zEm) / chord;

      Cell c;
      c.h_em = hEm;
      c.d_km = d_km;
      RefOut const geo = traceRay(tab, row400, xEm, zEm, ux, uz, StopMode::AtPoint, xDet,
                                  zDet, kDetectorAltitude_m, kStepRef_m);
      if (!geo.ok) {
        ++droppedUnreachable;
        c.status = geo.why;
        cells.push_back(c);
        continue;
      }
      // arrival zenith: photon travels along (uxEnd, uzEnd); source direction is -(that)
      double const cosArr = std::max(-1.0, std::min(1.0, -geo.uzEnd));
      c.theta_arr = std::acos(cosArr) * 180.0 / M_PI;
      c.L_m = geo.L_m;
      c.reachable = true;
      c.inFov = std::fabs(c.theta_arr - kZenithPointing_deg) <= kFovHalfAngle_deg;
      c.status = c.inFov ? "in FoV" : "outside FoV";
      // convergence re-check on this grid
      RefOut const half = traceRay(tab, row400, xEm, zEm, ux, uz, StopMode::AtPoint, xDet, zDet,
                                   kDetectorAltitude_m, 0.5 * kStepRef_m);
      if (half.ok && geo.tau > 1.0e-9) {
        double const hv = std::fabs(geo.tau - half.tau) / geo.tau;
        if (hv > maxHalving) {
          maxHalving = hv;
          halvingAt[0] = hEm;
          halvingAt[1] = d_km;
        }
      }
      maxInvariant = std::max(maxInvariant, geo.invariantDrift);
      for (std::size_t iw = 0; iw < 4; ++iw) {
        RefOut const r =
            (iw == 0) ? geo
                      : traceRay(tab, tab.rowOf(kLambdas_nm[iw]), xEm, zEm, ux, uz,
                                 StopMode::AtPoint, xDet, zDet, kDetectorAltitude_m, kStepRef_m);
        c.tau[iw] = r.tau;
        c.tRef[iw] = std::exp(-r.tau);
        if (!std::isfinite(c.tau[iw]) || !std::isfinite(c.tRef[iw])) ++nonFinite;
      }
      if (!c.inFov) ++droppedFov;
      cells.push_back(c);
    }
  }

  int reachableInFov = 0;
  for (auto const& c : cells) {
    if (c.reachable && c.inFov) ++reachableInFov;
  }
  std::printf("\n  geometry: %zu (h_em, d) pairs -> reachable %d, dropped unreachable %d, "
              "outside FoV %d, in FoV %d (x 4 lambda = %d cells)\n",
              cells.size(), static_cast<int>(cells.size()) - droppedUnreachable,
              droppedUnreachable, droppedFov, reachableInFov, reachableInFov * 4);
  std::printf("  reference convergence on this grid: max |dtau(1 m) - dtau(0.5 m)|/dtau = "
              "%.3e at h_em = %.0f m, d = %.0f km\n",
              maxHalving, halvingAt[0], halvingAt[1]);
  std::printf("  Bouguer invariant drift: max %.3e; non-finite truth values: %d\n",
              maxInvariant, nonFinite);

  std::printf("\n  dropped cells (h_em, d, status):\n");
  for (auto const& c : cells) {
    if (!c.reachable || !c.inFov) {
      std::printf("    %6.0f m %6.0f km  arrival %s  %s\n", c.h_em, c.d_km,
                  c.reachable ? "out-FoV" : "-", c.status);
    }
  }


  // -------------------------------------------------- table comparison (400 nm)
  struct TableStats {
    std::string name;
    std::size_t nodes{0};
    std::size_t bytes{0};
    double build_s{0.0};
    int pass{0}, flag{0}, fail{0}, lookupFail{0};
    double worstEps{0.0}, worstH{0.0}, worstD{0.0}, worstTheta{0.0};
    int horizontalPass{0}, horizontalTotal{0};
    double lookup_ns{0.0};
  };

  std::vector<SliceTable> tables;
  tables.push_back(buildSlice(tab, row400, 5, 5, false, "coarse (5 theta/block)"));
  tables.push_back(buildSlice(tab, row400, 13, 13, false, "fine (13 theta/block)"));
  tables.push_back(buildSlice(tab, row400, 13, 13, true, "cos-spaced (13/block, diag)"));

  std::printf("\n  table prototypes (slice at 3000 m, lambda = 400 nm, no +-25 m window):\n");
  std::vector<TableStats> stats;
  for (auto const& t : tables) {
    TableStats st;
    st.name = t.name;
    st.nodes = t.nodes();
    st.bytes = t.bytes();
    st.build_s = t.build_s;
    std::size_t const nBench = 500000;
    volatile double sink = 0.0;
    auto const tb0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < nBench; ++i) {
      double tau = 0.0;
      t.lookup(10000.0, 89.5, tau);
      sink = sink + tau;
    }
    auto const tb1 = std::chrono::steady_clock::now();
    st.lookup_ns = std::chrono::duration<double, std::nano>(tb1 - tb0).count() /
                   static_cast<double>(nBench);

    for (auto const& c : cells) {
      if (!c.reachable || !c.inFov) continue;
      double tauTable = 0.0;
      if (!t.lookup(c.h_em, c.theta_arr, tauTable)) {
        ++st.lookupFail;
        continue;
      }
      double const tProd = std::exp(-tauTable);
      double const tRef = c.tRef[1]; // 400 nm
      double const eps =
          (tRef >= kEpsPass) ? std::fabs(tProd - tRef) / tRef : std::fabs(tProd - tRef);
      bool const horizontal = (c.h_em <= 3100.0);
      if (horizontal) ++st.horizontalTotal;
      if (eps <= kEpsPass) {
        ++st.pass;
        if (horizontal) ++st.horizontalPass;
      } else if (eps <= kEpsAccept) {
        ++st.flag;
      } else {
        ++st.fail;
      }
      if (eps > st.worstEps) {
        st.worstEps = eps;
        st.worstH = c.h_em;
        st.worstD = c.d_km;
        st.worstTheta = c.theta_arr;
      }
    }
    stats.push_back(st);
    std::printf("  %-32s nodes %4zu  %5zu B  build %.2f s  lookup %.1f ns  empty nodes %zu\n",
                st.name.c_str(), st.nodes, st.bytes, st.build_s, st.lookup_ns, t.emptyNodes());
    std::printf("      pass %d, flagged %d, fail %d, lookup-fail %d (of %d in-FoV cells)\n",
                st.pass, st.flag, st.fail, st.lookupFail, reachableInFov);
    std::printf("      worst eps = %.3e at h_em = %.0f m, d = %.0f km, arrival %.3f deg; "
                "horizontal band %d/%d within 1 %%\n",
                st.worstEps, st.worstH, st.worstD, st.worstTheta, st.horizontalPass,
                st.horizontalTotal);
  }


  // per-cell detail for the finest theta-uniform table
  std::printf("\n  fine table, per in-FoV cell (400 nm):\n");
  std::printf("  %7s %7s %8s %9s %12s %12s %9s %s\n", "h_em", "d_km", "arr_deg", "L_m",
              "T_table", "T_ref", "eps", "class");
  for (auto const& c : cells) {
    if (!c.reachable || !c.inFov) continue;
    double tauTable = 0.0;
    if (!tables[1].lookup(c.h_em, c.theta_arr, tauTable)) continue;
    double const tProd = std::exp(-tauTable);
    double const tRef = c.tRef[1];
    double const eps =
        (tRef >= kEpsPass) ? std::fabs(tProd - tRef) / tRef : std::fabs(tProd - tRef);
    std::printf("  %7.0f %7.0f %8.3f %9.1f %12.6f %12.6f %9.2e %s\n", c.h_em, c.d_km,
                c.theta_arr, c.L_m, tProd, tRef, eps,
                (eps <= kEpsPass) ? "pass" : ((eps <= kEpsAccept) ? "flag" : "FAIL"));
  }

  // -------------------------------------------------- d = 0 clean-reject probe
  std::printf("\n  d = 0 probe (not part of the grid):\n");
  double probeTheta = -1.0;
  {
    double const hEm = 10000.0;
    RefOut const r =
        traceRay(tab, row400, 0.0, kEarthRadius_m + hEm, 0.0, -1.0, StopMode::AtPoint, 0.0,
                 kEarthRadius_m + kDetectorAltitude_m, kDetectorAltitude_m, kStepRef_m);
    probeTheta = r.ok ? std::acos(std::max(-1.0, std::min(1.0, -r.uzEnd))) * 180.0 / M_PI : -1.0;
    std::printf("    h_em = %.0f m, d = 0, reachable = %d, arrival %.2f deg, tau = %.6f -> %s\n",
                hEm, static_cast<int>(r.ok), probeTheta, r.tau,
                (std::fabs(probeTheta - kZenithPointing_deg) > kFovHalfAngle_deg)
                    ? "rejected: outside FoV"
                    : "in FoV");
    CHECK(r.ok);
    CHECK(std::isfinite(r.tau));
  }

  // ------------------------------------------------------------- acceptance
  TableStats const& best = *std::max_element(
      stats.begin(), stats.end(), [](TableStats const& x, TableStats const& y) {
        return (x.pass - x.fail) < (y.pass - y.fail);
      });
  double const passFraction =
      (reachableInFov > 0) ? static_cast<double>(best.pass) / reachableInFov : 0.0;
  std::printf("\n[CherenkovTable3km] acceptance (best table: %s)\n", best.name.c_str());
  std::printf("  1) >= 99 %% of reachable cells within 1 %%: %.1f %% (%d/%d) -> %s\n",
              100.0 * passFraction, best.pass, reachableInFov,
              (passFraction >= 0.99) ? "PASS" : "FAIL");
  std::printf("  2) no cell above 5 %%: worst eps = %.3e -> %s\n", best.worstEps,
              (best.worstEps <= kEpsAccept) ? "PASS" : "FAIL");
  std::printf("  3) lookup >= 100x faster than the tracer (254 us): %.1f ns = %.0fx -> %s\n",
              best.lookup_ns, 254000.0 / best.lookup_ns,
              (best.lookup_ns <= 2540.0) ? "PASS" : "FAIL");
  std::printf("  4) horizontal/upward cells covered (tracer returns T = 1 there): %d/%d "
              "within 1 %% -> %s\n",
              best.horizontalPass, best.horizontalTotal,
              (best.horizontalTotal > 0 && best.horizontalPass == best.horizontalTotal)
                  ? "PASS"
                  : "FAIL");
  std::printf("  5) d = 0 rejected with an explicit status and no NaN: arrival %.2f deg -> "
              "%s\n",
              probeTheta,
              (std::fabs(probeTheta - kZenithPointing_deg) > kFovHalfAngle_deg) ? "PASS"
                                                                               : "FAIL");
  std::printf("  table size: %.1f kB per lambda slice (%.1f kB for 4 lambdas)\n",
              static_cast<double>(best.bytes) / 1024.0,
              4.0 * static_cast<double>(best.bytes) / 1024.0);
  std::printf("  theta lever arm at h_em = h_det: dtheta/dd = 1/(2R) = %.3e rad/m -> one "
              "0.1 deg step = %.1f km of path\n",
              1.0 / (2.0 * kEarthRadius_m),
              0.1 * M_PI / 180.0 * 2.0 * kEarthRadius_m / 1000.0);

  CHECK(nonFinite == 0);
  CHECK(passFraction >= 0.99);
  CHECK(best.worstEps <= kEpsAccept);
  CHECK(best.lookup_ns <= 2540.0);
  CHECK(best.horizontalTotal > 0);
  CHECK(best.horizontalPass == best.horizontalTotal);
  CHECK(std::fabs(probeTheta - kZenithPointing_deg) > kFovHalfAngle_deg);
}

