/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * \file testCherenkovChordRoute.cpp
 *
 * Step 1 (2g-h2g): deterministic chord-route transmission against an
 * independent, test-only spherical reference over the full geometry grid, plus
 * one Monte-Carlo run each at Z30 and Z60 to confirm that the sampling path is
 * wired to the same chord model.
 *
 * PRODUCTION UNDER TEST (read-only):
 *   AtmosphericAbsorption::getTransmission  (AtmosphericAbsorption.hpp:171-196)
 *     tau_slant = |tau(l,h_emit) - tau(l,h_det)| / max(|cos_zenith|, 0.01)   (:181-192)
 *     T         = exp(-tau_slant)                                             (:195)
 *   fed with cos_zenith = |dz|/L of the PHOTON chord (Cherenkov.inl:2017-2022).
 *
 * REFERENCE (this file only, no production code):
 *   straight chord through a spherical Earth, R_E = 6371 km,
 *     r(s)^2 = r0^2 - 2 s r0 cos(theta) + s^2 ,  r0 = R_E + h_emit ,
 *     h(s)   = sqrt(r(s)^2) - R_E ,
 *     tau    = int alpha(h) ds ,  alpha = d tau/dh from the same MODTRAN table
 *              (400 nm row, cumulative from the ground up),
 *   integrated with steps <= 1 m that are split at the table's 1 km knots, so
 *   the quadrature is exact for the table's piecewise linear cumulative tau.
 *   The 1 m -> 0.5 m halving check is reported for every cell.
 *
 * PRE-REGISTERED THRESHOLDS: per photon |dT| <= 1e-3 and |dtau|/tau <= 1e-3.
 */

#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/Cherenkov.hpp>

#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/StraightTrajectory.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <catch2/catch_all.hpp>

#include <PhysicalUnitsCatch2.hpp>
#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>
#include <SetupTestTrajectory.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace corsika;

namespace {

  // ------------------------------------------------------------- constants
  constexpr double kEarthRadius_m = 6.371e6;
  constexpr double kLambda_nm = 400.0;
  constexpr double kDetectorAltitude_m = 2000.0; // base detector of the grid
  constexpr double kTolAbsT = 1.0e-3;            // |dT| per photon
  constexpr double kTolRelTau = 1.0e-3;          // |dtau|/tau per photon
  constexpr double kStep_m = 1.0;                // reference integration step
  constexpr double kKnot_m = 1000.0;             // table height spacing
  constexpr double kTauRefDet_m = 0.216105;      // tau(400 nm, 2000 m)
  constexpr std::size_t kBunchesMC = 1600000;    // 3 sigma < 1e-3 at T ~ 0.8
  constexpr std::uint64_t kCascadeSeed = 20261114;
  constexpr double kRefractiveIndex = 1.0003;
  constexpr double kLambdaMin_nm = 399.0;
  constexpr double kLambdaMax_nm = 401.0;

  // ------------------------------------------------------------------ table
  struct Table {
    std::size_t nLambda{0};
    std::size_t nHeight{0};
    std::vector<double> lambda_nm;
    std::vector<double> height_m;
    std::vector<std::vector<double>> tau; // [ilambda][iheight], cumulative from 0 m
    std::string path;
    bool loaded() const { return nLambda > 0 && nHeight > 0; }

    std::size_t slab(double const h_m) const {
      std::size_t ih = 0;
      while (ih + 1 < height_m.size() && height_m[ih + 1] <= h_m) ++ih;
      return ih;
    }
    /// Piecewise-constant extinction of the 400 nm row.
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

  constexpr double kBunchSize = 1.0;
  constexpr double kMCWindow_m = 200.0;          // +-200 m around the centre
  constexpr double kMCEmission_m = 10000.0;

  constexpr double kZeniths[] = {0.0, 30.0, 60.0, 75.0, 85.0, 88.0,
                                 89.0, 90.0, 91.0, 95.0, 120.0, 180.0};
  constexpr double kEmissions[] = {500.0, 2000.0, 5000.0, 10000.0, 20000.0};

  // ------------------------------------------------------------- geometry
  // Straight chord: r(s)^2 = r0^2 - 2 s r0 c + s^2, c = cos(theta), h = r - R_E.
  // Returns s at the first forward crossing of altitude hTarget, or -1 when the
  // chord never reaches that altitude (closest approach above it).
  double sAtAltitude(double const hStart_m, double const cosT, double const hTarget_m) {
    double const r0 = kEarthRadius_m + hStart_m;
    double const rT = kEarthRadius_m + hTarget_m;
    double const disc = rT * rT - r0 * r0 * (1.0 - cosT * cosT);
    if (disc < 0.0) return -1.0;
    double const sq = std::sqrt(disc);
    double const s = (cosT > 0.0) ? (r0 * cosT - sq) : (r0 * cosT + sq);
    return (s >= 0.0) ? s : -1.0;
  }

  /// Closest approach altitude of the forward chord (turning point).
  double closestAltitude(double const hStart_m, double const cosT) {
    double const r0 = kEarthRadius_m + hStart_m;
    return r0 * std::sqrt(std::max(0.0, 1.0 - cosT * cosT)) - kEarthRadius_m;
  }

  /// Minimum altitude along the segment [0, L] of the chord (not the global
  /// parabola minimum, which may lie beyond the segment).
  double pathMinAltitude(double const hStart_m, double const cosT, double const L_m) {
    double const r0 = kEarthRadius_m + hStart_m;
    double const sTurn = r0 * cosT;
    if (sTurn > 0.0 && sTurn < L_m) {
      return r0 * std::sqrt(std::max(0.0, 1.0 - cosT * cosT)) - kEarthRadius_m;
    }
    double const rEnd = std::sqrt(r0 * r0 - 2.0 * L_m * r0 * cosT + L_m * L_m);
    return std::min(hStart_m, rEnd - kEarthRadius_m);
  }

  struct RefIntegral {
    double tau{0.0};
    double hMin_m{0.0};
    double maxStep_m{0.0};
    std::size_t steps{0};
  };

  /// tau = int_0^L alpha(h(s)) ds with steps <= dsMax, split at the 1 km knots so
  /// that the quadrature is exact for the table's piecewise constant alpha.
  RefIntegral integrateChord(Table const& tab, double const hStart_m, double const cosT,
                             double const L_m, double const dsMax) {
    RefIntegral out;
    double const r0 = kEarthRadius_m + hStart_m;
    out.hMin_m = 1.0e30;
    double s = 0.0;
    for (std::size_t iter = 0; iter < 400000000ULL; ++iter) {
      if (!(s < L_m)) break;
      double const rs = std::sqrt(r0 * r0 - 2.0 * s * r0 * cosT + s * s);
      double const h = rs - kEarthRadius_m;
      out.hMin_m = std::min(out.hMin_m, h);
      double nextS = L_m;
      double const knot =
          (cosT > 0.0) ? std::floor(h / kKnot_m) * kKnot_m : std::ceil(h / kKnot_m) * kKnot_m;
      if (knot != h) {
        double const sk = sAtAltitude(hStart_m, cosT, knot);
        if (sk > s + 1.0e-9) nextS = std::min(L_m, sk);
      }
      double step = std::min(dsMax, nextS - s);
      if (!(step > 0.0)) break;
      double const sMid = s + 0.5 * step;
      double const rMid = std::sqrt(r0 * r0 - 2.0 * sMid * r0 * cosT + sMid * sMid);
      out.tau += tab.alpha(rMid - kEarthRadius_m) * step;
      out.maxStep_m = std::max(out.maxStep_m, step);
      s += step;
      ++out.steps;
    }
    if (L_m <= 0.0) out.hMin_m = hStart_m;
    return out;
  }

  // --------------------------------------------------------- MC cone helpers
  double cherenkovCosThetaC() {
    double const massGeV = get_mass(Code::Electron) / 1_GeV;
    double const E = 1.0;
    double const beta = std::sqrt(1.0 - (massGeV * massGeV) / (E * E));
    return 1.0 / (beta * kRefractiveIndex);
  }

  /// |dir_z| of a photon on the cone of a track at zenith thetaT.
  double dirZ(double const thetaT, double const cosC, double const sinC, double const phi) {
    return std::cos(thetaT) * cosC + std::sin(thetaT) * sinC * std::cos(phi);
  }

  /// Cone footprint radius (max distance from the mean-ray hit point) at the
  /// observation level, over the emission window.
  double coneFootprintRadius(double const thetaT, double const hCentre_m) {
    double const cosC = cherenkovCosThetaC();
    double const sinC = std::sqrt(std::max(0.0, 1.0 - cosC * cosC));
    double const tanT = std::tan(thetaT);
    double const shiftX_m = (hCentre_m - kDetectorAltitude_m) * tanT;
    double rMax = 0.0;
    for (int ih = 0; ih <= 32; ++ih) {
      double const h = hCentre_m - kMCWindow_m + 2.0 * kMCWindow_m * ih / 32.0;
      double const xEmit = (hCentre_m - h) * tanT;
      for (int ip = 0; ip < 360; ++ip) {
        double const phi = 2.0 * M_PI * ip / 360.0;
        double const dz = dirZ(thetaT, cosC, sinC, phi);
        double const dxDir = std::sin(thetaT) * cosC - std::cos(thetaT) * sinC * std::cos(phi);
        double const dyDir = -sinC * std::sin(phi);
        double const t = (h - kDetectorAltitude_m) / dz;
        rMax = std::max(rMax, std::hypot(xEmit + t * dxDir - shiftX_m, t * dyDir));
      }
    }
    return rMax;
  }

  /// Chord-model prediction averaged over the cone and the emission window
  /// (this is what the process implements with curvature off).
  double chordPrediction(Table const& tab, double const thetaT, double const hCentre_m) {
    double const cosC = cherenkovCosThetaC();
    double const sinC = std::sqrt(std::max(0.0, 1.0 - cosC * cosC));
    double const tauDet = tab.tauAt(kDetectorAltitude_m);
    double sum = 0.0;
    int n = 0;
    for (int ih = 0; ih <= 64; ++ih) {
      double const h = hCentre_m - kMCWindow_m + 2.0 * kMCWindow_m * ih / 64.0;
      double const dTau = std::abs(tab.tauAt(h) - tauDet);
      for (int ip = 0; ip < 360; ++ip) {
        double const phi = 2.0 * M_PI * ip / 360.0;
        double const dz = std::abs(dirZ(thetaT, cosC, sinC, phi));
        sum += std::exp(-dTau / std::max(dz, 0.01));
        ++n;
      }
    }
    return sum / n;
  }

  /// Spherical reference averaged over the same cone and emission window, on a
  /// coarse (h, phi) grid so that the 1 m integration stays affordable.
  double sphericalPrediction(Table const& tab, double const thetaT, double const hCentre_m) {
    double const cosC = cherenkovCosThetaC();
    double const sinC = std::sqrt(std::max(0.0, 1.0 - cosC * cosC));
    int const nh = 8;
    int const nphi = 16;
    double sum = 0.0;
    int n = 0;
    for (int ih = 0; ih <= nh; ++ih) {
      double const h = hCentre_m - kMCWindow_m + 2.0 * kMCWindow_m * ih / nh;
      for (int ip = 0; ip < nphi; ++ip) {
        double const phi = 2.0 * M_PI * (ip + 0.5) / nphi;
        double const dz = dirZ(thetaT, cosC, sinC, phi);
        double const L = sAtAltitude(h, dz, kDetectorAltitude_m);
        if (L < 0.0) continue;
        RefIntegral const r = integrateChord(tab, h, dz, L, kStep_m);
        sum += std::exp(-r.tau);
        ++n;
      }
    }
    return (n > 0) ? sum / n : 0.0;
  }


  // ------------------------------------------------------------------- MC run
  struct MCRun {
    std::size_t emitted{0};
    std::size_t recorded{0};
    std::size_t outside{0};
    double T{0.0};
    double wall_s{0.0};
    double radiusMean_m{0.0};
  };

  /// One Monte-Carlo run of the production chord route (curvature off), with the
  /// telescope radius sized from the cone footprint plus 20 %.
  MCRun runMC(double const zenithDeg, double const radius_m, std::string const& absFile) {
    logging::set_level(logging::level::warn);
    double const theta = zenithDeg * M_PI / 180.0;
    double const sinT = std::sin(theta);
    double const cosT = std::cos(theta);
    double const tanT = (cosT != 0.0) ? sinT / cosT : 0.0;
    double const earthRadiusM = kEarthRadius_m;
    double const groundZ_m = earthRadiusM + kDetectorAltitude_m;    // arg 5
    double const telescopeZ_m = groundZ_m;                          // arg 14
    double const shift_m = (kMCEmission_m - kDetectorAltitude_m) * tanT;
    double const stepLength_m = 2.0 * kMCWindow_m / cosT;
    double const xStart_m = -kMCWindow_m * tanT;
    double const zStart_m = earthRadiusM + kMCEmission_m + kMCWindow_m;

    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Electron);
    auto const& cs = *csPtr;
    [[maybe_unused]] auto const& env_dummy = env;
    auto [stack, viewPtr] = setup::testing::setup_stack(Code::Electron, 1_GeV, nodePtr, cs);
    [[maybe_unused]] test::StackView& view = *viewPtr;
    auto particle = stack->getNextParticle();

    Point const start(cs, xStart_m * 1_m, 0_m, zStart_m * 1_m);
    particle.setPosition(start);
    particle.setDirection(DirectionVector(cs, {sinT, 0., -cosT}));
    VelocityVector const velocity(cs, {sinT * constants::c, 0_m / second, -cosT * constants::c});
    Line const line(start, velocity);
    StraightTrajectory const track =
        setup::testing::make_track<StraightTrajectory>(line, stepLength_m * 1_m / constants::c);
    Step step(particle, track);

    std::size_t outside = 0;
    double radiusSum = 0.0;
    std::size_t calls = 0;
    // weighted transmission accumulators (deterministic absorption, see out.T below)
    double weightIn = 0.0, weightOutside = 0.0;
    std::size_t weightBunches = 0;
    cherenkov::CherenkovProcess::HitRecorder const recorder =
        [&](double const x_m, double const y_m, double /*z_m*/, double /*lambda_nm*/,
            double /*z_prod_m*/, double /*dx*/, double /*dy*/, double /*dz*/, double /*gx*/,
            double /*gy*/, double /*gz*/, double /*time*/, double const weight,
            double /*n_prod*/, double /*path*/) {
          double const r = std::hypot(x_m - shift_m, y_m);
          radiusSum += r;
          ++weightBunches;
          weightIn += weight;
          if (r > radius_m) { ++outside; weightOutside += weight; }
          ++calls;
        };

    cherenkov::CherenkovProcess proc(
        1.0003 /*n*/, kLambdaMin_nm, kLambdaMax_nm, kBunchSize, groundZ_m, recorder,
        cherenkov::ProjectionMode::Ground,
        std::shared_ptr<media::ShowerAxis const>{}, 0.0, 0.0, 0.0, shift_m, 0.0, telescopeZ_m,
        radius_m, std::function<std::pair<double, double>(double, double)>{}, earthRadiusM,
        -1.0, -1.0, 0.0, 0.0, true /*enableDispersion*/, false /*enableCurvature*/, absFile,
        0.0, std::shared_ptr<cherenkov::RefractiveIndexTable const>{}, false /*curvedTable*/);
    proc.setProfiling(true);

    auto const t0 = std::chrono::steady_clock::now();
    while (proc.getProfile().bunches < kBunchesMC) {
      proc.doContinuous(step, true);
      ++calls;
    }
    // Absorption is DETERMINISTIC since 2g-h2n wave 4c: every bunch is recorded and its
    // weight is multiplied by exp(-tau), so a COUNT ratio is 1 by construction.  The
    // transmission this test must compare against the chord prediction is the WEIGHT
    // ratio (measured: count ratio gave |dT| = 0.24-0.37 at Z60, i.e. it looked like a
    // 24 % physics failure when it was a stale metric).
    double const weightEmitted = double(proc.getProfile().bunches) * kBunchSize;
    MCRun out;
    out.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.emitted = proc.getProfile().bunches;
    out.recorded = proc.getProfile().recorded;
    out.outside = outside;
    out.T = (weightEmitted > 0.0) ? (weightIn - weightOutside) / weightEmitted : 0.0;
    out.radiusMean_m = (calls > 0) ? radiusSum / static_cast<double>(calls) : 0.0;
    return out;
  }


  // Detector altitude for the grid: 2 km below the emission point for downward
  // tracks, the same altitude for the horizontal case, and the next mirrored
  // altitude above the emission point for upward tracks.
  double detectorFor(double const hEmit_m, double const zenithDeg) {
    if (zenithDeg < 90.0) return kDetectorAltitude_m;
    if (zenithDeg == 90.0) return hEmit_m;
    if (hEmit_m < kDetectorAltitude_m) return kDetectorAltitude_m; // 0.5 -> 2
    if (hEmit_m < 5000.0) return 5000.0;                           // 2   -> 5
    if (hEmit_m < 10000.0) return 10000.0;                         // 5   -> 10
    return 35000.0;                                                // 10, 20 -> 35
  }

  struct Worst {
    double value{0.0};
    double zenithDeg{0.0};
    double hEmit_m{0.0};
    double hDet_m{0.0};
  };
  void update(Worst& w, double const v, double const z, double const he, double const hd) {
    if (v > w.value) {
      w.value = v;
      w.zenithDeg = z;
      w.hEmit_m = he;
      w.hDet_m = hd;
    }
  }


} // namespace

TEST_CASE("Cherenkov chord route vs spherical reference", "modules,cherenkov") {
  logging::set_level(logging::level::warn);

  // Fixed seed (documented open item: the process-internal survival stream is
  // hard-coded at 98765, so runs are reproducible but not independent).
  auto& rmng = RNGManager<>::getInstance();
  rmng.setSeed(kCascadeSeed);
  rmng.registerRandomStream("cascade");

  Table const tab = loadTable();
  REQUIRE(tab.loaded());
  cherenkov::AtmosphericAbsorption abs;
  abs.load(tab.path);
  REQUIRE(abs.isLoaded());

  std::printf("\n[CherenkovChord] production table %s\n", tab.path.c_str());
  std::printf("  tau(400 nm, 2000 m) = %.6f (anchor %.6f); isLoaded = %d\n", tab.tauAt(2000.0),
              kTauRefDet_m, static_cast<int>(abs.isLoaded()));
  CHECK(tab.tauAt(2000.0) == Approx(kTauRefDet_m).margin(1.0e-6));

  // =========================================================== deterministic
  // Main grid of the checkpoint: 12 zeniths x 5 emission heights.
  // ===========================================================================
  int cells = 0;
  int defined = 0;
  int reachable = 0;
  int passed = 0;
  int notDefined = 0;
  int unreachable = 0;
  int nonFinite = 0;
  int degenerate = 0;
  double maxHalving = 0.0;
  std::size_t lookupCalls = 0;
  Worst worstT;
  Worst worstTau;
  double thetaMaxAccuracy = -1.0;
  double thetaMaxGeometry = -1.0;
  bool accuracyStillClean = true;
  bool geometryStillClean = true;

  std::printf("  %6s %8s %7s %10s %9s %11s %11s %10s %11s %9s %9s %s\n", "theta", "h_emit",
              "h_det", "L_m", "h_min_m", "T_prod", "T_ref", "|dT|", "|dtau|/tau", "steps",
              "halving", "status");
  for (double const zenithDeg : kZeniths) {
    double const theta = zenithDeg * M_PI / 180.0;
    double const c = std::cos(theta);
    bool accClean = true;
    bool geoClean = true;
    for (double const hEmit : kEmissions) {
      ++cells;
      std::printf("  %6.1f %8.0f", zenithDeg, hEmit);
      double const hDet = detectorFor(hEmit, zenithDeg);
      double const hDetPrint = hDet;
      if (zenithDeg < 90.0 && !(hEmit > hDet)) {
        ++notDefined;
        std::printf(" %7.0f %10s %9s %11s %11s %10s %11s %9s %9s n/a: emission not above "
                    "the detector\n",
                    hDetPrint, "-", "-", "-", "-", "-", "-", "-", "-");
        continue;
      }
      ++defined;
      double const L = sAtAltitude(hEmit, c, hDet);
      if (zenithDeg == 90.0) {
        ++degenerate;
        double const Tp90 = abs.getTransmission(kLambda_nm, hEmit, hDet, 0.0);
        std::printf(" %7.0f %10.1f %9.1f %11.6f %11s %11s %10s %11s  %s\n", hDetPrint, L,
                    pathMinAltitude(hEmit, c, L), Tp90, "-", "-", "-", "-",
                    "degenerate: L = 0, |dz| = 0, T = 1 by construction (see block 2)");
        continue;
      }
      double const hMin = pathMinAltitude(hEmit, c, L);
      if (L < 0.0) {
        ++unreachable;
        geoClean = false;
        double const Tp0 = abs.getTransmission(kLambda_nm, hEmit, hDet, c);
        std::printf(" %7.0f %10s %9.1f %11.4f %11s %11s %10s %11s  %s\n", hDetPrint, "none",
                    closestAltitude(hEmit, c), Tp0, "-", "-", "-", "-",
                    "unreachable: chord closest approach is above h_det");
        continue;
      }
      ++reachable;
      // production: exactly the call the process makes on the chord route
      double const Tp = abs.getTransmission(kLambda_nm, hEmit, hDet, c);
      ++lookupCalls;

      RefIntegral const r1 = integrateChord(tab, hEmit, c, L, kStep_m);
      RefIntegral const r2 = integrateChord(tab, hEmit, c, L, 0.5 * kStep_m);
      double const Tr = std::exp(-r1.tau);
      double const dT = std::fabs(Tp - Tr);
      double const tauProd =
          (Tp > 0.0) ? -std::log(Tp) : std::numeric_limits<double>::infinity();
      double const relTau =
          (r1.tau > 1.0e-9) ? std::fabs(tauProd - r1.tau) / r1.tau : 0.0;
      double const halving = std::fabs(r1.tau - r2.tau) / std::max(r1.tau, 1.0e-30);
      maxHalving = std::max(maxHalving, halving);
      if (!std::isfinite(Tp) || !std::isfinite(Tr) || !std::isfinite(dT) ||
          !std::isfinite(relTau)) {
        ++nonFinite;
      }
      bool const pass = (dT <= kTolAbsT) && (relTau <= kTolRelTau);
      if (pass) {
        ++passed;
      } else {
        accClean = false;
      }
      update(worstT, dT, zenithDeg, hEmit, hDet);
      update(worstTau, relTau, zenithDeg, hEmit, hDet);
      std::printf(" %7.0f %10.1f %9.1f %11.3e %11.6f %10.3e %11.3e %9zu %9.1e %s\n", hDetPrint,
                  L, hMin, Tp, Tr, dT, relTau, r1.steps, halving, pass ? "pass" : "FAIL");
    }
    if (accuracyStillClean && accClean) {
      thetaMaxAccuracy = zenithDeg;
    } else {
      accuracyStillClean = false;
    }
    if (geometryStillClean && geoClean) {
      thetaMaxGeometry = zenithDeg;
    } else {
      geometryStillClean = false;
    }
  }


  // ------------------------------------------- (2) grazing horizontal chords
  // theta = 90 deg exactly: |dz| = 0 and h_emit = h_det, so the chord route has
  // no way to express the finite path.  Cold numbers for a 20 km horizontal leg.
  std::printf("\n[CherenkovChord] (2) grazing horizontal chords: theta = 90 deg, L = 20 km\n");
  std::printf("  %8s %11s %11s %10s\n", "h_m", "T_prod", "T_ref", "|dT|");
  double worstGrazing = 0.0;
  for (double const hEmit : kEmissions) {
    double const L = 20000.0;
    double const Tp = abs.getTransmission(kLambda_nm, hEmit, hEmit, 0.0);
    RefIntegral const r = integrateChord(tab, hEmit, 0.0, L, kStep_m);
    double const Tr = std::exp(-r.tau);
    worstGrazing = std::max(worstGrazing, std::fabs(Tp - Tr));
    std::printf("  %8.0f %11.6f %11.6f %10.3e\n", hEmit, Tp, Tr, std::fabs(Tp - Tr));
    // the chord route cannot leave T = 1 here (|dtau| = 0 for h_emit = h_det)
    CHECK(Tp == 1.0);
  }
  std::printf("  worst |dT| at 90 deg (20 km leg) = %.3e\n", worstGrazing);


  // ------------------------------------------------- (3) MC wiring checks
  // One run each, N = 1.6e6, radius = cone footprint + 20 %.  The gate is
  // against the deterministic chord prediction (the model the process implements),
  // because this run only has to show that the sampling is wired to it.
  std::printf("\n[CherenkovChord] (3) MC wiring check (N >= %zu, radius = footprint x 1.2)\n",
              kBunchesMC);
  std::printf("  %5s %9s %9s %11s %11s %11s %10s %10s %9s %8s\n", "case", "foot_r_m",
              "radius_m", "T_MC", "T_chord", "T_sph", "|dT|_chord", "sigma", "gate", "wall");
  double worstMC = 0.0;
  for (double const zenithDeg : {30.0, 60.0}) {
    double const theta = zenithDeg * M_PI / 180.0;
    double const foot = coneFootprintRadius(theta, kMCEmission_m);
    double const radius = 1.2 * foot;
    MCRun const mc = runMC(zenithDeg, radius, tab.path);
    double const tChord = chordPrediction(tab, theta, kMCEmission_m);
    double const tSph = sphericalPrediction(tab, theta, kMCEmission_m);
    double const sigma = std::sqrt(std::max(0.0, mc.T * (1.0 - mc.T)) /
                                   static_cast<double>(std::max<std::size_t>(1, mc.recorded)));
    double const gate = std::max(kTolAbsT, 3.0 * sigma);
    double const dChord = std::fabs(mc.T - tChord);
    double const dSph = std::fabs(mc.T - tSph);
    worstMC = std::max(worstMC, dChord);
    std::printf("  Z%-4.0f %9.2f %9.2f %11.6f %11.6f %11.6f %10.3e %10.3e %9.3e %8.3f s\n",
                zenithDeg, foot, radius, mc.T, tChord, tSph, dChord, sigma, gate, mc.wall_s);
    std::printf("        emitted %zu, recorded %zu, outside radius %zu; |T_MC - T_sph| = "
                "%.3e; N = %zu\n",
                mc.emitted, mc.recorded, mc.outside, dSph, mc.recorded);
    CHECK(mc.emitted >= kBunchesMC);
    CHECK(mc.outside == 0);
    CHECK(dChord <= gate);
  }

  // ------------------------------------------------------------- verdict
  // t_lookup: repeated calls at one representative geometry (clock overhead of a
  // single call is ~1000 ns and would swamp the number).
  volatile double sink = 0.0;
  std::size_t const kBenchCalls = 200000;
  auto const tb0 = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < kBenchCalls; ++i) {
    sink = sink + abs.getTransmission(kLambda_nm, 10000.0, 2000.0, 0.8660254037844386);
  }
  auto const tb1 = std::chrono::steady_clock::now();
  double const nsPerLookup =
      std::chrono::duration<double, std::nano>(tb1 - tb0).count() /
      static_cast<double>(kBenchCalls);
  std::printf("\n[CherenkovChord] verdict (thresholds: |dT| <= %.0e, |dtau|/tau <= %.0e)\n",
              kTolAbsT, kTolRelTau);
  std::printf("  cells %d | defined %d | reachable %d | degenerate %d | n/a %d | unreachable "
              "%d\n",
              cells, defined, reachable, degenerate, notDefined, unreachable);
  std::printf("  pass %d / %d reachable (%.1f %%)\n", passed, reachable,
              (reachable > 0) ? 100.0 * passed / reachable : 0.0);
  std::printf("  theta_max (accuracy, contiguous from 0 deg) = %.0f deg\n", thetaMaxAccuracy);
  std::printf("  theta_max (geometry, contiguous from 0 deg) = %.0f deg\n", thetaMaxGeometry);
  std::printf("  worst |dT|        = %.3e at theta = %.0f deg, h_emit = %.0f m, h_det = %.0f "
              "m\n",
              worstT.value, worstT.zenithDeg, worstT.hEmit_m, worstT.hDet_m);
  std::printf("  worst |dtau|/tau  = %.3e at theta = %.0f deg, h_emit = %.0f m, h_det = %.0f "
              "m\n",
              worstTau.value, worstTau.zenithDeg, worstTau.hEmit_m, worstTau.hDet_m);
  std::printf("  non-finite cells  = %d\n", nonFinite);
  std::printf("  reference convergence: max |dtau(1 m) - dtau(0.5 m)| / dtau = %.3e\n",
              maxHalving);
  std::printf("  t_lookup (getTransmission) = %.1f ns per call (%zu repeated calls, one "
              "geometry)\n",
              nsPerLookup, kBenchCalls);
  std::printf("  reference table anchors: tau(400,2000) = %.6f, tau(400,5000) = %.6f\n",
              tab.tauAt(2000.0), tab.tauAt(5000.0));

  CHECK(nonFinite == 0);
  CHECK(maxHalving < 1.0e-6);
  CHECK(lookupCalls == static_cast<std::size_t>(reachable));
}

