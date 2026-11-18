/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * \file testCherenkovE2E.cpp
 *
 * ON/OFF measurement of the production atmospheric transmission on the
 * chord route (retained L3 coarse check, 0.1 % target).  Numbers only.
 *
 * Fixed settings for every case: 1 GeV electron; lambda = 399 ... 401 nm;
 * bunchSize = 1; constant n = 1.000300 (arg 16 empty); dispersion flag true
 * (constant-index production branch); curvature off (args 23/27 = false);
 * projection mode Ground; telescope radius 500 m; absolute Earth-centred z
 * (z_abs = EarthRadius + altitude, as the applications use it); absorption
 * OFF and ON differ only in constructor argument 24.
 *
 * Cases (all with a +-200 m production window and 200 000 bunches per run):
 *   V10  vertical, production 10 km   (checkpoint 2g-h2e case, must reproduce)
 *   V5   vertical, production 5 km    (second emission height)
 *   Z30  track zenith 30 deg, 10 km   (1/cos factor)
 *   Z60  track zenith 60 deg, 10 km   (1/cos factor)
 *
 * ANGLE USED BY THE PRODUCTION CODE (checked before any run):
 *   Cherenkov.inl:2017-2022 forms the chord (det - emission) and sets
 *   cos_zenith = |dz_path| / path_length.  In Ground mode the hit is
 *   det = posEmit + t*d (Cherenkov.inl:1664-1666), so the chord is parallel to
 *   the PHOTON direction, not to the track.  Hence
 *       theta_used = theta_photon,  |dir_z| = cos(th_t)cos(th_C)
 *                                             + sin(th_t)sin(th_C)cos(phi)
 *   with th_C = arccos(1/(beta n)) the Cherenkov half-angle and phi the azimuth
 *   on the cone.  That value enters Cherenkov.inl:2029-2031 ->
 *   absorptionSurvives (Cherenkov.inl:745-747) ->
 *   AtmosphericAbsorption::survives -> getTransmission
 *   (AtmosphericAbsorption.hpp:171-196), which returns
 *       T = exp( -|tau(h_emit) - tau(h_det)| / |cos_zenith| ).
 *   The reference below therefore averages exp(-dtau/|dir_z|) over the cone
 *   azimuth phi and over the uniform emission window, and also quotes the
 *   single-angle anchor exp(-dtau/cos th_track).
 */

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
#include <string>
#include <utility>
#include <vector>

using namespace corsika;

namespace {

  // -------------------------------------------------------------- fixed case
  constexpr double kEarthRadius_m = 6.371e6;         // arg 17 (constructor default)
  constexpr double kObservationAltitude_m = 2000.0;  // observation level
  constexpr double kHalfWindow_m = 200.0;            // production window +-200 m
  constexpr double kTelescopeRadius_m = 500.0;       // arg 15
  constexpr double kLambdaMin_nm = 399.0;            // arg 2
  constexpr double kLambdaMax_nm = 401.0;            // arg 3
  constexpr double kBunchSize = 1.0;                 // arg 4
  constexpr double kRefractiveIndex = 1.0003;        // arg 1
  constexpr double kPrimaryEnergy_GeV = 1.0;
  constexpr std::uint64_t kCascadeSeed = 20261114;   // Poisson / wavelength stream
  constexpr std::size_t kBunchesPerRun = 200000;

  // Pre-registered reference anchors, tau at 400 nm (table md5
  // dfb2b7fd0d0c57ee0a704035a6a81e0b).
  constexpr double kTauRefDet_m = 0.216105;  // tau(400 nm, 2000 m)
  constexpr double kTauRefV10_m = 0.449487;  // tau(400 nm, 10000 m)
  constexpr double kTauRefV5_m = 0.331942;   // tau(400 nm,  5000 m)
  constexpr double kTRefV10Anchor = 0.791851;
  constexpr double kTauNodeTol = 1.0e-6;

  // Checkpoint 2g-h2e vertical result, reproduced by this build: the printed
  // T_prod values and the exact bunch counts of that run.
  // RE-FROZEN for 2g-h2n wave 4c (deterministic absorption).  The previous values were
  // bunch-COUNT ratios: 158579/200115 = 0.792439, i.e. the number of bunches that survived
  // a Bernoulli draw -- a SAMPLE with sigma = sqrt(T(1-T)/N) = 9.1e-4.  Absorption now
  // multiplies the recorded weight by exp(-tau) and drops nothing, so the count ratio is
  // 1 by construction and the transmission is the WEIGHT ratio, which is the exact
  // expectation rather than a sample:
  //
  //   previous (sampled): 0.792439 / 0.792386 / 0.792434   accepted 158579 / 158585 / 158593
  //   now (exact):        0.791840 / 0.791840 / 0.791839   accepted 200115 / 200136 / 200134
  //
  // The 6e-4 shift is 0.66 sigma of the OLD estimator -- the two agree within its own
  // noise, which is exactly what an unbiased estimator change must do.  The physics gate
  // below (|T_prod - T_ref| <= 1e-3 against the truth) is untouched and is the real test.
  constexpr double kPrevV10[3] = {0.791840, 0.791840, 0.791839};
  constexpr std::size_t kPrevV10Emitted[3] = {200115, 200136, 200134};
  constexpr std::size_t kPrevV10Accepted[3] = {200115, 200136, 200134};
  constexpr double kReproduceTol = 1.0e-6; // printed to 6 decimals, weight-derived

  // Pre-registered gates.
  constexpr double kGateAbsolute = 1.0e-3;   // |T_prod - T_ref|
  constexpr double kGateRelative = 5.0e-3;   // 0.5 % wherever T_ref > 0.05
  constexpr double kOffAcceptanceMin = 0.999;
  constexpr double kStopPointSigmaFactor = 5.0;

  /// One measurement case: a straight track of given zenith angle whose
  /// production window is centred on a given emission altitude.
  struct CaseSpec {
    char const* name;
    double zenithDeg;
    double centreAltitude_m;
  };


  // ---------------------------------------------------------------------------
  // Minimal independent reader for the shipped MODTRAN table.  Used only to
  // extract the reference values tau(400 nm, 2000 m) and tau(400 nm, 10000 m);
  // both are exact grid nodes.  Format: AtmosphericAbsorption.hpp:13-21,
  // 105 wavelengths x 51 heights, cumulative tau from the ground upwards.
  // ---------------------------------------------------------------------------
  struct ReferenceTable {
    std::size_t nLambda{0};
    std::size_t nHeight{0};
    std::vector<double> lambda_nm;
    std::vector<double> height_m;
    std::vector<std::vector<double>> tau; // [ilambda][iheight]
    std::string path;

    bool loaded() const { return nLambda > 0 && nHeight > 0; }

    /// Bilinear interpolation in the same order as
    /// AtmosphericAbsorption::getOpticalDepth (height first, then wavelength).
    double tauAt(double const lambdaNm, double const h_m) const {
      double const lam = std::max(lambda_nm.front(), std::min(lambda_nm.back(), lambdaNm));
      double const h = std::max(height_m.front(), std::min(height_m.back(), h_m));
      std::size_t iw = 0;
      for (std::size_t i = 0; i + 1 < lambda_nm.size(); ++i) {
        iw = i;
        if (lambda_nm[i + 1] >= lam) break;
      }
      std::size_t ih = 0;
      for (std::size_t i = 0; i + 1 < height_m.size(); ++i) {
        ih = i;
        if (height_m[i + 1] >= h) break;
      }
      double const fw = (lambda_nm[iw + 1] != lambda_nm[iw])
                            ? (lam - lambda_nm[iw]) / (lambda_nm[iw + 1] - lambda_nm[iw])
                            : 0.0;
      double const fh = (height_m[ih + 1] != height_m[ih])
                            ? (h - height_m[ih]) / (height_m[ih + 1] - height_m[ih])
                            : 0.0;
      double const tau_w0 = tau[iw][ih] * (1.0 - fh) + tau[iw][ih + 1] * fh;
      double const tau_w1 = tau[iw + 1][ih] * (1.0 - fh) + tau[iw + 1][ih + 1] * fh;
      return tau_w0 * (1.0 - fw) + tau_w1 * fw;
    }
  };

  ReferenceTable loadReferenceTable() {
    std::vector<std::string> candidates;
    // Absolute fallback (verified in this workspace).
    candidates.push_back(
        "/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/"
        "abstable_MODTRAN_bottomup.dat");
    // Source-relative: <repo>/tests/modules/testCherenkovE2E.cpp ->
    // <repo>/modules/data/CHERENKOV/atmosphere/...
    std::string const here(__FILE__);
    auto const marker = here.rfind("/tests/modules/");
    if (marker != std::string::npos) {
      candidates.push_back(here.substr(0, marker) +
                           "/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat");
    }
    // Working-directory relatives (ctest runs from the build tree).
    candidates.push_back("modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat");
    candidates.push_back("../modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat");
    candidates.push_back("../../modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat");

    for (auto const& path : candidates) {
      std::ifstream in(path);
      if (!in.is_open()) continue;
      std::size_t nl = 0;
      std::size_t nh = 0;
      if (!(in >> nl >> nh) || nl == 0 || nh == 0) continue;
      ReferenceTable r;
      r.nLambda = nl;
      r.nHeight = nh;
      r.lambda_nm.resize(nl);
      r.height_m.resize(nh);
      bool ok = true;
      for (std::size_t i = 0; i < nl && ok; ++i) ok = static_cast<bool>(in >> r.lambda_nm[i]);
      for (std::size_t i = 0; i < nh && ok; ++i) ok = static_cast<bool>(in >> r.height_m[i]);
      r.tau.assign(nl, std::vector<double>(nh, 0.0));
      for (std::size_t i = 0; i < nl && ok; ++i) {
        for (std::size_t j = 0; j < nh && ok; ++j) ok = static_cast<bool>(in >> r.tau[i][j]);
      }
      if (!ok) continue;
      r.path = path;
      return r;
    }
    return ReferenceTable{};
  }

  // ---------------------------------------------------------------------------
  // Cherenkov half-angle exactly as the process computes it (Cherenkov.inl:1323-1324,
  // :1589): cos(theta_C) = 1/(beta*n) for the primary energy of this case.
  // ---------------------------------------------------------------------------
  double cherenkovCosThetaC() {
    double const massGeV = get_mass(Code::Electron) / 1_GeV;
    double const beta =
        std::sqrt(1.0 - (massGeV * massGeV) / (kPrimaryEnergy_GeV * kPrimaryEnergy_GeV));
    return 1.0 / (beta * kRefractiveIndex);
  }

  /// Pre-registered reference for one case, all from tau(400 nm, h) of the table.
  struct CaseReference {
    double deltaTau{0.0};            // tau(h_centre) - tau(2000 m)
    double anchor{0.0};             // exp(-deltaTau / cos(theta_track))
    double coneAtCentre{0.0};       // <exp(-deltaTau / |dir_z|)> over the cone
    double heightCorrected{0.0};    // <exp(-[tau(h) - tau(2000)] / |dir_z|)> over h and phi
    double predictedAcceptance{0.0}; // cone fraction within radius 500 m of the shift
    double predictedPath_m{0.0};     // mean chord length
  };

  CaseReference referenceFor(ReferenceTable const& ref, CaseSpec const& spec) {
    constexpr int kPhiSamples = 720;
    constexpr int kHeightSamples = 128;
    double const theta = spec.zenithDeg * M_PI / 180.0;
    double const sinT = std::sin(theta);
    double const cosT = std::cos(theta);
    double const tanT = sinT / cosT;
    double const cosC = cherenkovCosThetaC();
    double const sinC = std::sqrt(std::max(0.0, 1.0 - cosC * cosC));
    double const tauDet = ref.tauAt(400.0, kObservationAltitude_m);
    double const tauCentre = ref.tauAt(400.0, spec.centreAltitude_m);
    // Telescope centre: where the mean ray from the window centre lands.
    double const shiftX_m = (spec.centreAltitude_m - kObservationAltitude_m) * tanT;

    CaseReference out;
    out.deltaTau = tauCentre - tauDet;
    out.anchor = std::exp(-out.deltaTau / cosT);

    double sumCone = 0.0;
    for (int i = 0; i < kPhiSamples; ++i) {
      double const phi = 2.0 * M_PI * (i + 0.5) / kPhiSamples;
      double const dz = cosT * cosC + sinT * sinC * std::cos(phi); // |dir_z|
      sumCone += std::exp(-out.deltaTau / dz);
    }
    out.coneAtCentre = sumCone / kPhiSamples;

    double sumWindow = 0.0;
    double sumPath = 0.0;
    double sumInside = 0.0;
    for (int ih = 0; ih < kHeightSamples; ++ih) {
      double const h = spec.centreAltitude_m - kHalfWindow_m +
                       2.0 * kHalfWindow_m * (ih + 0.5) / kHeightSamples;
      // x of the emission point on the track, 0 at h = centreAltitude.
      double const xEmit = (spec.centreAltitude_m - h) * tanT;
      double const dTau = ref.tauAt(400.0, h) - tauDet;
      for (int ip = 0; ip < kPhiSamples; ++ip) {
        double const phi = 2.0 * M_PI * (ip + 0.5) / kPhiSamples;
        double const dz = cosT * cosC + sinT * sinC * std::cos(phi);
        double const dxDir = sinT * cosC - cosT * sinC * std::cos(phi);
        double const dyDir = -sinC * std::sin(phi);
        double const t = (h - kObservationAltitude_m) / dz;
        double const xHit = xEmit + t * dxDir;
        double const yHit = t * dyDir;
        sumWindow += std::exp(-dTau / dz);
        sumPath += t;
        sumInside += (std::hypot(xHit - shiftX_m, yHit) <= kTelescopeRadius_m) ? 1.0 : 0.0;
      }
    }
    double const nSamples = static_cast<double>(kPhiSamples) * kHeightSamples;
    out.heightCorrected = sumWindow / nSamples;
    out.predictedAcceptance = sumInside / nSamples;
    out.predictedPath_m = sumPath / nSamples;
    return out;
  }

  // ---------------------------------------------------------------------------
  // One measurement: emitted bunches, recorded (= accepted) bunches and the
  // bookkeeping needed to audit the geometry.
  // ---------------------------------------------------------------------------
  struct RunResult {
    std::string label;
    bool absorptionOn{false};
    std::size_t emitted{0};      // bunches generated      (Cherenkov.inl:1568)
    std::size_t recorded{0};     // bunches written out    (Cherenkov.inl:2037-2047)
    std::size_t recorderCalls{0};
    std::size_t outsideRadius{0};
    double radiusSum_m{0.0};
    double weightSum{0.0};  ///< sum of recorded bunch weights (photons)
    double radiusMax_m{0.0};
    std::size_t zProdCount{0};
    double zProdMin_m{0.0};
    double zProdMax_m{0.0};
    double zProdSum_m{0.0};
    double pathMin_m{0.0};
    double pathMax_m{0.0};
    double pathSum_m{0.0};
    double lambdaMin_nm{0.0};
    double lambdaMax_nm{0.0};
    double nProdMin{0.0};
    double nProdMax{0.0};
    std::size_t calls{0};
    std::size_t absorptionDraws{0}; // CherenkovProcess::absorptionDrawCount()
    double cosZenithSum{0.0};       // sum of |dir_ground_z| over recorded bunches
    double wall_s{0.0};
    std::array<std::size_t, 11> radiusHist{}; // 50 m bins, [10] = above 500 m

    double acceptance() const {
      return (emitted > 0) ? static_cast<double>(insideRadius()) / static_cast<double>(emitted)
                           : 0.0;
    }
    /// Bunches that both survived (or, OFF, were never absorbed) and landed
    /// inside the telescope radius.  Ground mode performs no radius cut of its
    /// own (Cherenkov.inl:1629-1673), so the cut is applied by the recorder.
    std::size_t insideRadius() const { return recorded - outsideRadius; }
    double zProdMean_m() const {
      return (zProdCount > 0) ? zProdSum_m / static_cast<double>(zProdCount) : 0.0;
    }
    double pathMean_m() const {
      return (zProdCount > 0) ? pathSum_m / static_cast<double>(zProdCount) : 0.0;
    }
    double radiusMean_m() const {
      return (recorderCalls > 0) ? radiusSum_m / static_cast<double>(recorderCalls) : 0.0;
    }
    double cosZenithMean() const {
      return (recorderCalls > 0) ? cosZenithSum / static_cast<double>(recorderCalls) : 0.0;
    }
    /// Weighted transmission: sum of recorded weights / sum of emitted weights.
    /// Absorption is DETERMINISTIC since 2g-h2n wave 4c -- every bunch is recorded and
    /// its weight carries exp(-tau) -- so a COUNT ratio is 1 by construction (measured:
    /// 200115/200115 = 1.0 against the frozen 0.7924, i.e. |dT| = 0.2076).
    double transmission() const {
      double const wEmitted = static_cast<double>(emitted) * kBunchSize;
      return (wEmitted > 0.0) ? weightSum / wEmitted : acceptance();
    }
  };

  /// Build the process for one case, emit until kBunchesPerRun bunches have been
  /// generated and report the numbers for this run.  OFF and ON differ only in
  /// `absorptionFile` (constructor argument 24).
  RunResult runCase(CaseSpec const& spec, bool const absorptionOn,
                    std::string const& absorptionFile, std::string const& label,
                    bool const printGeometry) {
    logging::set_level(logging::level::warn);

    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Electron);
    auto const& cs = *csPtr;
    [[maybe_unused]] auto const& env_dummy = env;
    [[maybe_unused]] auto const& node_dummy = nodePtr;

    auto [stack, viewPtr] = setup::testing::setup_stack(Code::Electron, 1_GeV, nodePtr, cs);
    [[maybe_unused]] test::StackView& view = *viewPtr;
    auto particle = stack->getNextParticle();

    // ---- geometry: one production window along a straight inclined track ----
    double const theta = spec.zenithDeg * M_PI / 180.0;
    double const sinTheta = std::sin(theta);
    double const cosTheta = std::cos(theta);
    double const tanTheta = (cosTheta != 0.0) ? sinTheta / cosTheta : 0.0;

    double const earthRadiusM = kEarthRadius_m;                        // arg 17
    double const groundZ_m = earthRadiusM + kObservationAltitude_m;    // arg  5 (absolute)
    double const telescopeZ_m = earthRadiusM + kObservationAltitude_m; // arg 14 (absolute)
    double const telescopeRadius_m = kTelescopeRadius_m;               // arg 15

    // Telescope centre: where the mean ray from the window centre lands at the
    // observation level.  Zero (exactly) for a vertical track.
    double const telescopeShift_m =
        (tanTheta != 0.0) ? (spec.centreAltitude_m - kObservationAltitude_m) * tanTheta : 0.0;
    double const telescopeX_m = telescopeShift_m;                      // arg 12
    double const telescopeY_m = 0.0;                                   // arg 13

    // Emission window: centre altitude +- kHalfWindow m along the track, with
    // x = 0 at the window centre.
    double const stepLength_m = 2.0 * kHalfWindow_m / cosTheta;
    double const xStart_m = (tanTheta != 0.0) ? -kHalfWindow_m * tanTheta : 0.0;
    double const zStart_m = earthRadiusM + spec.centreAltitude_m + kHalfWindow_m;
    double const zEnd_m = earthRadiusM + spec.centreAltitude_m - kHalfWindow_m;

    Point const start(cs, xStart_m * 1_m, 0_m, zStart_m * 1_m);
    particle.setPosition(start);
    DirectionVector const direction(cs, {sinTheta, 0., -cosTheta});
    particle.setDirection(direction);

    VelocityVector const velocity(
        cs, {sinTheta * constants::c, 0_m / second, -cosTheta * constants::c});
    Line const line(start, velocity);
    StraightTrajectory const track = setup::testing::make_track<StraightTrajectory>(
        line, stepLength_m * 1_m / constants::c);
    Step step(particle, track);


    RunResult res;
    res.label = label;
    res.absorptionOn = absorptionOn;

    // ---- the 27 constructor arguments, one named local each ----------------
    //  Cherenkov.hpp:240-266, in that order.  Arguments 9-15 and 17-21 are all
    //  plain doubles, so they are listed one per line to make a swap impossible
    //  to miss.
    double const refractiveIndex = kRefractiveIndex;               //  1
    double const lambdaMin = kLambdaMin_nm;                        //  2
    double const lambdaMax = kLambdaMax_nm;                        //  3
    double const bunchSize = kBunchSize;                           //  4
    //  groundZ_m                                                  //  5 (above)
    //  recorder                                                   //  6 (below)
    cherenkov::ProjectionMode const mode = cherenkov::ProjectionMode::Ground; // 7
    std::shared_ptr<media::ShowerAxis const> const primaryShowerAxis = nullptr; // 8
    double const showerCoreX_m = 0.0;                              //  9
    double const showerCoreY_m = 0.0;                              // 10
    double const showerCoreZ_m = 0.0;                              // 11
    //  telescopeX_m, telescopeY_m, telescopeZ_m                  // 12-14 (above)
    //  telescopeRadius_m                                          // 15 (above)
    std::function<std::pair<double, double>(double, double)> const refractiveIndexFunc =
        {};                                                        // 16
    //  earthRadiusM                                               // 17 (above)
    double const telescopeAzimuth_rad = -1.0;                      // 18
    double const telescopeZenith_rad = -1.0;                       // 19
    double const showerAzimuth_rad = 0.0;                          // 20
    double const showerZenith_rad = 0.0;                           // 21
    bool const enableDispersion = true;                            // 22
    bool const enableCurvature = false;                            // 23
    //  absorptionFile                                             // 24 (parameter)
    double const telescopePositionCorrectionFraction = 0.0;         // 25
    std::shared_ptr<cherenkov::RefractiveIndexTable const> const
        refractiveIndexTable = nullptr;                            // 26
    bool const enableCurvedTable = false;                          // 27

    // ---- recorder (argument 6): counts, radius, heights, path lengths ------
    //  HitRecorder signature: x, y, z, lambda, z_prod, (dir x3), (dir_ground
    //  x3), time, weight, n_prod, path_length  -- 15 doubles.
    //  The radius is measured from the shifted telescope centre; the chord
    //  length is rebuilt from the photon's own |dir_z|, exactly as the process
    //  does at Cherenkov.inl:2017-2022.
    double radiusSum = 0.0;
    double radiusMax = 0.0;
    double zSum = 0.0;
    double zMin = 1.0e30;
    double zMax = -1.0e30;
    double pathSum = 0.0;
    double pathMin = 1.0e30;
    double pathMax = -1.0e30;
    double cosZenithSum = 0.0;
    double lamMinSeen = 1.0e30;
    double lamMaxSeen = -1.0e30;
    double nMinSeen = 1.0e30;
    double nMaxSeen = -1.0e30;
    std::size_t recorderCalls = 0;
    double weightSum = 0.0;  // photons represented by the recorded bunches
    std::size_t outsideRadius = 0;
    std::array<std::size_t, 11> radiusHist{};

    cherenkov::CherenkovProcess::HitRecorder const recorder =
        [&](double const x_m, double const y_m, double const z_m, double const lambda_nm,
            double const z_prod_m, double /*dir_x*/, double /*dir_y*/, double /*dir_z*/,
            double /*dir_ground_x*/, double /*dir_ground_y*/, double const dir_ground_z,
            double /*time_ns*/, double const weight, double const n_prod,
            double /*path_length_m*/) {
          ++recorderCalls;
          weightSum += weight;
          double const r = std::hypot(x_m - telescopeX_m, y_m - telescopeY_m);
          radiusSum += r;
          radiusMax = std::max(radiusMax, r);
          if (r > telescopeRadius_m) ++outsideRadius;
          radiusHist[std::min<std::size_t>(10, static_cast<std::size_t>(r / 50.0))] += 1;
          zMin = std::min(zMin, z_prod_m);
          zMax = std::max(zMax, z_prod_m);
          zSum += z_prod_m;
          double const cosZenith = std::fabs(dir_ground_z);
          double const path = (cosZenith > 1.0e-12)
                                  ? (z_prod_m - (z_m - earthRadiusM)) / cosZenith
                                  : 0.0;
          cosZenithSum += cosZenith;
          pathMin = std::min(pathMin, path);
          pathMax = std::max(pathMax, path);
          pathSum += path;
          lamMinSeen = std::min(lamMinSeen, lambda_nm);
          lamMaxSeen = std::max(lamMaxSeen, lambda_nm);
          nMinSeen = std::min(nMinSeen, n_prod);
          nMaxSeen = std::max(nMaxSeen, n_prod);
          (void)weight;
        };

    cherenkov::CherenkovProcess proc(refractiveIndex, lambdaMin, lambdaMax, bunchSize,
                                     groundZ_m, recorder, mode, primaryShowerAxis,
                                     showerCoreX_m, showerCoreY_m, showerCoreZ_m,
                                     telescopeX_m, telescopeY_m, telescopeZ_m,
                                     telescopeRadius_m, refractiveIndexFunc, earthRadiusM,
                                     telescopeAzimuth_rad, telescopeZenith_rad,
                                     showerAzimuth_rad, showerZenith_rad, enableDispersion,
                                     enableCurvature, absorptionFile,
                                     telescopePositionCorrectionFraction, refractiveIndexTable,
                                     enableCurvedTable);
    proc.setProfiling(true); // bunches emitted / recorded (Cherenkov.inl:1568, :2047)

    if (printGeometry) {
      std::printf("\n[CherenkovE2E] %s: geometry and medium (identical for OFF and ON; only "
                  "argument 24 differs)\n",
                  spec.name);
      std::printf("  track               zenith %.2f deg, azimuth 0; direction (%.6f, 0, %.6f)\n",
                  spec.zenithDeg, sinTheta, -cosTheta);
      std::printf("  observation level   %.1f m -> groundZ_m = %.6f (arg 5), telescopeZ_m = "
                  "%.6f (arg 14)\n",
                  kObservationAltitude_m, groundZ_m, telescopeZ_m);
      std::printf("  production window   z = %.6f ... %.6f (abs) = %.1f ... %.1f m altitude, "
                  "x = %.3f ... %.3f m, step %.3f m along track\n",
                  zStart_m, zEnd_m, spec.centreAltitude_m + kHalfWindow_m,
                  spec.centreAltitude_m - kHalfWindow_m, xStart_m,
                  (tanTheta != 0.0) ? kHalfWindow_m * tanTheta : 0.0, stepLength_m);
      std::printf("  telescope           (%.3f, %.1f, %.6f) abs, radius %.1f m; shift from "
                  "(0,0) = %.3f m; mode Ground\n",
                  telescopeX_m, telescopeY_m, telescopeZ_m, telescopeRadius_m,
                  telescopeShift_m);
      std::printf("                      tele az/zenith = %.1f/%.1f rad, shower az/zenith = "
                  "%.1f/%.1f rad (unused in Ground mode)\n",
                  telescopeAzimuth_rad, telescopeZenith_rad, showerAzimuth_rad,
                  showerZenith_rad);
      std::printf("  refractive index    n = %.6f constant, no z/lambda dependence "
                  "(arg 16 empty), n_prod printed per run\n",
                  refractiveIndex);
      std::printf("  medium              lambda %.1f ... %.1f nm, bunchSize %.1f, "
                  "enableDispersion = %s\n",
                  lambdaMin, lambdaMax, bunchSize, enableDispersion ? "true" : "false");
      std::printf("                      -> constant-index production branch "
                  "(Cherenkov.inl:1409-1416): dispersion has no effect here\n");
      std::printf("  curvature           enableCurvature = %s, enableCurvedTable = %s "
                  "(production default); absorption path = chord model\n",
                  enableCurvature ? "true" : "false", enableCurvedTable ? "true" : "false");
      std::printf("                      (Cherenkov.inl:739-747 -> "
                  "AtmosphericAbsorption::survives, T = exp(-|dtau|/cos_zenith))\n");
    }


    // ---- emit until the required number of bunches has been generated ------
    auto const t0 = std::chrono::steady_clock::now();
    while (proc.getProfile().bunches < kBunchesPerRun) {
      proc.doContinuous(step, true);
      ++res.calls;
    }
    res.wall_s =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    res.emitted = proc.getProfile().bunches;
    res.recorded = proc.getProfile().recorded;
    res.absorptionDraws = proc.absorptionDrawCount();
    res.recorderCalls = recorderCalls;
    res.weightSum = weightSum;
    res.outsideRadius = outsideRadius;
    res.radiusSum_m = radiusSum;
    res.radiusMax_m = radiusMax;
    res.zProdCount = recorderCalls;
    res.zProdMin_m = (recorderCalls > 0) ? zMin : 0.0;
    res.zProdMax_m = (recorderCalls > 0) ? zMax : 0.0;
    res.zProdSum_m = zSum;
    res.pathMin_m = (recorderCalls > 0) ? pathMin : 0.0;
    res.pathMax_m = (recorderCalls > 0) ? pathMax : 0.0;
    res.pathSum_m = pathSum;
    res.cosZenithSum = cosZenithSum;
    res.lambdaMin_nm = (recorderCalls > 0) ? lamMinSeen : 0.0;
    res.lambdaMax_nm = (recorderCalls > 0) ? lamMaxSeen : 0.0;
    res.nProdMin = (recorderCalls > 0) ? nMinSeen : 0.0;
    res.nProdMax = (recorderCalls > 0) ? nMaxSeen : 0.0;
    res.radiusHist = radiusHist;

    std::printf("  [%s/%s] absorptionFile = \"%s\"  (argument 24)\n", spec.name, label.c_str(),
                absorptionFile.c_str());
    std::printf("    emitted %zu bunches in %zu doContinuous calls, recorded %zu "
                "(acceptance %.4f %%), outside %.1f m: %zu\n",
                res.emitted, res.calls, res.recorded, 100.0 * res.acceptance(),
                telescopeRadius_m, res.outsideRadius);
    std::printf("    radius from shifted centre (%.3f, %.1f): mean %.3f m, max %.3f m\n",
                telescopeX_m, telescopeY_m, res.radiusMean_m(), res.radiusMax_m);
    std::printf("    emission height %.3f ... %.3f m (mean %.3f, n = %zu); chord length "
                "%.3f ... %.3f m (mean %.3f); |cos_zenith| mean %.6f\n",
                res.zProdMin_m, res.zProdMax_m, res.zProdMean_m(), res.zProdCount,
                res.pathMin_m, res.pathMax_m, res.pathMean_m(), res.cosZenithMean());

    std::printf("    lambda %.4f ... %.4f nm; n_prod %.6f ... %.6f\n", res.lambdaMin_nm,
                res.lambdaMax_nm, res.nProdMin, res.nProdMax);
    std::printf("    survival draws %zu (contract: exactly 1 per bunch, 0 when absorption is "
                "off)\n",
                res.absorptionDraws);
    std::printf("    wall %.3f s\n", res.wall_s);

    if (!printGeometry) return res; // histograms only for the control run

    std::printf("    distance-from-centre histogram (50 m bins):");
    for (std::size_t i = 0; i < res.radiusHist.size(); ++i) {
      if (i < 10) {
        std::printf(" %zu-%zu:%zu", i * 50, (i + 1) * 50, res.radiusHist[i]);
      } else {
        std::printf(" >500:%zu", res.radiusHist[i]);
      }
    }
    std::printf("\n");

    return res;
  }

} // namespace

TEST_CASE("Cherenkov ON/OFF atmospheric transmission (end to end)", "modules,cherenkov") {
  logging::set_level(logging::level::warn);

  // The Poisson photon count and the wavelengths come from the "cascade"
  // stream; seed it explicitly so the runs are reproducible.  The absorption
  // survival draw uses the process-internal stream (seed 98765,
  // Cherenkov.inl:71, taken once per bunch at Cherenkov.inl:745-747), which is
  // not reachable from outside: a hard-coded seed, logged as an open item.  The
  // stream is registered once, so every run consumes a disjoint part of the
  // same sequence and repeats are not independent.
  auto& rmng = RNGManager<>::getInstance();
  rmng.setSeed(kCascadeSeed);
  rmng.registerRandomStream("cascade");

  // ---------------------------------------------------------- reference data
  ReferenceTable const ref = loadReferenceTable();
  REQUIRE(ref.loaded());
  double const cosC = cherenkovCosThetaC();
  std::printf("\n[CherenkovE2E] reference table %s\n", ref.path.c_str());
  std::printf("  grid %zu wavelengths x %zu heights; tau(400 nm, 2000 m) = %.6f\n",
              ref.nLambda, ref.nHeight, ref.tauAt(400.0, kObservationAltitude_m));
  CHECK(ref.tauAt(400.0, kObservationAltitude_m) == Approx(kTauRefDet_m).margin(kTauNodeTol));
  CHECK(ref.tauAt(400.0, 10000.0) == Approx(kTauRefV10_m).margin(kTauNodeTol));
  CHECK(ref.tauAt(400.0, 5000.0) == Approx(kTauRefV5_m).margin(kTauNodeTol));

  // ------------------------------------------------------- pre-registration
  // All reference values are formed from the table BEFORE any emission.
  //   theta_used = theta_photon (see the file header): |dir_z| =
  //   cos(th_t)cos(th_C) + sin(th_t)sin(th_C)cos(phi), th_C = 1.402973 deg.
  //   anchor           = exp(-[tau(h_c) - tau(2000)] / cos(th_t))
  //   height-corrected = <exp(-[tau(h) - tau(2000)] / |dir_z|)> over h and phi
  std::vector<CaseSpec> const cases = {
      {"V10", 0.0, 10000.0}, {"V5", 0.0, 5000.0}, {"Z30", 30.0, 10000.0}, {"Z60", 60.0, 10000.0}};
  std::vector<CaseReference> refs;
  refs.reserve(cases.size());
  std::printf("\n[CherenkovE2E] PRE-REGISTERED REFERENCES (table md5 dfb2b7fd..., 400 nm, "
              "theta_C = %.6f deg)\n",
              std::acos(cosC) * 180.0 / M_PI);
  std::printf("  %-4s %7s %8s %10s %10s %10s %10s %10s\n", "case", "zenith", "h_centre",
              "dtau", "T_ref_anch", "T_ref_hc", "cone_cov", "pred_path");
  for (std::size_t i = 0; i < cases.size(); ++i) {
    refs.push_back(referenceFor(ref, cases[i]));
    CaseReference const& r = refs.back();
    std::printf("  %-4s %6.1f^ %7.0f m %10.6f %10.6f %10.6f %10.4f %10.1f m\n", cases[i].name,
                cases[i].zenithDeg, cases[i].centreAltitude_m, r.deltaTau, r.anchor,
                r.heightCorrected, r.predictedAcceptance, r.predictedPath_m);
  }


  // One case: OFF control, then (only if the control passes) one ON run.
  struct CaseOutcome {
    std::string name;
    CaseSpec spec{};
    CaseReference ref{};
    RunResult off;
    RunResult on;
    bool controlPassed{false};
    bool onRun{false};
    double tProd{0.0};
    double sigma{0.0};
    double dTAnchor{0.0};
    double dTHc{0.0};
    double rel{0.0};
    double gate{0.0};
    double stopPoint{0.0};
  };

  std::vector<CaseOutcome> outcomes;
  outcomes.reserve(cases.size());

  // ---------------------------------------------------- V10 (must reproduce)
  {
    CaseOutcome c;
    c.name = cases[0].name;
    c.spec = cases[0];
    c.ref = refs[0];
    c.off = runCase(cases[0], false, std::string(""), "OFF", true);
    std::printf("  OFF acceptance = %.6f %% (gate >= %.1f %%)\n", 100.0 * c.off.acceptance(),
                100.0 * kOffAcceptanceMin);
    c.controlPassed = c.off.acceptance() >= kOffAcceptanceMin;
    CHECK(c.controlPassed);
    std::printf("\n  [V10] reproduction of checkpoint 2g-h2e (printed T_prod, exact "
                "counts):\n");
    if (c.controlPassed) {
      double const nTrials = static_cast<double>(c.off.insideRadius());
      double sumT = 0.0;
      for (int i = 0; i < 3; ++i) {
        RunResult const r =
            runCase(cases[0], true, ref.path, "ON" + std::to_string(i + 1), false);
        double const T = r.transmission();
        sumT += T;
        std::printf("  [V10/ON%d] emitted %zu, accepted %zu, T_prod = %.6f, sigma = %.3e, "
                    "wall %.3f s | previous %.6f / %zu / %zu\n",
                    i + 1, r.emitted, r.insideRadius(), T, std::sqrt(T * (1.0 - T) / nTrials),
                    r.wall_s, kPrevV10[i], kPrevV10Emitted[i], kPrevV10Accepted[i]);
        CHECK(r.emitted == kPrevV10Emitted[i]);
        CHECK(r.insideRadius() == kPrevV10Accepted[i]);
        CHECK(std::fabs(T - kPrevV10[i]) <= kReproduceTol);
        if (i == 2) c.on = r;
      }
      c.tProd = sumT / 3.0;
      c.sigma = std::sqrt(c.tProd * (1.0 - c.tProd) / nTrials);
      c.onRun = true;
    }
    outcomes.push_back(c);
  }


  // ------------------------------------------------------- V5, Z30, Z60
  for (std::size_t i = 1; i < cases.size(); ++i) {
    CaseOutcome c;
    c.name = cases[i].name;
    c.spec = cases[i];
    c.ref = refs[i];
    c.off = runCase(cases[i], false, std::string(""), "OFF", true);
    c.controlPassed = c.off.acceptance() >= kOffAcceptanceMin;
    std::printf("  OFF acceptance = %.6f %% (gate >= %.1f %%), predicted cone coverage "
                "within %.1f m = %.4f\n",
                100.0 * c.off.acceptance(), 100.0 * kOffAcceptanceMin, kTelescopeRadius_m,
                c.ref.predictedAcceptance);
    if (!c.controlPassed) {
      std::printf("  STOP: %s OFF acceptance below %.1f %%; ON not run (pre-registered "
                  "rule); radius histogram above.\n",
                  c.name.c_str(), 100.0 * kOffAcceptanceMin);
      outcomes.push_back(c);
      continue;
    }
    c.on = runCase(cases[i], true, ref.path, "ON", false);
    c.tProd = c.on.transmission();
    double const nTrials = static_cast<double>(c.off.insideRadius());
    c.sigma = std::sqrt(c.tProd * (1.0 - c.tProd) / nTrials);
    c.onRun = true;
    outcomes.push_back(c);
  }


  // --------------------------------------------------------------- summary
  std::printf("\n[CherenkovE2E] SUMMARY (sigma with N = accepted_OFF of the same case)\n");
  std::printf("  %-4s %9s %9s %9s %10s %10s %10s %10s %10s %9s %9s\n", "case", "emitted",
              "accOFF", "accON", "T_prod", "sigma", "Tref_anch", "Tref_hc", "|dT|", "gate",
              "status");
  for (auto& c : outcomes) {
    if (c.onRun) {
      c.dTAnchor = std::fabs(c.tProd - c.ref.anchor);
      c.dTHc = std::fabs(c.tProd - c.ref.heightCorrected);
      c.rel = c.dTHc / c.ref.heightCorrected;
      c.gate = std::max(kGateAbsolute, 3.0 * c.sigma);
      c.stopPoint = kStopPointSigmaFactor * c.sigma + 1.0e-3;
    }
    std::printf("  %-4s %9zu %8.4f%% ", c.name.c_str(), c.off.emitted,
                100.0 * c.off.acceptance());
    if (c.onRun) {
      std::printf("%8.4f%% %10.6f %10.3e %10.6f %10.6f %10.3e %9.2e %9s\n",
                  100.0 * c.on.acceptance(), c.tProd, c.sigma, c.ref.anchor,
                  c.ref.heightCorrected, c.dTHc, c.gate,
                  (c.dTHc <= c.gate && c.rel <= kGateRelative) ? "PASS" : "FAIL");
    } else {
      std::printf("%8s %10s %10s %10.6f %10.6f %10s %9s %9s\n", "-", "-", "-",
                  c.ref.anchor, c.ref.heightCorrected, "-", "-", "STOPPED");
    }
  }

  std::printf("\n[CherenkovE2E] gates: absolute 1e-3, relative 5e-3, max(1e-3, 3 sigma), "
              "stop point 5 sigma + 1e-3\n");
  for (auto const& c : outcomes) {
    if (!c.onRun) {
      std::printf("  %-4s no T_prod: OFF acceptance %.4f %% < %.1f %% (pre-registered stop "
                  "rule)\n",
                  c.name.c_str(), 100.0 * c.off.acceptance(), 100.0 * kOffAcceptanceMin);
      continue;
    }
    std::printf("  %-4s |dT| = %.3e (vs anchor %.3e), rel = %.3e, gate = %.3e, stop = %.3e -> "
                "absolute:%s relative:%s gate:%s stop:%s\n",
                c.name.c_str(), c.dTHc, c.dTAnchor, c.rel, c.gate, c.stopPoint,
                (c.dTHc <= kGateAbsolute) ? "pass" : "FAIL",
                (c.rel <= kGateRelative) ? "pass" : "FAIL",
                (c.dTHc <= c.gate) ? "pass" : "FAIL",
                (c.dTHc <= c.stopPoint) ? "pass" : "FAIL");
    CHECK(c.dTHc <= kGateAbsolute);
    CHECK(c.rel <= kGateRelative);
    CHECK(c.dTHc <= c.gate);
    CHECK(c.dTHc <= c.stopPoint);
  }

  // ---------------------------------------------- corrected sanity lines
  std::printf("\n[CherenkovE2E] corrected sanity lines:\n");
  for (auto const& c : outcomes) {
    if (c.onRun) {
      std::printf("  %-4s tau_BoE ~ 0.19 | T_ref(anchor) = %.6f | T_ref(height-corr) = %.6f "
                  "| T_prod = %.6f\n",
                  c.name.c_str(), c.ref.anchor, c.ref.heightCorrected, c.tProd);
    } else {
      std::printf("  %-4s tau_BoE ~ 0.19 | T_ref(anchor) = %.6f | T_ref(height-corr) = %.6f "
                  "| T_prod = not measured (OFF acceptance %.4f %%)\n",
                  c.name.c_str(), c.ref.anchor, c.ref.heightCorrected,
                  100.0 * c.off.acceptance());
    }
  }

  // ------------------------------------------- geometry / stream contracts
  for (auto const& c : outcomes) {
    CHECK(c.off.absorptionDraws == 0);
    CHECK(c.off.emitted >= kBunchesPerRun);
    CHECK((c.off.outsideRadius == 0 || !c.controlPassed));
    CHECK(std::fabs(c.off.zProdMin_m - (c.spec.centreAltitude_m - kHalfWindow_m)) < 1.0);
    CHECK(std::fabs(c.off.zProdMax_m - (c.spec.centreAltitude_m + kHalfWindow_m)) < 1.0);
    CHECK(std::fabs(c.off.pathMean_m() - c.ref.predictedPath_m) < 5.0);
    if (c.onRun) {
      // Absorption is DETERMINISTIC since wave 4c: the transmission is applied as a
      // weight factor exp(-tau) and NO randomness is consumed on any route.  The old
      // contract here was "exactly one draw per emitted bunch"; it is now "no draws at
      // all", which is the stronger and checkable statement (and it is what keeps two
      // configurations with the same seed perfectly paired).
      CHECK(c.on.absorptionDraws == 0);
      CHECK(std::fabs(c.on.lambdaMin_nm - kLambdaMin_nm) < 0.01);
      CHECK(std::fabs(c.on.lambdaMax_nm - kLambdaMax_nm) < 0.01);
      CHECK(std::fabs(c.on.nProdMin - kRefractiveIndex) < 1.0e-9);
    }
  }
}

