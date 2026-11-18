/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * \file testCherenkovAbsorptionChord.cpp
 *
 * Wave 2b: the EXACT chord optical depth of AtmosphericAbsorption.
 *
 * The absorption table stores the cumulative column tau(h) = int_0^h alpha dh',
 * piecewise linear in altitude, so the extinction is constant inside each slab and the
 * integral of alpha along a straight chord has a closed form.  getOpticalDepthChord()
 * evaluates it; this test pins it three ways:
 *
 *   1. a vertical chord (impact parameter 0) must reproduce |tau(hE) - tau(hD)| EXACTLY
 *      -- that is the definition of the stored quantity, so any error in the slab
 *      bookkeeping shows up immediately;
 *   2. an independent midpoint quadrature along the same chord must converge TO the
 *      closed form as the step halves (and land within 1e-5 at the finer step);
 *   3. the Chapman air-mass approximation that getTransmissionSpherical used before is
 *      compared against it, so the bias the replacement removed stays measurable:
 *      negligible at the zenith, 3.5 % at 80 deg, 64 % at 88 deg, order unity at 89.9.
 */

#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace corsika::cherenkov;

namespace {

  constexpr double kR = 6.371e6; // the tests' earth radius
  constexpr double kLam = 400.0;

  std::string tablePath() {
    std::vector<std::string> const candidates = {
        "/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/"
        "abstable_MODTRAN_bottomup.dat",
        "modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat",
        "../modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat",
        "../../modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat"};
    for (auto const& p : candidates) {
      std::ifstream in(p);
      if (in.is_open()) return p;
    }
    return {};
  }

  /// The table's own tau row at 400 nm, read independently of the class under test.
  struct Grid {
    std::vector<double> h, tau;
    bool ok{false};
  };
  Grid loadGrid(std::string const& path) {
    Grid g;
    std::ifstream in(path);
    if (!in.is_open()) return g;
    std::size_t nl = 0, nh = 0;
    if (!(in >> nl >> nh) || nl == 0 || nh == 0) return g;
    std::vector<double> lam(nl);
    g.h.resize(nh);
    for (auto& v : lam) in >> v;
    for (auto& v : g.h) in >> v;
    std::vector<double> all(nl * nh);
    for (auto& v : all) in >> v;
    std::size_t row = 0;
    for (std::size_t i = 0; i < nl; ++i)
      if (std::fabs(lam[i] - kLam) < 1e-9) row = i;
    g.tau.assign(nh, 0.0);
    for (std::size_t j = 0; j < nh; ++j) g.tau[j] = all[row * nh + j];
    g.ok = true;
    return g;
  }

  double alphaAt(Grid const& g, double h) {
    std::size_t k = 0;
    while (k + 2 < g.h.size() && g.h[k + 1] <= h) ++k;
    if (h >= g.h.back()) return 0.0;
    return std::max(0.0, (g.tau[k + 1] - g.tau[k]) / (g.h[k + 1] - g.h[k]));
  }

} // namespace

TEST_CASE("AtmosphericAbsorption exact chord optical depth", "modules,cherenkov") {
  std::string const path = tablePath();
  REQUIRE(!path.empty());
  Grid const g = loadGrid(path);
  REQUIRE(g.ok);

  AtmosphericAbsorption abs;
  abs.load(path);
  REQUIRE(abs.isLoaded());

  SECTION("a vertical chord equals the column difference exactly") {
    double worst = 0.0;
    for (double hE : {0.0, 500.0, 1000.0, 3000.0, 10000.0, 30000.0}) {
      for (double hD : {0.0, 1000.0, 3000.0, 20000.0}) {
        if (hE == hD) continue;
        // a radial chord: its length IS the height difference, and the geometry is
        // unambiguous
        double const L = std::fabs(hD - hE);
        double const chord = abs.getOpticalDepthChord(kLam, hE, hD, L, kR);
        double const column =
            std::fabs(abs.getOpticalDepth(kLam, hE) - abs.getOpticalDepth(kLam, hD));
        double const rel = (column > 1e-12) ? std::fabs(chord - column) / column
                                           : std::fabs(chord - column);
        worst = std::max(worst, rel);
      }
    }
    std::printf("  vertical chords: worst relative difference from the column = %.2e\n", worst);
    CHECK(worst <= 1.0e-12);
  }

  SECTION("an independent midpoint quadrature converges to the closed form") {
    // The cases are given by the CHORD LENGTH, not by an angle: that is the contract
    // (r_emit, r_det, L) fixes the triangle.  The quadrature runs along the real
    // segment [0, L], so it is independent of the dip test -- if the closed form
    // picked the wrong intersection of the line with r = r_det it would disagree by
    // up to a factor 2, which is exactly how the near-horizontal descending branch
    // regressed before.
    std::printf("  %7s %7s %10s %9s %9s %12s %12s %10s %10s\n", "hE", "hD", "L [m]", "cos_zE",
                "cos_zD", "closed", "quadrature", "rel(ds)", "rel(ds/2)");
    double worstFine = 0.0;
    struct Case {
      double hE, hD, L;
      char const* note;
    };
    std::vector<Case> const cases = {
        {1000.0, 3000.0, 2001.0, "vertical-ish"},
        {1000.0, 3000.0, 5000.0, "slant"},
        {500.0, 3000.0, 20000.0, ""},
        {500.0, 3000.0, 200000.0, "grazing, dips to ~470 m"},
        {0.0, 3000.0, 100000.0, "from the ground, no dip"},
        {10000.0, 3000.0, 207283.0, "descending, NEAR side (no dip) <- regression"},
        {10000.0, 3000.0, 430817.0, "descending, FAR side (dips)"},
        {3100.0, 3000.0, 10005.0, "acceptance row 3100/10 km <- regression"},
    };
    for (auto const& c : cases) {
      double const rE = kR + c.hE, rD = kR + c.hD, L = c.L;
      // the same triangle the class must use; the law of cosines gives the angle from
      // the INWARD radial, so the zenith cosine is its negative
      double const czE = std::max(-1.0, std::min(1.0, -(rE * rE + L * L - rD * rD) / (2.0 * rE * L)));
      double const czD = std::max(-1.0, std::min(1.0, -(rD * rD + L * L - rE * rE) / (2.0 * rD * L)));
      double const st = std::sqrt(std::max(0.0, 1.0 - czE * czE));
      double const closed = abs.getOpticalDepthChord(kLam, c.hE, c.hD, L, kR);
      auto quadrature = [&](double ds) {
        double acc = 0.0;
        long const n = static_cast<long>(L / ds) + 1;
        for (long i = 0; i < n; ++i) {
          double const s = (double(i) + 0.5) * ds;
          double const z = rE + s * czE, x = s * st;
          acc += alphaAt(g, std::sqrt(z * z + x * x) - kR) * ds;
        }
        return acc;
      };
      double const ds = std::max(0.02, L / 2000000.0);
      double const q1 = quadrature(ds), q2 = quadrature(0.5 * ds);
      double const r1 = (closed > 1e-9) ? std::fabs(q1 - closed) / closed : std::fabs(q1 - closed);
      double const r2 = (closed > 1e-9) ? std::fabs(q2 - closed) / closed : std::fabs(q2 - closed);
      worstFine = std::max(worstFine, r2);
      std::printf("  %7.0f %7.0f %10.0f %9.5f %9.5f %12.6f %12.6f %10.2e %10.2e  %s\n", c.hE, c.hD,
                  L, czE, czD, closed, q2, r1, r2, c.note);
      CHECK(r2 <= 1.0e-5); // the finer step must match the closed form.  (No monotonicity
                           // assertion: the integrand is piecewise CONSTANT, so the
                           // midpoint error is set by where the samples fall relative to
                           // the slab boundaries and is not monotone in ds -- measured
                           // 2.9e-8 -> 1.9e-7 on the 500 m / cos 0.1 ray.)
    }
    std::printf("  worst at ds/2 = %.2e\n", worstFine);
  }

  SECTION("the removed Chapman approximation: size of the bias") {
    std::printf("  %8s %10s %12s %12s %10s\n", "zenith", "cos_z", "T_exact", "T_chapman", "rel");
    double relZenith = 0.0, relAt80 = 0.0, relAt89 = 0.0;
    for (double zen : {0.0, 30.0, 60.0, 80.0, 85.0, 89.0, 89.9}) {
      double const cz = std::cos(zen * M_PI / 180.0);
      double const st = std::sin(zen * M_PI / 180.0);
      double const rE = kR + 1000.0, rD = kR + 3000.0;
      double const p = rE * st;
      double const disc = rD * rD - p * p;
      if (!(disc >= 0.0)) continue;
      double const L = std::sqrt(disc) - rE * cz;
      double const Te = abs.getTransmissionSpherical(kLam, 1000.0, 3000.0, L, kR);
      double const Tc = abs.getTransmissionSphericalChapman(kLam, 1000.0, 3000.0, L, kR);
      double const rel = (Tc > 1e-12) ? std::fabs(Te - Tc) / Tc : std::fabs(Te - Tc);
      std::printf("  %8.1f %10.5f %12.6e %12.6e %10.2e\n", zen, cz, Te, Tc, rel);
      if (zen == 0.0) relZenith = rel;
      if (zen == 80.0) relAt80 = rel;
      if (zen == 89.0) relAt89 = rel;
      // The two models CROSS OVER: the Chapman air mass slightly OVER-absorbs around
      // 80 deg (T_exact 0.5112 > T_chapman 0.4938) and under-absorbs above ~85 deg.
      // What is asserted is the size of the disagreement, not its sign.
    }
    // regression guards on the measured bias (see diagnostics/b_table/README.md)
    CHECK(relZenith <= 1.0e-5); // 4.9e-8 measured
    CHECK(relAt80 >= 1.0e-2);   // 3.5e-2 measured (Chapman over-absorbs)
    CHECK(relAt89 >= 5.0e-1);   // 9.2e-1 measured (Chapman under-absorbs)
  }
}
