/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <catch2/catch_all.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace corsika;

namespace {
  // Synthetic exponential atmosphere on a uniform 10 m grid (mirrors the fine
  // grid the application builds): N(h) = N0*exp(-h/H).
  struct SyntheticAtmosphere {
    static constexpr double R = 6371000.0; // sea level absolute z [m]
    static constexpr double H = 8000.0;    // scale height [m]
    static constexpr double N0 = 3.27e-4;  // sea-level refractivity
    std::vector<double> alt_km;
    std::vector<double> lnN;
    SyntheticAtmosphere() {
      for (double h_km = 0.0; h_km <= 120.0 + 1e-9; h_km += 0.01) {
        alt_km.push_back(h_km);
        lnN.push_back(std::log(N0 * std::exp(-(h_km * 1000.0) / H)));
      }
    }
    double N(double h) const { return N0 * std::exp(-h / H); }
  };
} // namespace

TEST_CASE("RefractiveIndexTable", "modules,cherenkov") {
  SyntheticAtmosphere atm;
  cherenkov::RefractiveIndexTable t(atm.alt_km, atm.lnN, SyntheticAtmosphere::R);
  REQUIRE(t.valid());

  SECTION("N(z) reproduces the analytic profile") {
    double maxErr = 0.0;
    for (double h = 1000.0; h < 110000.0; h += 1234.0) {
      double const exact = atm.N(h);
      maxErr = std::max(maxErr, std::abs(t.N(SyntheticAtmosphere::R + h) - exact) / exact);
    }
    INFO("max rel N error = " << maxErr);
    REQUIRE(maxErr < 1e-9);
  }

  SECTION("A'(z) == N(z) (central finite differences)") {
    double maxErr = 0.0;
    double const dz = 0.5;
    for (double h = 1000.0; h < 110000.0; h += 4999.0) {
      double const z = SyntheticAtmosphere::R + h + 0.137; // avoid cell edges
      double const dA = (t.A(z + dz) - t.A(z - dz)) / (2.0 * dz);
      maxErr = std::max(maxErr, std::abs(dA - t.N(z)) / t.N(z));
    }
    INFO("max rel |A'(z)-N(z)| = " << maxErr);
    REQUIRE(maxErr < 1e-6);
  }

  SECTION("closed form matches a 10000-segment integral at full precision") {
    double maxRel = 0.0;
    for (double const dz_m : {5.0, 10.0, 50.0, 1000.0, 37000.0}) {
      double const z0 = SyntheticAtmosphere::R + 37000.0;
      double const z1 = z0 - dz_m;
      double const uz = -1.0; // vertical downward, path length = dz_m
      double const closed = (t.A(z1) - t.A(z0)) / uz;
      int const n = 10000;
      double num = 0.0;
      double const dstep = dz_m / n;
      for (int i = 0; i < n; ++i) num += t.N(z0 + (i + 0.5) * dstep * uz) * dstep;
      double const rel = std::abs(closed - num) / num;
      std::printf("  dz=%8.3f m  closed=%.17g  num=%.17g  rel=%.3e\n", dz_m, closed, num,
                  rel);
      maxRel = std::max(maxRel, rel);
    }
    INFO("max rel closed-vs-10000seg = " << maxRel);
    REQUIRE(maxRel < 1e-6);
  }

  SECTION("out-of-range z: clamp below, saturate above") {
    REQUIRE(t.N(t.zMin() - 100.0) == t.N(t.zMin()));
    REQUIRE(t.A(t.zMax() + 50000.0) >= t.A(t.zMax()));
    REQUIRE(t.N(t.zMax() + 50000.0) == 0.0); // 1e-10 floor -> 0
  }
}
