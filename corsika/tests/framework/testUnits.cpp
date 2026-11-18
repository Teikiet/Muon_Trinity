/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <array>
#include <sstream>

#include <catch2/catch_all.hpp>

using namespace corsika;
using Catch::Approx;

TEST_CASE("PhysicalUnits", "[Units]") {

  logging::set_level(logging::level::info);

  SECTION("Consistency") {
    CHECK(1_m / 1_m == Approx(1));
    // CHECK_FALSE( 1_m/1_s == 1 ); // static assert
  }

  SECTION("Constructors") {
    [[maybe_unused]] auto E1 = 10_GeV;
    CHECK(E1 == 10_GeV);

    LengthType l1 = 10_nm;
    [[maybe_unused]] auto l2 = l1;
    CHECK(l2 == l1);

    LengthType arr0[5];
    arr0[0] = 5_m;

    [[maybe_unused]] LengthType arr1[2] = {{1_mm}, {2_cm}};

    [[maybe_unused]] std::array<HEPEnergyType, 4> arr2; // empty array

    [[maybe_unused]] std::array<HEPEnergyType, 4> arr3 = {1_GeV, 1_eV, 5_MeV};

    auto p1 = 10_s * newton;
    CHECK(p1 == 10_s * newton);
  }

  SECTION("Powers in literal units") {
    CHECK(1_s / 1_ds == Approx(1e1));
    CHECK(1_m / 1_cm == Approx(1e2));
    CHECK(1_m / 1_mm == Approx(1e3));
    CHECK(1_V / 1_uV == Approx(1e6));
    CHECK(1_s / 1_ns == Approx(1e9));
    CHECK(1_eV / 1_peV == Approx(1e12));
    CHECK(1_A / 1_fA == Approx(1e15));
    CHECK(1_mol / 1_amol == Approx(1e18));
    CHECK(1_K / 1_zK == Approx(1e21));
    CHECK(1_K / 1_yK == Approx(1e24));
    CHECK(1_V / 1_rV == Approx(1e27));
    CHECK(1_V / 1_qV == Approx(1e30));
    CHECK(1_b / 1_mb == Approx(1e3));

    CHECK(1_A / 1_hA == Approx(1e-2));
    CHECK(1_m / 1_km == Approx(1e-3));
    CHECK(1_m / 1_Mm == Approx(1e-6));
    CHECK(1_V / 1_GV == Approx(1e-9));
    CHECK(1_s / 1_Ts == Approx(1e-12));
    CHECK(1_eV / 1_PeV == Approx(1e-15));
    CHECK(1_A / 1_EA == Approx(1e-18));
    CHECK(1_K / 1_ZK == Approx(1e-21));
    CHECK(1_mol / 1_Ymol == Approx(1e-24));
    CHECK(1_V / 1_RV == Approx(1e-27));
    CHECK(1_V / 1_QV == Approx(1e-30));

    CHECK(std::min(1_A, 2_A) == 1_A);
  }

  SECTION("Powers and units") {
    CHECK(1 * ampere / 1_A == Approx(1e0));
    CHECK(mega * bar / bar == Approx(1e6));
  }

  SECTION("Formulas") {
    const HEPEnergyType E2 = 20_GeV * 2;
    CHECK(E2 == 40_GeV);
    CHECK(E2 / 1_GeV == Approx(40));

    const MassType m = 1_kg;
    const SpeedType v = 1_m / 1_s;
    CHECK(m * v == 1_s * newton);

    const double lgE = log10(E2 / 1_GeV);
    CHECK(lgE == Approx(log10(40.)));

    const auto E3 = E2 + 100_GeV + pow(10, lgE) * 1_GeV;
    CHECK(E3 == 180_GeV);

    CHECK(sqrt(5_GeV * 5_GeV) / 5_GeV == Approx(1));
    CHECK(cbrt(static_pow<3>(5_GeV)) / 5_GeV == Approx(1));
  }

  SECTION("Output") {
    {
      const HEPEnergyType E = 5_eV;
      std::stringstream stream;
      stream << E;
      CHECK(stream.str() == std::string("5 eV"));
    }
    {
      const HEPEnergyType E = 5_EeV;
      std::stringstream stream;
      stream << E;
      CHECK(stream.str() == std::string("5e+18 eV"));
    }
  }

  SECTION("Special") {

    const LengthType farAway = std::numeric_limits<double>::infinity() * meter;
    CHECK(farAway > 100000_m);
    CHECK_FALSE(farAway < 1e19 * meter);
  }

  SECTION("static_pow") {
    double x = 235.7913;
    CHECK(1 == static_pow<0, double>(x));
    CHECK(x == static_pow<1, double>(x));
    CHECK(x * x == static_pow<2, double>(x));
    CHECK(1 / x == static_pow<-1, double>(x));
    CHECK(1 / x / x == static_pow<-2, double>(x));
  }

  SECTION("HEP/SI conversion") {
    auto const invEnergy = 1 / 197.326978_MeV; // should be convertible to length or time

    LengthType const length = convert_HEP_to_SI<LengthType::dimension_type>(invEnergy);
    CHECK((length / 1_fm) == Approx(1));

    TimeType const time = convert_HEP_to_SI<TimeType::dimension_type>(invEnergy);
    CHECK((time / (1_fm / constants::c)) == Approx(1));

    auto const protonMass = 938.272'081'3_MeV; // convertible to mass or SI energy
    MassType protonMassSI = convert_HEP_to_SI<MassType::dimension_type>(protonMass);
    CHECK((protonMassSI / 1.672'621'898e-27_kg) == Approx(1));
    CHECK((protonMassSI / (1.007'276 * constants::u)) == Approx(1));
  }

  SECTION("SI/HEP conversion") {
    CHECK(convert_SI_to_HEP(constants::c) == Approx(1));
    CHECK(convert_SI_to_HEP(constants::hBar) == Approx(1));

    {
      auto const invLength = 1 / 197.326978_fm; // should be convertible to HEPEnergy
      HEPEnergyType const energy = convert_SI_to_HEP(invLength);
      CHECK(energy / 1_MeV == Approx(1));
    }

    CHECK(convert_SI_to_HEP(6.5823e-25_s) * 1_GeV == Approx(1).epsilon(1e-4));

    CHECK(convert_SI_to_HEP(3.8938e-32 * meter * meter) * 1_GeV * 1_GeV ==
          Approx(1).epsilon(1e-4));
  }
}
