/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#define protected public // to also test the internal state of objects

#include <corsika/stack/VectorStack.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace std;
using Catch::Approx;

TEST_CASE("VectorStack", "stack") {

  logging::set_level(logging::level::info);

  SECTION("read+write") {

    VectorStack s;
    s.addParticle(std::make_tuple(Code::Electron, 1.5_GeV, DirectionVector({1, 0, 0}),
                                  Point(1 * meter, 1 * meter, 1 * meter), 100_s));

    // read
    CHECK(s.getEntries() == 1);
    CHECK(s.getSize() == 1);
    auto pout = s.getNextParticle();
    CHECK(pout.getPID() == Code::Electron);
    CHECK(pout.getEnergy() == 1.5_GeV + Electron::mass);
    CHECK(pout.getKineticEnergy() == 1.5_GeV);
    CHECK(pout.getTime() == 100_s);
    CHECK(pout.getChargeNumber() == -1);

    s.clear();
    CHECK(s.getEntries() == 0);
    CHECK(s.getSize() == 0);
  }

  SECTION("write+delete") {

    VectorStack s;
    for (int i = 0; i < 99; ++i)

      s.addParticle(std::make_tuple(Code::Electron, 1_GeV, DirectionVector({1, 0, 0}),
                                    Point(1 * meter, 1 * meter, 1 * meter), 100_s));

    CHECK(s.getSize() == 99);

    for (int i = 0; i < 99; ++i) s.getNextParticle().erase();

    CHECK(s.getEntries() == 0);
    CHECK(s.getSize() == 1);
  }

  SECTION("stack operations") {

    VectorStack s;
    // add 99 particles, each 10th particle is a positron
    // i=9, 19, 29, etc. are positron
    for (int i = 0; i < 99; ++i) {
      Code const pid = ((i + 1) % 10 == 0) ? Code::Positron : Code::Electron;
      s.addParticle(std::make_tuple(pid, i * 1.5_GeV - Electron::mass,
                                    DirectionVector({1, 0, 0}),
                                    Point(1 * meter, 1 * meter, 1 * meter), 100_s));
    }

    // copy
    {
      s.copy(s.begin() + 9, s.begin() + 10); // positron to electron
      const auto& p9 = s.cbegin() + 9;
      const auto& p10 = s.cbegin() + 10;

      CHECK(p9.getPID() == Code::Positron);
      CHECK(p9.getEnergy() == 9 * 1.5_GeV);
      CHECK(p9.getTime() == 100_s);

      CHECK(p10.getPID() == Code::Positron);
      CHECK(p10.getEnergy() == 9 * 1.5_GeV);
      CHECK(p10.getTime() == 100_s);
    }

    // copy
    {
      s.copy(s.begin() + 93, s.begin() + 9); // electron to positron
      const auto& p93 = s.cbegin() + 93;
      const auto& p9 = s.cbegin() + 9;

      CHECK(p9.getPID() == Code::Electron);
      CHECK(p9.getEnergy() == 93 * 1.5_GeV);
      CHECK(p9.getTime() == 100_s);

      CHECK(p93.getPID() == Code::Electron);
      CHECK(p93.getEnergy() == 93 * 1.5_GeV);
      CHECK(p93.getTime() == 100_s);
    }

    // copy
    {
      s.copy(s.begin() + 89, s.begin() + 79); // positron to positron
      const auto& p89 = s.cbegin() + 89;
      const auto& p79 = s.cbegin() + 79;

      CHECK(p89.getPID() == Code::Positron);
      CHECK(p89.getEnergy() == 89 * 1.5_GeV);
      CHECK(p89.getTime() == 100_s);

      CHECK(p79.getPID() == Code::Positron);
      CHECK(p79.getEnergy() == 89 * 1.5_GeV);
      CHECK(p79.getTime() == 100_s);
    }

    // invalid copy
    { CHECK_THROWS(s.copy(s.begin(), s.end() + 1000)); }

    // swap
    {
      s.swap(s.begin() + 11, s.begin() + 10);
      const auto& p11 = s.cbegin() + 11; // now: positron
      const auto& p10 = s.cbegin() + 10; // now: electron

      CHECK(p11.getPID() == Code::Positron);
      CHECK(p11.getEnergy() == 9 * 1.5_GeV);
      CHECK(p11.getTime() == 100_s);

      CHECK(p10.getPID() == Code::Electron);
      CHECK(p10.getEnergy() == 11 * 1.5_GeV);
      CHECK(p10.getTime() == 100_s);
    }

    // swap two positrons
    {
      s.swap(s.begin() + 29, s.begin() + 59);
      const auto& p29 = s.cbegin() + 29;
      const auto& p59 = s.cbegin() + 59;

      CHECK(p29.getPID() == Code::Positron);
      CHECK(p29.getEnergy() == 59 * 1.5_GeV);
      CHECK(p29.getTime() == 100_s);

      CHECK(p59.getPID() == Code::Positron);
      CHECK(p59.getEnergy() == 29 * 1.5_GeV);
      CHECK(p59.getTime() == 100_s);
    }

    for (int i = 0; i < 99; ++i) s.last().erase();
    CHECK(s.getEntries() == 0);
  }

  SECTION("secondary generation") {

    // check inheritance of time and location

    VectorStack s;
    auto p1 = s.addParticle(std::make_tuple(
        Code::Photon, 1.5_GeV - Electron::mass, DirectionVector({1, 0, 0}),
        Point(1 * meter, 1 * meter, 1 * meter), 10_s));
    auto p2 = p1.addSecondary(
        std::make_tuple(Code::Photon, 10_GeV, DirectionVector({1, 0, 0})));
    auto const delta_pos = Vector<length_d>(1_m, 0_m, 0_m);
    auto const delta_time = 1_s;
    auto p3 = p1.addSecondary(std::make_tuple(
        Code::Photon, 10_GeV, DirectionVector({1, 0, 0}), delta_pos, delta_time));

    CHECK(p1.getTime() / 10_s == Approx(1));
    CHECK(p1.getPosition().getX() / 1_m == Approx(1));
    CHECK(p1.getPosition().getY() / 1_m == Approx(1));
    CHECK(p1.getPosition().getZ() / 1_m == Approx(1));

    CHECK(p2.getTime() / 10_s == Approx(1));
    CHECK(p2.getPosition().getX() / 1_m == Approx(1));
    CHECK(p2.getPosition().getY() / 1_m == Approx(1));
    CHECK(p2.getPosition().getZ() / 1_m == Approx(1));

    CHECK(p3.getTime() / 11_s == Approx(1));
    CHECK(p3.getPosition().getX() / 1_m == Approx(2));
    CHECK(p3.getPosition().getY() / 1_m == Approx(1));
    CHECK(p3.getPosition().getZ() / 1_m == Approx(1));
  }
}
