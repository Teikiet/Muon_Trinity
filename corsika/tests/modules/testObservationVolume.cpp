/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>
#include <corsika/modules/ObservationVolume.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/Box.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/modules/writers/WriterOff.hpp>

#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>
#include <SetupTestTrajectory.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

using namespace corsika;
using Catch::Approx;

// some geometries such as box do not support curved track yet.
using Tracking = tracking_line::Tracking;

TEST_CASE("ObservationVolume", "interface") {

  auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Oxygen);
  auto const& cs = *csPtr; // Déréférencer le pointeur pour obtenir une référence
  [[maybe_unused]] auto const& env_dummy = env;
  [[maybe_unused]] auto const& node_dummy = nodePtr;

  test::Stack stack;

  Point center(cs, {0_m, 0_m, 0_m});

  Point pointInit(cs, {-10_m, 0_m, 0_m});
  DirectionVector direction(cs, {1., 0., 0.});
  VelocityVector velocity = constants::c * direction;
  auto particle =
      stack.addParticle(std::make_tuple(Code::NuE, 1_GeV, direction, pointInit, 0_ns));
  Line line(pointInit, velocity);
  StraightTrajectory trajectory(line, 1_s);

  SECTION("sphere volume") {
    Sphere sphere(center, 5_m);
    ObservationVolume<tracking_line::Tracking, Sphere, WriterOff> obs(sphere);

    Step step(particle, trajectory);
    LengthType const length = obs.getMaxStepLength(particle, trajectory);
    ProcessReturn const ret = obs.doContinuous(step, true);

    CHECK(length / 5_m == Approx(1).margin(1e-4));
    CHECK(ret == ProcessReturn::ParticleAbsorbed);

    // particle does not reach volume:
    {
      setup::Trajectory no_hit_track =
          setup::testing::make_track<setup::Trajectory>(line, 1_nm / constants::c);
      LengthType const no_hit = obs.getMaxStepLength(particle, no_hit_track);
      CHECK(no_hit == std::numeric_limits<double>::infinity() * 1_m);
    }

    // particle inside volume:
    {
      particle.setPosition({cs, {-1_m, 0_m, 0_m}});
      setup::Trajectory inside_track =
          setup::testing::make_track<setup::Trajectory>(line, 1_nm / constants::c);
      LengthType const hit_length = obs.getMaxStepLength(particle, inside_track);
      ProcessReturn const ret = obs.doContinuous(step, true);
      CHECK(hit_length / 1_m == Approx(0).margin(1e-12));
      CHECK(ret == ProcessReturn::ParticleAbsorbed);
    }
  }

  SECTION("box volume") {
    Box box(cs, 10_m);
    ObservationVolume<tracking_line::Tracking, Box, WriterOff> obs(box);

    Step step(particle, trajectory);
    LengthType const length = obs.getMaxStepLength(particle, trajectory);
    ProcessReturn const ret = obs.doContinuous(step, true);

    CHECK(length / 5_m == Approx(1).margin(1e-4));
    CHECK(ret == ProcessReturn::ParticleAbsorbed);
  }

  SECTION("output") {
    Sphere sphere(center, 5_m);
    // default writer will be used here
    ObservationVolume<tracking_line::Tracking, Sphere> obs(sphere);
    auto const cfg = obs.getConfig();
    CHECK(cfg["type"]);
  }
}