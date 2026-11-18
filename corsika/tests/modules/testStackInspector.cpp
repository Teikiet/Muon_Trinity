/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/modules/StackInspector.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <../framework/testCascade.hpp> //! \todo fix this

using namespace corsika;

TEST_CASE("StackInspector", "modules") {

  logging::set_level(logging::level::info);

  Point const origin(0_m, 0_m, 0_m);
  VelocityVector v(0_m / second, 0_m / second, 1_m / second);
  Line line(origin, v);
  StraightTrajectory track(line, 10_s);

  TestCascadeStack stack;
  stack.clear();
  HEPEnergyType E0 = 100_GeV;
  stack.addParticle(std::make_tuple(Code::Electron, 1_GeV, DirectionVector({0, 0, -1}),
                                    Point(0_m, 0_m, 10_km), 0_ns));
  stack.addParticle(std::make_tuple(get_nucleus_code(16, 8), 1_GeV,
                                    DirectionVector({0, 0, -1}), Point(0_m, 0_m, 10_km),
                                    0_ns));

  SECTION("interface") {

    StackInspector<TestCascadeStack> model(1, true, E0);
    model.doStack(stack);
    // there are no actions, nothing to check...
  }
}
