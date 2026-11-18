/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/media/density/DensityFunction.hpp>
#include <corsika/media/density_and_composition/FlatExponential.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/VolumeTreeNode.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <catch2/catch_all.hpp>

// using namespace
using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

const auto density = 1_kg / (1_m * 1_m * 1_m);

auto setupEnvironmentHomogenous(Code vTargetCode) {
  // setup environment, geometry
  auto env = std::make_unique<media::Environment<media::IMediumModel>>();
  auto& universe = *(env->getUniverse());
  auto const& cs = env->getCoordinateSystem();

  auto theMedium = media::Environment<media::IMediumModel>::createNode<Sphere>(
      Point(cs, 0_m, 0_m, 0_m), 1_km * std::numeric_limits<double>::infinity());

  using MyHomogeneousModel = HomogeneousMedium<media::IMediumModel>;
  theMedium->setModelProperties<MyHomogeneousModel>(
      density, NuclearComposition({vTargetCode}, {1.}));

  auto const* nodePtr = theMedium.get();
  universe.addChild(std::move(theMedium));

  return std::make_tuple(std::move(env), std::ref(cs), nodePtr);
}

TEST_CASE("Homogeneous Density") {

  logging::set_level(logging::level::info);

  auto [env, cs, nodePtr] = setupEnvironmentHomogenous(Code::Nitrogen);
  [[maybe_unused]] auto const& env_dummy = env;
  [[maybe_unused]] auto const& node_dummy = nodePtr;

  auto const observationHeight = 0_km;
  auto const injectionHeight = 10_km;
  auto const t = -observationHeight + injectionHeight;
  Point const showerCore(cs, 0_m, 0_m, observationHeight);
  Point const injectionPos = showerCore + Vector<dimensionless_d>{cs, {0, 0, 1}} * t;

  media::ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos), *env,
                                     true, // -> throw exceptions
                                     20};  // -> number of bins

  CHECK(showerAxis.getSteplength() == 500_m);

  CHECK(showerAxis.getMaximumX() / (10_km * density) == Approx(1).epsilon(1e-8));

  CHECK(showerAxis.getMinimumX() == 0_g / square(1_cm));

  const Point p{cs, 10_km, 20_km, 8.3_km};
  CHECK(showerAxis.getProjectedX(p) / (1.7_km * density) == Approx(1).epsilon(1e-8));

  const LengthType d = 6.789_km;
  CHECK(showerAxis.getX(d) / (d * density) == Approx(1).epsilon(1e-8));

  const Vector<dimensionless_d> dir{cs, {0, 0, -1}};
  CHECK(showerAxis.getDirection().getComponents(cs) == dir.getComponents(cs));
  CHECK(showerAxis.getStart().getCoordinates() == injectionPos.getCoordinates());

  CHECK_THROWS(showerAxis.getX(-1_m));
  CHECK_THROWS(showerAxis.getX((injectionPos - showerCore).getNorm() + 1_m));

  media::ShowerAxis const showerAxisNoThrow{injectionPos, (showerCore - injectionPos),
                                            *env,
                                            false, // -> do not throw exceptions
                                            20};   // -> number of bins
  CHECK(showerAxisNoThrow.getX(-1_m) == showerAxis.getMinimumX());
  CHECK(showerAxisNoThrow.getX((injectionPos - showerCore).getNorm() + 1_m) ==
        showerAxis.getMaximumX());
}

auto setupEnvironmentExponential() {
  // setup environment, geometry
  auto env = std::make_unique<media::Environment<media::IMediumModel>>();
  auto& universe = *(env->getUniverse());
  auto const& cs = env->getCoordinateSystem();

  auto theMedium = media::Environment<media::IMediumModel>::createNode<Sphere>(
      Point{cs, 0_m, 0_m, 0_m}, 1_km * std::numeric_limits<double>::infinity());

  auto const point = Point{cs, 0_m, 0_m, 0_m};
  auto const axis = DirectionVector{cs, {0, 0, 1}};
  auto const lambda = 2_m;

  using MyExponentialModel = FlatExponential<media::IMediumModel>;
  theMedium->setModelProperties<MyExponentialModel>(
      point, axis, density, lambda, NuclearComposition({Code::Nitrogen}, {1.}));

  auto const* nodePtr = theMedium.get();
  universe.addChild(std::move(theMedium));

  return std::make_tuple(std::move(env), std::ref(cs), nodePtr, point, axis, lambda);
}

TEST_CASE("Exponential Density") {
  logging::set_level(logging::level::trace);

  // rho_0 Exp[ (r - x0) \cdot n / lambda ].
  auto [env, cs, nodePtr, start, axis, lambda] = setupEnvironmentExponential();

  auto const totalLength = 5 * abs(lambda);
  auto const end = start + axis * totalLength;
  media::ShowerAxis const showerAxis{start, end, *env,
                                     true, // -> throw exceptions
                                     50};  // -> number of bins

  CHECK(showerAxis.getSteplength() == totalLength / 50.0);
  CHECK(showerAxis.getMinimumX() == 0_g / square(1_m));
  CHECK(showerAxis.getMaximumX() ==
        density * lambda * (std::exp(totalLength / lambda) - 1.0));
  auto const units = 1_kg / square(1_m);
  CHECK(showerAxis.getMaximumX() / units ==
        Approx(showerAxis.getX(totalLength - 1_nm) / units).epsilon(1e-6));
  CHECK(showerAxis.getX(abs(lambda)) / units ==
        Approx(density * abs(lambda) * (std::exp(1.0) - 1.0) / units).epsilon(1e-8));
}

auto setupEnvironmentHardBoundary(Code vTargetCode) {

  // setup environment, geometry
  auto env = std::make_unique<media::Environment<media::IMediumModel>>();
  auto& universe = *(env->getUniverse());
  auto const& cs = env->getCoordinateSystem();

  auto theMedium = media::Environment<media::IMediumModel>::createNode<Sphere>(
      Point(cs, 0_m, 0_m, 0_m), 1_km * std::numeric_limits<double>::infinity());

  using MyHomogeneousModel = HomogeneousMedium<media::IMediumModel>;
  theMedium->setModelProperties<MyHomogeneousModel>(
      density, NuclearComposition({vTargetCode}, {1.}));

  auto notDenseMedium = media::Environment<media::IMediumModel>::createNode<Sphere>(
      Point(cs, 0_m, 0_m, 0_m), 1_m);
  using MyHomogeneousModel = HomogeneousMedium<media::IMediumModel>;
  notDenseMedium->setModelProperties<MyHomogeneousModel>(
      density * 1e-9, NuclearComposition({vTargetCode}, {1.}));

  theMedium->addChild(std::move(notDenseMedium));
  universe.addChild(std::move(theMedium));

  return std::make_tuple(std::move(env), std::ref(cs));
}

TEST_CASE("Hard Boundary") {
  auto [env, cs] = setupEnvironmentHardBoundary(Code::Nitrogen);

  Point const start(cs, 0_m, 0_m, 0_m);
  Point const end(cs, 0_m, 0_m, 2_m);
  media::ShowerAxis const showerAxis{start, end, *env,
                                     true, // -> throw exceptions
                                     50};  // -> number of bins

  auto const units = 1_kg / square(1_m);
  // Integrate out to the hard boundary to confirm not-dense region's existance
  CHECK(showerAxis.getX(1_m - 1_nm) / units ==
        Approx(density * 1e-9 * 1_m / units).epsilon(1e-8));

  // Integrating just beyond the boundary essentially equal to just the dense portion
  auto delta = 0.1 * 1_m;
  CHECK(showerAxis.getX(1_m + delta) / units ==
        Approx(density * delta / units).epsilon(1e-6));
  delta = 0.99 * 1_m;
  CHECK(showerAxis.getX(1_m + delta) / units ==
        Approx(density * delta / units).epsilon(1e-6));
}
