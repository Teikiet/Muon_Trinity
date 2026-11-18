/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <boost/filesystem.hpp>

#include <catch2/catch_all.hpp>

#include <corsika/modules/writers/InteractionWriter.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>

#include <corsika/setup/SetupTrajectory.hpp>

#include <SetupTestStack.hpp>
#include <SetupTestTrajectory.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

using DummyEnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

TEST_CASE("InteractionWriter", "process") {
  logging::set_level(logging::level::info);

  auto env = std::make_unique<media::Environment<media::IMediumModel>>();
  auto& universe = *(env->getUniverse());
  auto theMedium = media::Environment<media::IMediumModel>::createNode<Sphere>(
      Point{0_m, 0_m, 0_m}, 1_km * std::numeric_limits<double>::infinity());
  using MyHomogeneousModel = HomogeneousMedium<media::IMediumModel>;
  const auto density = 1_kg / (1_m * 1_m * 1_m);
  theMedium->setModelProperties<MyHomogeneousModel>(
      density, NuclearComposition({Code::Nitrogen}, {1.}));
  universe.addChild(std::move(theMedium));

  std::string const outputDir = "./output_interactions";

  // preparation
  if (boost::filesystem::exists(outputDir)) { boost::filesystem::remove_all(outputDir); }
  boost::filesystem::create_directory(outputDir);

  // setup empty particle stack
  test::Stack stack;
  stack.clear();

  // Make shower axis
  auto start = Point(-1_m, 0_m, 0_m);
  auto stop = Point(1_m, 0_m, 0_m);
  auto length = stop - start;
  media::ShowerAxis axis(start, length, *env);

  Plane const plane(Point(0_m, 0_m, 0_m), DirectionVector({0., 0., 1.}));

  ObservationPlane<setup::Tracking, ParticleWriterParquet> obsPlane{
      plane, DirectionVector({1., 0., 0.}),
      true,       // plane should "absorb" particles
      1e-6 * 1_m, // particles should be placed slightly past the plane
      false};

  InteractionWriter<setup::Tracking, ParticleWriterParquet> writer(axis, obsPlane);
  writer.startOfLibrary(outputDir);
  writer.startOfShower(0);

  CHECK(!writer.getInteractionCounter());

  // add primary particle to stack
  auto particle = stack.addParticle(std::make_tuple(
      Code::Proton, 10_GeV, DirectionVector({1, 0, 0}), Point(0_m, 0_m, 0_m), 0_ns));

  test::StackView view(particle);
  auto projectile = view.getProjectile();

  std::vector<Code> const particleList = {Code::PiPlus,   Code::PiMinus, Code::KPlus,
                                          Code::KMinus,   Code::K0Long,  Code::K0Short,
                                          Code::Electron, Code::MuPlus,  Code::NuE,
                                          Code::Neutron,  Code::NuMu};

  for (auto const proType : particleList)
    projectile.addSecondary(std::make_tuple(proType, 1_GeV, DirectionVector({1, 0, 0})));

  // pre-check on stack size
  CHECK(view.getEntries() == 11);
  CHECK(stack.getEntries() == 12);

  // correctly initialized counter
  CHECK(!writer.getInteractionCounter());

  // run the process and check for consistency
  writer.doSecondaries(view);
  CHECK(1 == writer.getInteractionCounter());
  CHECK(view.getEntries() == 11);
  CHECK(stack.getEntries() == 12);

  // ensure the counter bumps
  writer.doSecondaries(view);
  CHECK(2 == writer.getInteractionCounter());

  writer.endOfShower(0);
  writer.endOfLibrary();

  CHECK(boost::filesystem::exists(outputDir + "/interactions.parquet"));

  auto const config = writer.getConfig();
  CHECK(config["type"].as<std::string>() == "Interactions");

  auto const summary = writer.getSummary();
  CHECK(summary["shower_0"]["n_secondaries"].as<int>() == 11);

  // clean up
  if (boost::filesystem::exists(outputDir)) { boost::filesystem::remove_all(outputDir); }
}
