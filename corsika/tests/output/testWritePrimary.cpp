/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <boost/filesystem.hpp>
#include <catch2/catch_all.hpp>

#include <corsika/framework/core/Logging.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

using namespace corsika;
using Catch::Approx;

auto setupWritePrimaryEnv() {
  // Set up all of the background
  auto env = std::make_unique<media::Environment<media::IMediumModel>>();
  auto const& cs = env->getCoordinateSystem();

  Plane const plane(Point(cs, {0_m, 0_m, 0_m}), DirectionVector(cs, {0., 0., 1.}));

  ObservationPlane<setup::Tracking, ParticleWriterParquet> obsPlane{
      plane, DirectionVector(cs, {1., 0., 0.}),
      true,       // plane should "absorb" particles
      1e-6 * 1_m, // ignored for absorbing planes
      false};

  return std::make_tuple(std::move(env), std::ref(cs), obsPlane);
}

std::string resetPrimaryWriterOutputDir() {
  // Helper function to test the output directory
  std::string const file_dir = "./output_dir_write_prim";
  if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
  boost::filesystem::create_directory(file_dir);

  return file_dir;
}

TEST_CASE("PrimaryWriter") {
  logging::set_level(logging::level::info);

  SECTION("basic_functionality") {

    auto [env, cs, obsPlane] = setupWritePrimaryEnv();

    PrimaryWriter<setup::Tracking, ParticleWriterParquet> primaryWriter(obsPlane);

    YAML::Node summary = primaryWriter.getSummary();
    YAML::Node config = primaryWriter.getConfig();
  }

  SECTION("empty_output") {

    auto [env, cs, obsPlane] = setupWritePrimaryEnv();

    PrimaryWriter<setup::Tracking, ParticleWriterParquet> primaryWriter(obsPlane);

    auto const file_dir = resetPrimaryWriterOutputDir(); // Make output dir

    primaryWriter.startOfLibrary(file_dir);
    primaryWriter.startOfShower(0);
    primaryWriter.endOfShower(0);

    YAML::Node summary = primaryWriter.getSummary();
    CHECK(!summary.size());
  }

  SECTION("single_shower") {
    auto [env, cs, obsPlane] = setupWritePrimaryEnv();

    // Make the primary information
    auto const primaryProps =
        std::make_tuple(Code::Proton, 1_GeV, DirectionVector(cs, {0., 0., -1.}),
                        Point(cs, {0_m, 0_m, 1_m}), 0_ns);

    PrimaryWriter<setup::Tracking, ParticleWriterParquet> primaryWriter(obsPlane);

    auto const file_dir = resetPrimaryWriterOutputDir(); // Make output dir

    primaryWriter.startOfLibrary(file_dir);
    primaryWriter.recordPrimary(primaryProps);
    primaryWriter.startOfShower(0);
    primaryWriter.endOfShower(0);
    primaryWriter.endOfLibrary();

    YAML::Node summary = primaryWriter.getSummary();
    CHECK(1 == summary.size()); // A particle was recorded

    CHECK(2212 == summary["shower_0"]["pdg"].as<int>());
    CHECK(summary["shower_0"]["kinetic_energy"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["x"].as<double>() == Approx(0));
    CHECK(summary["shower_0"]["y"].as<double>() == Approx(0));
    CHECK(summary["shower_0"]["z"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["nx"].as<double>() == Approx(0));
    CHECK(summary["shower_0"]["ny"].as<double>() == Approx(0));
    CHECK(summary["shower_0"]["nz"].as<double>() == Approx(-1));
    CHECK(summary["shower_0"]["time"].as<double>() == Approx(0));
  }

  SECTION("multiple_showers") {
    auto [env, cs, obsPlane] = setupWritePrimaryEnv();

    // Make the primary information
    auto const primaryProps =
        std::make_tuple(Code::Proton, 1_GeV, DirectionVector(cs, {0., 0., -1.}),
                        Point(cs, {0_m, 0_m, 1_m}), 0_ns);

    PrimaryWriter<setup::Tracking, ParticleWriterParquet> primaryWriter(obsPlane);

    auto const file_dir = resetPrimaryWriterOutputDir(); // Make output dir

    unsigned int const n_shower = 5;
    primaryWriter.startOfLibrary(file_dir);
    for (unsigned int ishower = 0; ishower < n_shower; ishower++) {
      primaryWriter.recordPrimary(primaryProps);
      primaryWriter.startOfShower(0);
      primaryWriter.endOfShower(0);
    }
    primaryWriter.endOfLibrary();

    YAML::Node summary = primaryWriter.getSummary();
    CHECK(n_shower == summary.size()); // There was an entry for each shower
  }
}
