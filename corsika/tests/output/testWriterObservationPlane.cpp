/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <boost/filesystem.hpp>

#include <corsika/modules/writers/ParticleWriterParquet.hpp>

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>

using namespace corsika;
using Catch::Approx;

struct TestWriterPlane : public ParticleWriterParquet {

  TestWriterPlane(bool const printZ)
      : ParticleWriterParquet(printZ) {}

  YAML::Node getConfig() const { return YAML::Node(); }

  void checkWrite() {
    ParticleWriterParquet::write(Code::Unknown, 1_GeV, 2_m, 3_m, 0_m, 1, 0, 0, 1_us, 1.0);
    ParticleWriterParquet::write(Code::Proton, 1_GeV, 2_m, 3_m, 0_m, -1, 0, 0, 1_us, 1.0);
    ParticleWriterParquet::write(Code::MuPlus, 1_GeV, 2_m, 3_m, 0_m, 0, 1, 0, 1_us, 1.0);
    ParticleWriterParquet::write(Code::MuMinus, 1_GeV, 2_m, 3_m, 0_m, 0, -1, 0, 1_us,
                                 1.0);
    ParticleWriterParquet::write(Code::Photon, 1_GeV, 2_m, 3_m, 0_m, 0, 0, -1, 1_us, 1.0);
  }
};

TEST_CASE("ObservationPlaneWriterParquet") {

  logging::set_level(logging::level::info);

  SECTION("standard") {

    std::string const file_dir = "./output_dir_obs_plane";

    // preparation
    if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
    boost::filesystem::create_directory(file_dir);

    TestWriterPlane test(true);
    test.startOfLibrary(file_dir);

    test.startOfShower(0);
    test.checkWrite();
    test.endOfShower(0);

    test.endOfLibrary();

    CHECK(boost::filesystem::exists(file_dir + "/particles.parquet"));

    auto const summary = test.getSummary();

    CHECK(summary["shower_0"]["hadron"]["count"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["muon"]["count"].as<double>() == Approx(2));
    CHECK(summary["shower_0"]["em"]["count"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["other"]["count"].as<double>() == Approx(1));

    CHECK(summary["shower_0"]["hadron"]["kinetic_energy"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["muon"]["kinetic_energy"].as<double>() == Approx(2));
    CHECK(summary["shower_0"]["em"]["kinetic_energy"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["other"]["kinetic_energy"].as<double>() == Approx(1));

    // clean things up
    if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
  }

  SECTION("without_z") {

    std::string const file_dir = "./output_dir_obs_plane_no_z";

    // preparation
    if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
    boost::filesystem::create_directory(file_dir);

    TestWriterPlane test(false); // do not print the z-coord
    test.startOfLibrary(file_dir);

    test.startOfShower(0);
    test.checkWrite();
    test.endOfShower(0);

    test.endOfLibrary();

    CHECK(boost::filesystem::exists(file_dir + "/particles.parquet"));

    auto const summary = test.getSummary();

    CHECK(summary["shower_0"]["hadron"]["count"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["muon"]["count"].as<double>() == Approx(2));
    CHECK(summary["shower_0"]["em"]["count"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["other"]["count"].as<double>() == Approx(1));

    CHECK(summary["shower_0"]["hadron"]["kinetic_energy"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["muon"]["kinetic_energy"].as<double>() == Approx(2));
    CHECK(summary["shower_0"]["em"]["kinetic_energy"].as<double>() == Approx(1));
    CHECK(summary["shower_0"]["other"]["kinetic_energy"].as<double>() == Approx(1));

    // clean things up
    if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
  }

  SECTION("reset_between_showers") {
    // Check that the counters are reset between showers

    std::string const file_dir = "./output_dir_obs_plane_reset";

    // preparation
    if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
    boost::filesystem::create_directory(file_dir);

    TestWriterPlane test(false); // do not print the z-coord
    test.startOfLibrary(file_dir);

    test.startOfShower(0);
    test.checkWrite();
    test.endOfShower(0);

    test.startOfShower(1);
    test.checkWrite();
    test.endOfShower(1);

    auto const summary = test.getSummary();

    CHECK(summary["shower_0"]["hadron"]["count"].as<double>() ==
          summary["shower_1"]["hadron"]["count"].as<double>());
    CHECK(summary["shower_0"]["muon"]["count"].as<double>() ==
          summary["shower_1"]["muon"]["count"].as<double>());
    CHECK(summary["shower_0"]["em"]["count"].as<double>() ==
          summary["shower_1"]["em"]["count"].as<double>());
    CHECK(summary["shower_0"]["other"]["count"].as<double>() ==
          summary["shower_1"]["other"]["count"].as<double>());

    CHECK(boost::filesystem::exists(file_dir + "/particles.parquet"));

    test.endOfLibrary();

    // clean things up
    if (boost::filesystem::exists(file_dir)) { boost::filesystem::remove_all(file_dir); }
  }
}
