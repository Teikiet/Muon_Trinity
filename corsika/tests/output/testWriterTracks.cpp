/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <boost/filesystem.hpp>

#include <corsika/modules/writers/TrackWriterParquet.hpp>

#include <corsika/framework/core/Logging.hpp>
#include <SetupTestEnvironment.hpp>

using namespace corsika;
using Catch::Approx;

struct TestWriterTrack : public TrackWriterParquet {

  YAML::Node getConfig() const { return YAML::Node(); }

  void checkWrite() {
    TrackWriterParquet::write(Code::Unknown, 1_eV, 1.0, {2_m, 3_m, 4_m}, 1_ns,
                              {5_m, 6_m, 7_m}, 0_eV, 2_ns);
  }
};

TEST_CASE("TrackWriterParquet") {

  logging::set_level(logging::level::info);

  SECTION("standard") {

    // preparation
    if (boost::filesystem::exists("./output_dir_tracks")) {
      boost::filesystem::remove_all("./output_dir_tracks");
    }
    boost::filesystem::create_directory("./output_dir_tracks");

    TestWriterTrack test;
    test.startOfLibrary("./output_dir_tracks");
    test.startOfShower(0);
    test.checkWrite();
    test.endOfShower(0);
    test.endOfLibrary();

    CHECK(boost::filesystem::exists("./output_dir_tracks/tracks.parquet"));
  }
}
