/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <boost/filesystem.hpp>

#include <corsika/modules/writers/WriterOff.hpp>

#include <corsika/framework/core/Logging.hpp>

using namespace corsika;
using Catch::Approx;

/*
class TestEnergyLoss : public corsika::EnergyLossWriter<> {
public:
  TestEnergyLoss(corsika::media::ShowerAxis const& axis)
      : EnergyLossWriter(axis) {}

  YAML::Node getConfig() const { return YAML::Node(); }
};
*/

TEST_CASE("WriterOff") {

  logging::set_level(logging::level::info);

  WriterOff test{"irrelevant", 3};
  WriterOff test2{};

  test.startOfLibrary("./output_dir_eloss");
  test.startOfShower(0);
  test.endOfShower(0);
  test.endOfLibrary();

  auto const config = test.getConfig();

  auto const summary = test.getSummary();
}
