/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/process/NullModel.hpp>

using namespace corsika;

TEST_CASE("NullModel", "[processes]") {

  logging::set_level(logging::level::info);

  SECTION("interface") {
    [[maybe_unused]] NullModel model; // nothing to test...
  }
}
