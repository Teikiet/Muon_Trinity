/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/process/NullModel.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;

/*
 * The NullModel can do really nothing, so we can basically test
 * nothing.
 */

TEST_CASE("NullModel", "[processes]") {

  logging::set_level(logging::level::info);

  SECTION("interface") {
    [[maybe_unused]] NullModel nm;

    CHECK(is_process_v<decltype(nm)>);
    CHECK(count_processes<decltype(nm)>::count == 0);
  }
}
