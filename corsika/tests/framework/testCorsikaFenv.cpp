/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <cmath>
#include <csignal>
#include <iostream>

int gRESULT = 1;

extern "C" {
static void handle_fpe(int /*signo*/) { gRESULT = 0; }
}

TEST_CASE("CorsikaFenv", "[fenv]") {

  logging::set_level(logging::level::info);

  SECTION("Enable all exceptions") { feenableexcept(FE_ALL_EXCEPT); }

  signal(SIGFPE, handle_fpe);

  SECTION("exception") {

    [[maybe_unused]] auto trigger = std::log(0.);
    std::cout << "trigger: " << trigger << std::endl;
    CHECK(gRESULT == 0);
  }
}
}
