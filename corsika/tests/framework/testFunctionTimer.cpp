/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/analytics/FunctionTimer.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <catch2/catch_all.hpp>

#include <chrono>
#include <thread>

using namespace corsika;

int testFunc() {
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  return 31415;
}

class TestClass {
public:
  int operator()() {
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    return 31415;
  }
};

TEST_CASE("FunctionTimer", "[Timer]") {

  logging::set_level(logging::level::info);

  SECTION("Measure runtime of a free function") {

    auto test = corsika::FunctionTimer(testFunc);

    CORSIKA_LOG_DEBUG(test());
    CORSIKA_LOG_DEBUG(test.getTime().count());
  }

  SECTION("Measure runtime of a class functor") {
    TestClass testC;
    auto test = corsika::FunctionTimer(testC);

    CORSIKA_LOG_DEBUG(test());
    CORSIKA_LOG_DEBUG(test.getTime().count());
  }
}
