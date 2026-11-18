/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <random>

namespace bh = boost::histogram;

TEST_CASE("SaveHistogram") {

  corsika::logging::set_level(corsika::logging::level::info);

  std::mt19937 rng;
  std::normal_distribution n1{2., 6.};
  std::exponential_distribution e{1.};
  std::uniform_int_distribution u{1, 10};
  std::uniform_real_distribution<double> r{-3, 3};

  auto h = bh::make_histogram(
      bh::axis::regular{5, 0., 10., "normal"},
      bh::axis::regular{3, 0., 4., "exponential"},
      bh::axis::category<int>{{2, 3, 5, 7}, "integer category"},
      bh::axis::regular<double, bh::use_default, bh::use_default,
                        bh::axis::option::growth_t>{10, -1, 1, "integer category"});

  for (int i{0}; i < 100'000; ++i) {
    auto const a = n1(rng);
    auto const b = e(rng);
    auto const c = u(rng);
    auto const d = r(rng);

    h(a, b, c, d);
  }

  REQUIRE_NOTHROW(corsika::save_hist(h, "hist.npz", true));
  REQUIRE_THROWS(corsika::save_hist(h, "hist.npz", false));
}
