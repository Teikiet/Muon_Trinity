/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <functional>
#include <random>
#include <string_view>

#include <corsika/framework/random/RNGManager.hpp>

namespace corsika {
  using rng_function_type = std::function<void(double*, std::size_t)>;

  namespace detail {
    inline void rng_func(corsika::default_prng_type& rng, double* dest, std::size_t N) {
      std::uniform_real_distribution<double> udist(0.0, 1.0);
      std::generate(dest, std::next(dest, N), std::bind(udist, std::ref(rng)));
    }
  } // namespace detail

  inline void connect_random_stream(corsika::default_prng_type& rng,
                                    void (*injection_func)(rng_function_type)) {
    using namespace std::placeholders;
    injection_func(std::bind(detail::rng_func, rng, _1, _2));
  }

  inline void connect_random_stream(std::string_view stream_name,
                                    void (*injection_func)(rng_function_type)) {
    auto& rng = RNGManager<>::getInstance().getRandomStream(std::string{stream_name});
    connect_random_stream(rng, injection_func);
  }
} // namespace corsika
