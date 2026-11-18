/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <iterator>
#include <random>
#include <sstream>
#include <tuple>
#include <utility>

namespace corsika {

  template <typename CBPRNG>
  inline void RNGManager<CBPRNG>::registerRandomStream(string_type const& pStreamName) {

    auto const& it = rngs_.find(pStreamName);

    if (it == rngs_.end()) // key not in container, so create one and initialize the value
      rngs_.emplace(std::piecewise_construct, std::forward_as_tuple(pStreamName.c_str()),
                    std::forward_as_tuple(seed_, uint32_t(rngs_.size())));
  }

  template <typename CBPRNG>
  inline typename RNGManager<CBPRNG>::prng_type& RNGManager<CBPRNG>::getRandomStream(
      string_type const& pStreamName) {

    if (isRegistered(pStreamName)) {
      return rngs_.at(pStreamName);
    } else { // this stream name is not in the map
      throw std::runtime_error(
          "'" + pStreamName +
          "' is not a registered stream. Ensure that you have called\n"
          " RNGManager<>::getInstance().registerRandomStream(\"" +
          pStreamName +
          "\");\nbefore initializing instances of the class "
          "that need this stream.");
    }
  }

  template <typename CBPRNG>
  inline bool RNGManager<CBPRNG>::isRegistered(string_type const& pStreamName) const {
    return rngs_.count(pStreamName) > 0;
  }

  template <typename CBPRNG>
  inline std::stringstream RNGManager<CBPRNG>::dumpState() const {
    std::stringstream buffer;
    for (auto const& [streamName, rng] : rngs_) {
      buffer << '"' << streamName << "\" = \"" << rng << '"' << std::endl;
    }

    return buffer;
  }

} // namespace corsika
