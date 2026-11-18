/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * EngineCaller.hpp
 *
 *  Created on: 23/02/2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

#include <stdint.h>

namespace random_iterator {

  namespace detail {

    template <typename DistibutionType, typename EngineType>
    struct EngineCaller {
      typedef EngineType engine_type;
      typedef typename engine_type::state_type state_type;
      typedef typename engine_type::seed_type seed_type;
      typedef typename engine_type::advance_type advance_type;
      typedef typename engine_type::init_type init_type;

      typedef DistibutionType distribution_type;
      typedef typename distribution_type::result_type result_type;

      EngineCaller() = delete;

      EngineCaller(distribution_type const& dist, seed_type seed, uint32_t stream)
          : distribution_(dist)
          , seed_(seed)
          , stream_(stream) {}

      EngineCaller(EngineCaller<DistibutionType, EngineType> const& other)
          : distribution_(other.distribution_)
          , seed_(other.seed_)
          , stream_(other.stream_) {}

      EngineCaller<DistibutionType, EngineType>& operator=(
          EngineCaller<DistibutionType, EngineType> const& other) {

        if (this == &other) return *this;

        distribution_ = other.distribution_;
        seed_ = other.seed_;
        stream_ = other.stream_;

        return *this;
      }

      inline result_type operator()(advance_type n) const {

        EngineType eng(seed_, stream_);
        eng.discard(n);

        return static_cast<distribution_type>(distribution_)(eng);
      }

      distribution_type distribution_;
      seed_type seed_;
      uint32_t stream_;
    };

  } // namespace detail

} // namespace random_iterator
