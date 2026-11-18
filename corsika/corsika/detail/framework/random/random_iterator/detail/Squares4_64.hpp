/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * Squares4.hpp
 *
 *  Created on: 25/02/2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

#include <stdint.h>
#include "SquaresKeys.hpp"

namespace random_iterator {

  namespace detail {

    class Squares4_64 {

    public:
      typedef uint64_t init_type;
      typedef uint64_t state_type;
      typedef uint64_t seed_type;
      typedef uint64_t advance_type;
      typedef uint32_t result_type;

      Squares4_64() = delete;

      Squares4_64(seed_type s, uint32_t)
          : state_(0)
          , seed_(seed_type{splitmix<seed_type>(s)}) {}

      Squares4_64(Squares4_64 const& other)
          : state_(other.getState())
          , seed_(other.getSeed()) {}

      inline Squares4_64& operator=(Squares4_64 const& other) {
        if (this == &other) return *this;

        state_ = other.getState();
        seed_ = other.getSeed();

        return *this;
      }

      inline result_type operator()(void) {
        uint64_t x, y, z;

        y = x = seed_ * state_;
        z = y + seed_;

        x = x * x + y;
        x = (x >> 32) | (x << 32); /* round 1 */

        x = x * x + z;
        x = (x >> 32) | (x << 32); /* round 2 */

        x = x * x + y;
        x = (x >> 32) | (x << 32); /* round 3 */

        ++state_; /* advance state */

        return (x * x + z) >> 32; /* round 4 */
      }

      inline void reset(void) { state_ = 0; }

      inline void discard(advance_type n) { state_ += n; }

      inline seed_type getSeed() const { return seed_; }

      inline void setSeed(seed_type seed) { seed_ = seed; }

      inline state_type getState() const { return state_; }

      inline void setState(state_type state) { state_ = state; }

      inline static uint64_t generateSeed(size_t i) { return keys[i]; }

      static constexpr result_type min() { return 0; }

      static constexpr result_type max() {
        return std::numeric_limits<result_type>::max();
      }

      friend inline std::ostream& operator<<(std::ostream& os, const Squares4_64& be) {
        return os << "state: " << be.getState() << " seed: " << be.getSeed();
      }

    private:
      state_type state_;
      seed_type seed_;
    };

  } // namespace detail

} // namespace random_iterator
