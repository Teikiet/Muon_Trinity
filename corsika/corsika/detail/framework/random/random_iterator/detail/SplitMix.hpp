/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * SplitMix.hpp
 *
 *  Created on: 25/02/2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

namespace random_iterator {

  namespace detail {

    template <typename UIntType>
    inline UIntType splitmix(UIntType&);

    template <>
    inline uint32_t splitmix<uint32_t>(uint32_t& x) {
      uint32_t z = (x += 0x6D2B79F5UL);
      z = (z ^ (z >> 15)) * (z | 1UL);
      z ^= z + (z ^ (z >> 7)) * (z | 61UL);
      return z ^ (z >> 14);
    }

    template <>
    inline uint64_t splitmix<uint64_t>(uint64_t& x) {
      uint64_t z = (x += 0x9e3779b97f4a7c15);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
      z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
      return z ^ (z >> 31);
    }

  } // namespace detail

} // namespace random_iterator
