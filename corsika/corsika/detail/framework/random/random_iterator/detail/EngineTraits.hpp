/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * EngineTraits.hpp
 *
 *  Created on: 23 de fev. de 2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

#include <stdint.h>
#include "Random123/array.h"
#include "Random123/philox.h"
#include "Random123/threefry.h"
#include "Random123/ars.h"
#include "Random123/ReinterpretCtr.hpp"

namespace random_iterator {

  namespace detail {

    template <typename Engine>
    struct random_traits;
    /*
     * random_traits<T>::state_type { counter, state}
     * random_traits<T>::advance_type;
     * random_traits<T>::init_type;
     * random_traits<T>::result_type;
     */

    // philox
    template <>
    struct random_traits<random_iterator_r123::Philox2x64> {
      typedef typename random_iterator_r123::Philox2x64::ctr_type state_type;
      typedef typename random_iterator_r123::Philox2x64::key_type seed_type;
      typedef typename random_iterator_r123::Philox2x64::ukey_type init_type;
      typedef uint64_t advance_type;
      typedef state_type::value_type result_type;

      enum { arity = 2 };
    };

    template <>
    struct random_traits<random_iterator_r123::Philox4x64> {
      typedef typename random_iterator_r123::Philox4x64::ctr_type state_type;
      typedef typename random_iterator_r123::Philox4x64::key_type seed_type;
      typedef typename random_iterator_r123::Philox4x64::ukey_type init_type;
      typedef uint64_t advance_type;
      typedef state_type::value_type result_type;

      enum { arity = 4 };
    };
    //
    template <>
    struct random_traits<random_iterator_r123::Threefry2x64> {
      typedef typename random_iterator_r123::Threefry2x64::ctr_type state_type;
      typedef typename random_iterator_r123::Threefry2x64::key_type seed_type;
      typedef typename random_iterator_r123::Threefry2x64::ukey_type init_type;
      typedef uint64_t advance_type;
      typedef state_type::value_type result_type;

      enum { arity = 2 };
    };

    //
    template <>
    struct random_traits<random_iterator_r123::Threefry4x64> {
      typedef typename random_iterator_r123::Threefry4x64::ctr_type state_type;
      typedef typename random_iterator_r123::Threefry4x64::key_type seed_type;
      typedef typename random_iterator_r123::Threefry4x64::ukey_type init_type;
      typedef uint64_t advance_type;
      typedef state_type::value_type result_type;

      enum { arity = 4 };
    };

#if RANDOM_ITERATOR_R123_USE_AES_NI
    template <>
    struct random_traits<random_iterator_r123::ARS4x32> {
      typedef typename random_iterator_r123::ARS4x32::ctr_type state_type;
      typedef typename random_iterator_r123::ARS4x32::key_type seed_type;
      typedef typename random_iterator_r123::ARS4x32::ukey_type init_type;
      typedef uint64_t advance_type;
      typedef state_type::value_type result_type;

      enum { arity = 4 };
    };

    template <>
    struct random_traits<random_iterator_r123::ARS2x64> {
      typedef typename random_iterator_r123::ARS2x64::ctr_type state_type;
      typedef typename random_iterator_r123::ARS2x64::key_type seed_type;
      typedef typename random_iterator_r123::ARS2x64::ukey_type init_type;
      typedef uint64_t advance_type;
      typedef state_type::value_type result_type;

      enum { arity = 2 };
    };

#endif

  } // namespace detail

} // namespace random_iterator
