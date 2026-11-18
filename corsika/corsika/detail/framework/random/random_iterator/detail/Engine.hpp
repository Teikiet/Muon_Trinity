/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * Engine.hpp
 *
 *  Created on: 22 de fev. de 2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

#include <limits>

#include "EngineTraits.hpp"
#include "SplitMix.hpp"
#include "Squares3_128.hpp"
#include "Squares4_128.hpp"
#include "Squares3_64.hpp"
#include "Squares4_64.hpp"

namespace random_iterator {

  namespace detail {

    template <typename Engine123>
    class Engine {

      typedef unsigned trigger_type;

    public:
      typedef Engine123 engine_type;
      typedef typename detail::random_traits<engine_type>::state_type state_type;
      typedef typename detail::random_traits<engine_type>::seed_type seed_type;
      typedef typename detail::random_traits<engine_type>::advance_type advance_type;
      typedef typename detail::random_traits<engine_type>::init_type init_type;
      typedef typename detail::random_traits<engine_type>::result_type result_type;

      static const unsigned arity = detail::random_traits<engine_type>::arity;

      Engine() = delete;

      Engine(result_type seed)
          : engine_(engine_type{})
          , cache_(state_type{})
          , state_(state_type{})
          , seed_(seed_type{})
          , trigger_(arity) {
        init_type temp{{}};

        for (unsigned i = 0; i < temp.size(); ++i) {
          temp[i] = detail::splitmix<result_type>(seed);
        }

        seed_ = temp;
      }

      Engine(result_type seed, uint32_t stream)
          : engine_(engine_type{})
          , cache_(state_type{})
          , state_(state_type{})
          , seed_(seed_type{})
          , trigger_(arity) {
        init_type temp{{}};

        for (unsigned i = 0; i < temp.size(); ++i) {
          temp[i] = detail::splitmix<result_type>(seed);
        }

        state_[arity - 1] = stream;
        state_[arity - 1] = state_[arity - 1] << 32;

        seed_ = temp;
      }

      Engine(init_type seed)
          : engine_(engine_type{})
          , cache_(state_type{})
          , state_(state_type{})
          , seed_(seed)
          , trigger_(arity) {}

      Engine(init_type seed, uint32_t stream)
          : engine_(engine_type{})
          , cache_(state_type{})
          , state_(state_type{})
          , seed_(seed)
          , trigger_(arity) {
        state_[arity - 1] = stream;
        state_[arity - 1] << 32;
      }

      Engine(Engine<Engine123> const& other)
          : engine_(engine_type{})
          , cache_(other.getCache())
          , state_(other.getState())
          , seed_(other.getSeed())
          , trigger_(other.getTrigger()) {}

      inline Engine<Engine123>& operator=(Engine<Engine123> const& other) {

        if (this == &other) return *this;

        engine_ = engine_type{};
        cache_ = other.getCache();
        state_ = other.getState();
        seed_ = other.getSeed();
        trigger_ = other.getTrigger();
        return *this;
      }

      inline result_type operator()(void) {
        result_type result = 0;

        if (trigger_ == arity) {

          trigger_ = 0;
          cache_ = engine_(state_.incr(), seed_);
          result = cache_[trigger_];
          ++trigger_;
        } else {
          result = cache_[trigger_];
          ++trigger_;
        }

        return result;
      }

      inline void discard(advance_type n) {

        state_.incr(n);
        trigger_ = arity;
      }

      inline void reset(void) {
        trigger_ = arity;
        state_ = state_type{};
        cache_ = state_type{};
      }

      inline const seed_type& getSeed() const { return seed_; }

      inline void setSeed(seed_type seed) { seed_ = seed; }

      inline void setSeed(result_type seed) {

        init_type temp{{}};
        for (unsigned i = 0; i < temp.size(); ++i) {
          temp[i] = detail::splitmix<result_type>(seed);
        }
        seed_ = temp;
      }

      inline const state_type& getState() const { return state_; }

      inline void setState(const state_type& state) { state_ = state; }

      inline trigger_type getTrigger() const { return trigger_; }

      inline void setTrigger(trigger_type trigger) { trigger_ = trigger; }
      static constexpr result_type min() { return 0; }

      static constexpr result_type max() {
        return std::numeric_limits<result_type>::max();
      }

      friend inline std::ostream& operator<<(std::ostream& os,
                                             const Engine<Engine123>& be) {
        return os << "state: " << be.getState() << " seed: " << be.getSeed()
                  << " trigger: " << be.getTrigger();
      }

    private:
      inline state_type getCache() const { return cache_; }

      engine_type engine_;
      state_type cache_;
      state_type state_;
      seed_type seed_;
      trigger_type trigger_;
    };

    typedef Engine<random_iterator_r123::Philox2x64> philox;
    typedef Engine<random_iterator_r123::Threefry2x64> threefry;
    typedef detail::Squares3_64 squares3_64;
    typedef detail::Squares4_64 squares4_64;
    typedef detail::Squares3_128 squares3_128;
    typedef detail::Squares4_128 squares4_128;

#if RANDOM_ITERATOR_R123_USE_AES_NI
    typedef Engine<random_iterator_r123::ARS2x64> ars;
#endif

  } // namespace detail

} // namespace random_iterator
