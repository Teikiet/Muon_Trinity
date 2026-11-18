/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/stack/Stack.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  template <template <typename> class TParticleInterfaceA,
            template <typename> class TParticleInterfaceB, typename TStackIterator>
  template <typename... TArgs1>
  inline void CombinedParticleInterface<
      TParticleInterfaceA, TParticleInterfaceB,
      TStackIterator>::setParticleData(std::tuple<TArgs1...> const vA) {
    pi_a_type::setParticleData(vA);
    pi_b_type::setParticleData();
  }

  template <template <typename> class TParticleInterfaceA,
            template <typename> class TParticleInterfaceB, typename TStackIterator>
  template <typename... TArgs1, typename... TArgs2>
  inline void CombinedParticleInterface<
      TParticleInterfaceA, TParticleInterfaceB,
      TStackIterator>::setParticleData(std::tuple<TArgs1...> const vA,
                                       std::tuple<TArgs2...> const vB) {
    pi_a_type::setParticleData(vA);
    pi_b_type::setParticleData(vB);
  }

  template <template <typename> class TParticleInterfaceA,
            template <typename> class TParticleInterfaceB, typename TStackIterator>
  template <typename... TArgs1>
  inline void CombinedParticleInterface<
      TParticleInterfaceA, TParticleInterfaceB,
      TStackIterator>::setParticleData(pi_a_type& p, std::tuple<TArgs1...> const vA) {
    // static_assert(MT<I>::has_not, "error");
    pi_a_type::setParticleData(static_cast<pi_a_type&>(p), vA); // original stack
    pi_b_type::setParticleData(static_cast<pi_b_type&>(p));     // addon stack
  }

  template <template <typename> class TParticleInterfaceA,
            template <typename> class TParticleInterfaceB, typename TStackIterator>
  template <typename... TArgs1, typename... TArgs2>
  inline void CombinedParticleInterface<
      TParticleInterfaceA, TParticleInterfaceB,
      TStackIterator>::setParticleData(pi_c_type& p, std::tuple<TArgs1...> const vA,
                                       std::tuple<TArgs2...> const vB) {

    pi_a_type::setParticleData(static_cast<pi_a_type&>(p), vA);
    pi_b_type::setParticleData(static_cast<pi_b_type&>(p), vB);
  }

  template <template <typename> class TParticleInterfaceA,
            template <typename> class TParticleInterfaceB, typename TStackIterator>
  inline std::string CombinedParticleInterface<TParticleInterfaceA, TParticleInterfaceB,
                                               TStackIterator>::asString() const {
    return fmt::format("[[{}][{}]]", pi_a_type::asString(), pi_b_type::asString());
  }

  template <typename Stack1Impl, typename Stack2Impl>
  inline void CombinedStackImpl<Stack1Impl, Stack2Impl>::clear() {
    Stack1Impl::clear();
    Stack2Impl::clear();
  }
  template <typename Stack1Impl, typename Stack2Impl>
  inline void CombinedStackImpl<Stack1Impl, Stack2Impl>::copy(const unsigned int i1,
                                                              const unsigned int i2) {
    if (i1 >= getSize() || i2 >= getSize()) {
      std::ostringstream err;
      err << "CombinedStack: trying to access data beyond size of stack!";
      throw std::runtime_error(err.str());
    }
    Stack1Impl::copy(i1, i2);
    Stack2Impl::copy(i1, i2);
  }

  template <typename Stack1Impl, typename Stack2Impl>
  inline void CombinedStackImpl<Stack1Impl, Stack2Impl>::swap(const unsigned int i1,
                                                              const unsigned int i2) {
    if (i1 >= getSize() || i2 >= getSize()) {
      std::ostringstream err;
      err << "CombinedStack: trying to access data beyond size of stack!";
      throw std::runtime_error(err.str());
    }
    Stack1Impl::swap(i1, i2);
    Stack2Impl::swap(i1, i2);
  }

  template <typename Stack1Impl, typename Stack2Impl>
  inline void CombinedStackImpl<Stack1Impl, Stack2Impl>::incrementSize() {
    Stack1Impl::incrementSize();
    Stack2Impl::incrementSize();
  }

  template <typename Stack1Impl, typename Stack2Impl>
  inline void CombinedStackImpl<Stack1Impl, Stack2Impl>::decrementSize() {
    Stack1Impl::decrementSize();
    Stack2Impl::decrementSize();
  }

} // namespace corsika
