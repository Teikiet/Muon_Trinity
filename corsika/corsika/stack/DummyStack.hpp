/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/stack/Stack.hpp>

#include <string>
#include <tuple>

namespace corsika::dummy_stack {

  /**
   * Example of a particle object on the stack, with NO DATA.
   */

  /**
   * However, conceptually we need to provide fake data. A stack without data does not
   * work...
   */

  struct NoData { /* nothing */
    int nothing = 0;
  };

  template <typename StackIteratorInterface>
  struct ParticleInterface : public corsika::ParticleBase<StackIteratorInterface> {

    typedef corsika::ParticleBase<StackIteratorInterface> super_type;

  public:
    void setParticleData(const std::tuple<NoData>& /*v*/) {}
    void setParticleData(super_type& /*parent*/, const std::tuple<NoData>& /*v*/) {}

    std::string asString() const { return "dummy-data"; }

    // unfortunately we need those dummy getter
    // for some more complex tests with "history"
    HEPEnergyType getEnergy() const { return 0_GeV; }
    MomentumVector getMomentum() const {
      return MomentumVector(get_root_CoordinateSystem(), {0_GeV, 0_GeV, 0_GeV});
    }
    Code getPID() const { return Code::Unknown; }
  };

  /**
   *
   * Memory implementation of the most simple (no-data) particle stack object.
   */

  class DummyStackImpl {

  public:
    DummyStackImpl() = default;

    DummyStackImpl(DummyStackImpl const&) = default;

    DummyStackImpl(DummyStackImpl&&) = default;

    DummyStackImpl& operator=(DummyStackImpl const&) = default;
    DummyStackImpl& operator=(DummyStackImpl&&) = default;

    void init() { entries_ = 0; }

    void clear() { entries_ = 0; }

    int getSize() const { return entries_; }
    int getCapacity() const { return entries_; }

    /**
     *   Function to copy particle at location i2 in stack to i1.
     */
    void copy(const int /*i1*/, const int /*i2*/) {}

    void swap(const int, const int) {}

    void incrementSize() { entries_++; }
    void decrementSize() { entries_--; }

    int getEntries() const { return entries_; }
    void setEntries(int entries = 0) { entries_ = entries; }

  private:
    int entries_ = 0;

  }; // end class DummyStackImpl

  typedef Stack<DummyStackImpl, ParticleInterface> DummyStack;

} // namespace corsika::dummy_stack
