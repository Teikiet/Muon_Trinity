/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/stack/Stack.hpp>

#include <tuple>
#include <utility>
#include <vector>

namespace corsika::weights {

  // default version for particle-creation from input data
  template <typename TParentStack>
  inline void WeightDataInterface<TParentStack>::setParticleData(
      std::tuple<double> const v) {
    setWeight(std::get<0>(v));
  }

  template <typename TParentStack>
  inline void WeightDataInterface<TParentStack>::setParticleData(
      WeightDataInterface const& parent, std::tuple<double> const) {
    setWeight(parent.getWeight()); // copy Weight from parent particle!
  }

  template <typename TParentStack>
  inline void WeightDataInterface<TParentStack>::setParticleData() {
    setWeight(1);
  } // default weight

  template <typename TParentStack>
  inline void WeightDataInterface<TParentStack>::setParticleData(
      WeightDataInterface const& parent) {
    setWeight(parent.getWeight()); // copy Weight from parent particle!
  }

  template <typename TParentStack>
  inline std::string WeightDataInterface<TParentStack>::asString() const {
    return fmt::format("weight={}", getWeight());
  }

  template <typename TParentStack>
  inline void WeightDataInterface<TParentStack>::setWeight(double const v) {
    super_type::getStackData().setWeight(super_type::getIndex(), v);
  }

  template <typename TParentStack>
  inline double WeightDataInterface<TParentStack>::getWeight() const {
    return super_type::getStackData().getWeight(super_type::getIndex());
  }

  // definition of stack-data object to store geometry information

  // these functions are needed for the Stack interface
  inline void WeightData::clear() { weight_vector_.clear(); }

  inline unsigned int WeightData::getSize() const { return weight_vector_.size(); }

  inline unsigned int WeightData::getCapacity() const { return weight_vector_.size(); }

  inline void WeightData::copy(int const i1, int const i2) {
    weight_vector_[i2] = weight_vector_[i1];
  }

  inline void WeightData::swap(int const i1, int const i2) {
    std::swap(weight_vector_[i1], weight_vector_[i2]);
  }

  // custom data access function
  inline void WeightData::setWeight(int const i, double const v) {
    weight_vector_[i] = v;
  }

  inline double WeightData::getWeight(int const i) const { return weight_vector_[i]; }

  // these functions are also needed by the Stack interface
  inline void WeightData::incrementSize() {
    weight_vector_.push_back(1);
  } // default weight

  inline void WeightData::decrementSize() {
    if (weight_vector_.size() > 0) { weight_vector_.pop_back(); }
  }

} // namespace corsika::weights
