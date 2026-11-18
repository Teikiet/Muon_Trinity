/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/stack/Stack.hpp>
#include <tuple>
#include <vector>

/**
 * definition of stack-data object for unit tests
 *
 * TestStackData contain only a single variable "Data" stored in fData
 * with get/setData functions.
 */
class TestStackData {

public:
  // these functions are needed for the Stack interface
  void clear() { data_.clear(); }
  unsigned int getSize() const { return data_.size(); }
  unsigned int getCapacity() const { return data_.size(); }
  void copy(const unsigned int i1, const unsigned int i2) { data_[i2] = data_[i1]; }
  void swap(const unsigned int i1, const unsigned int i2) {
    double tmp0 = data_[i1];
    data_[i1] = data_[i2];
    data_[i2] = tmp0;
  }

  // custom data access function
  void setData(const unsigned int i, const double v) { data_[i] = v; }
  double getData(const unsigned int i) const { return data_[i]; }

  // these functions are also needed by the Stack interface
  void incrementSize() { data_.push_back(0.); }
  void decrementSize() {
    if (data_.size() > 0) { data_.pop_back(); }
  }

  // custom private data section
private:
  std::vector<double> data_;
};

/**
 * From static_cast of a StackIterator over entries in the
 * TestStackData class you get and object of type
 * TestParticleInterface defined here
 *
 * It provides Get/Set methods to read and write data to the "Data"
 * storage of TestStackData obtained via
 * "StackIteratorInterface::getStackData()", given the index of the
 * iterator "StackIteratorInterface::getIndex()"
 *
 */
template <typename StackIteratorInterface>
class TestParticleInterface : public corsika::ParticleBase<StackIteratorInterface> {

  typedef corsika::ParticleBase<StackIteratorInterface> super_type;

public:
  /*
     The SetParticleData methods are called for creating new entries
     on the stack. You can specifiy various parametric versions to
     perform this task:
  */

  // default version for particle-creation from input data
  void setParticleData(std::tuple<double> v) { setData(std::get<0>(v)); }
  void setParticleData(TestParticleInterface<StackIteratorInterface>& /*parent*/,
                       std::tuple<double> v) {
    setData(std::get<0>(v));
  }

  // here are the fundamental methods for access to TestStackData data
  void setData(const double v) {
    super_type::getStackData().setData(super_type::getIndex(), v);
  }
  double getData() const {
    return super_type::getStackData().getData(super_type::getIndex());
  }
};
