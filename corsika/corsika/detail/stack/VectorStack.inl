/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/stack/Stack.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>

#include <corsika/framework/geometry/Vector.hpp>

#include <string>
#include <tuple>
#include <vector>

namespace corsika {

  template <typename StackIteratorInterface>
  inline void ParticleInterface<StackIteratorInterface>::setParticleData(
      particle_data_type const& v) {
    this->setPID(std::get<0>(v));
    this->setKineticEnergy(std::get<1>(v));
    this->setDirection(std::get<2>(v));
    this->setPosition(std::get<3>(v));
    this->setTime(std::get<4>(v));
  }

  template <typename StackIteratorInterface>
  inline void ParticleInterface<StackIteratorInterface>::setParticleData(
      ParticleInterface<StackIteratorInterface> const& parent,
      secondary_data_type const& v) {
    this->setPID(std::get<0>(v));
    this->setKineticEnergy(std::get<1>(v));
    this->setDirection(std::get<2>(v));
    this->setPosition(parent.getPosition()); // position
    this->setTime(parent.getTime());         // parent time
  }

  template <typename StackIteratorInterface>
  inline void ParticleInterface<StackIteratorInterface>::setParticleData(
      ParticleInterface<StackIteratorInterface> const& parent,
      secondary_extended_data_type const& v) {
    this->setPID(std::get<0>(v));
    this->setKineticEnergy(std::get<1>(v));
    this->setDirection(std::get<2>(v));
    this->setPosition(parent.getPosition() + std::get<3>(v)); // + position
    this->setTime(parent.getTime() + std::get<4>(v));         // + parent time
  }

  template <typename StackIteratorInterface>
  inline std::string ParticleInterface<StackIteratorInterface>::asString() const {
    return fmt::format("particle: i={}, PID={}, Ekin={}GeV", super_type::getIndex(),
                       get_name(this->getPID()), this->getKineticEnergy() / 1_GeV);
  }

  inline void VectorStackImpl::clear() {
    dataPID_.clear();
    dataEkin_.clear();
    direction_.clear();
    position_.clear();
    time_.clear();
  }

  inline void VectorStackImpl::copy(size_t const i1, size_t const i2) {
    // index range check
    if (i1 >= getSize() || i2 >= getSize()) {
      std::ostringstream err;
      err << "VectorStackImpl: trying to access data beyond size of stack !";
      throw std::runtime_error(err.str());
    }
    dataPID_[i2] = dataPID_[i1];
    dataEkin_[i2] = dataEkin_[i1];
    direction_[i2] = direction_[i1];
    position_[i2] = position_[i1];
    time_[i2] = time_[i1];
  }

  inline void VectorStackImpl::swap(size_t const i1, size_t const i2) {
    // index range check
    if (i1 >= getSize() || i2 >= getSize()) {
      std::ostringstream err;
      err << "VectorStackImpl: trying to access data beyond size of stack !";
      throw std::runtime_error(err.str());
    }
    std::swap(dataPID_[i2], dataPID_[i1]);
    std::swap(dataEkin_[i2], dataEkin_[i1]);
    std::swap(direction_[i2], direction_[i1]);
    std::swap(position_[i2], position_[i1]);
    std::swap(time_[i2], time_[i1]);
  }

  inline void VectorStackImpl::incrementSize() {
    dataPID_.push_back(Code::Unknown);
    dataEkin_.push_back(0 * electronvolt);

    direction_.push_back(DirectionVector({0, 0, 0}));

    position_.push_back(Point(0 * meter, 0 * meter, 0 * meter));
    time_.push_back(0 * second);
  }

  inline void VectorStackImpl::decrementSize() {
    if (dataEkin_.size() > 0) {
      dataPID_.pop_back();
      dataEkin_.pop_back();
      direction_.pop_back();
      position_.pop_back();
      time_.pop_back();
    }
  }

} // namespace corsika
