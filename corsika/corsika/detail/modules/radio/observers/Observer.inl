/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/observers/Observer.hpp>

namespace corsika {

  template <typename TObserverImpl>
  inline Observer<TObserverImpl>::Observer(std::string const& name, Point const& location,
                                           CoordinateSystem const& coordinateSystem)
      : name_(name)
      , location_(location)
      , coordinateSystem_(coordinateSystem) {}

  template <typename TObserverImpl>
  inline Point const& Observer<TObserverImpl>::getLocation() const {
    return location_;
  }

  template <typename TObserverImpl>
  inline std::string const& Observer<TObserverImpl>::getName() const {
    return name_;
  }

  template <typename TObserverImpl>
  inline TObserverImpl& Observer<TObserverImpl>::implementation() {
    return static_cast<TObserverImpl&>(*this);
  }

} // namespace corsika
