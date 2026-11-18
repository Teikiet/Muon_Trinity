/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

namespace corsika {

  inline DummyOutputManager::DummyOutputManager() {}

  inline DummyOutputManager::~DummyOutputManager() {}

  template <typename TOutput>
  inline void DummyOutputManager::add(std::string const&, TOutput&) {}

  inline void DummyOutputManager::startOfLibrary() {}

  inline void DummyOutputManager::startOfShower() {}

  inline void DummyOutputManager::endOfShower() {}

  inline void DummyOutputManager::endOfLibrary() {}

} // namespace corsika
