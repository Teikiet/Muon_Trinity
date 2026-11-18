/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

/**
 * \file CoordinateSystemRegistry.inl
 *
 * Inline implementations for CoordinateSystemRegistry.
 **/

namespace corsika {

  inline CoordinateSystemRegistry::CoordinateSystemRegistry()
      : effectiveOrigin_(Eigen::Vector3d::Zero()) {
    // Create root CS at index 0 with identity transform
    EigenTransform id = EigenTransform::Identity();
    systems_.emplace_back(nullptr, id);
  }

  inline CoordinateSystemRegistry& CoordinateSystemRegistry::getInstance() {
    static CoordinateSystemRegistry instance;
    return instance;
  }

  inline const CoordinateSystem& CoordinateSystemRegistry::registerCS(
      const CoordinateSystem* parent, EigenTransform const& transform) {

    if (parent == nullptr) {
      throw std::runtime_error("Cannot register coordinate system with nullptr parent.");
    }

    std::unique_lock<std::shared_mutex> lock(mutex_);

    systems_.emplace_back(parent, transform);

    return systems_.back();
  }

  inline const CoordinateSystem& CoordinateSystemRegistry::getByIndex(
      size_t index) const {
    if (index >= systems_.size()) {
      throw std::out_of_range("CoordinateSystem index out of range");
    }
    return systems_[index];
  }

  inline const CoordinateSystem& CoordinateSystemRegistry::getRoot() const {
    return systems_[0];
  }

  inline void CoordinateSystemRegistry::setEffectiveOrigin(
      Eigen::Vector3d const& origin) {
    std::unique_lock<std::shared_mutex> lock(originMutex_);
    effectiveOrigin_ = origin;
  }

  inline Eigen::Vector3d CoordinateSystemRegistry::getEffectiveOrigin() const {
    std::shared_lock<std::shared_mutex> lock(originMutex_);
    return effectiveOrigin_;
  }

  inline size_t CoordinateSystemRegistry::size() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return systems_.size();
  }

  inline const CoordinateSystem& get_root_CoordinateSystem() {
    return CoordinateSystemRegistry::getInstance().getRoot();
  }

  // CoordinateSystem constructors (need to be defined after Registry)
  inline CoordinateSystem::CoordinateSystem(const CoordinateSystem* parent,
                                            EigenTransform const& transf)
      : parent_(parent)
      , transformToParent_(transf) {
    // Cache will be computed on first use (lazy initialization)
  }

  inline CoordinateSystem::CoordinateSystem()
      : parent_(nullptr)
      , transformToParent_(EigenTransform::Identity()) {
    // Root CS: identity transformation
  }

} // namespace corsika