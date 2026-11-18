/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

/**
 * \file CoordinateSystemRegistry.hpp
 *
 * Singleton registry for coordinate system management.
 **/

#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <vector>
#include <shared_mutex>
#include <Eigen/Dense>
#include <deque>

namespace corsika {

  /**
   * Singleton registry for coordinate system management.
   *
   * All coordinate systems are stored in a contiguous vector and accessed by const
   * reference. Coordinate systems are anonymous with auto-generated identifiers for
   * debugging purposes only.
   *
   * Thread-safety: concurrent reads are allowed via shared_mutex, modifications require
   * exclusive lock.
   */
  class CoordinateSystemRegistry {
  public:
    /**
     * Get the singleton instance.
     */
    static CoordinateSystemRegistry& getInstance();

    /**
     * Register a new coordinate system.
     *
     * @param parent Parent coordinate system (must not be nullptr)
     * @param transform Transformation from this coordinate system to parent
     * @return const reference to the registered coordinate system
     *
     * Thread-safe. The counter for generating identifiers is protected by the same
     * mutex used for modifying the storage vector.
     */
    const CoordinateSystem& registerCS(const CoordinateSystem* parent,
                                       EigenTransform const& transform);

    /**
     * Get the root coordinate system.
     *
     * The root coordinate system is always at index 0 with identity transform and
     * nullptr parent. No locking required.
     */
    const CoordinateSystem& getRoot() const;

    /**
     * Get a coordinate system by index.
     *
     * @param index Index in the registry (0 is root)
     * @return const reference to coordinate system
     * @throws std::out_of_range if index is invalid
     *
     * No locking required after initialization.
     */
    const CoordinateSystem& getByIndex(size_t index) const;

    /**
     * Get number of registered coordinate systems.
     */
    size_t size() const;

    /**
     * Set effective origin for numerical precision.
     *
     * @param origin Origin point for coordinate storage
     *
     * Vector coordinates are stored relative to this origin to maintain precision
     * when working at large distances from the global origin.
     */
    void setEffectiveOrigin(Eigen::Vector3d const& origin);

    /**
     * Get current effective origin.
     */
    Eigen::Vector3d getEffectiveOrigin() const;

    // Singleton: delete copy/move
    CoordinateSystemRegistry(const CoordinateSystemRegistry&) = delete;
    CoordinateSystemRegistry& operator=(const CoordinateSystemRegistry&) = delete;
    CoordinateSystemRegistry(CoordinateSystemRegistry&&) = delete;
    CoordinateSystemRegistry& operator=(CoordinateSystemRegistry&&) = delete;

  private:
    /**
     * Private constructor creates root coordinate system at index 0.
     */
    CoordinateSystemRegistry();

    // Storage in a deque to ensure that all references are valid even after reallocation
    // (vector may invalidate references)
    std::deque<CoordinateSystem, Eigen::aligned_allocator<CoordinateSystem>> systems_;

    // Thread-safety
    mutable std::shared_mutex mutex_;

    // Effective origin
    Eigen::Vector3d effectiveOrigin_;
    mutable std::shared_mutex originMutex_;
  };

  /**
   * Convenience function to access root coordinate system.
   */
  inline const CoordinateSystem& get_root_CoordinateSystem();

} // namespace corsika

#include <corsika/detail/framework/geometry/CoordinateSystemRegistry.inl>