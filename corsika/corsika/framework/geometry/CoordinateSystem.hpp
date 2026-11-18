/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

/**
 * \file CoordinateSystem.hpp
 **/

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <Eigen/Dense>
#include <stdexcept>
#include <memory>

namespace corsika {

  typedef Eigen::Transform<double, 3, Eigen::Affine> EigenTransform;
  typedef Eigen::Translation<double, 3> EigenTranslation;

  template <typename T>
  class Vector; // fwd decl

  class CoordinateSystem; // fwd decl

  /**
   * Get the root coordinate system from the registry.
   */
  const CoordinateSystem& get_root_CoordinateSystem();

  /**
   * Creates new CoordinateSystem by translation along \a vector
   */
  const CoordinateSystem& make_translation(const CoordinateSystem& cs,
                                           QuantityVector<length_d> const& vector);

  /**
   * creates a new CoordinateSystem in which vVec points in direction of the new z-axis
   */
  template <typename TDim>
  const CoordinateSystem& make_rotationToZ(const CoordinateSystem& cs,
                                           Vector<TDim> const& vVec);

  /**
   * creates a new CoordinateSystem, rotated around axis by angle.
   */
  template <typename TDim>
  const CoordinateSystem& make_rotation(const CoordinateSystem& cs,
                                        QuantityVector<TDim> const& axis,
                                        double const angle);

  /**
   * creates a new CoordinateSystem, translated by \a translation and rotated around \a
   * axis by \a angle.
   */
  template <typename TDim>
  const CoordinateSystem& make_translationAndRotation(
      const CoordinateSystem& cs, QuantityVector<length_d> const& translation,
      QuantityVector<TDim> const& axis, double const angle);

  /**
   * A class to store the reference coordinate system for a geometric object
   *
   * A CoordinateSystem is stored in the CoordinateSystemRegistry and can only be created
   * relative to other CoordinateSystems. The geometric transformation between all
   * CoordinateSystems is always known and cached.
   *
   * The static (singleton) function \ref get_root_CoordinateSystem is the way to access
   * the global root CoordinateSystem. CoordinateSystem objects are accessed via const
   * references from the registry.
   *
   * New CoordinateSystem are created using: \ref make_rotationToZ, \ref make_rotation,
   * or \ref make_translationAndRotation.
   *
   * CoordinateSystem objects are immutable after construction and stored in the registry.
   */

  class CoordinateSystem {

  public:
    /**
     * Constructor from parent CS and transformation matrix
     */
    CoordinateSystem(const CoordinateSystem* parent, EigenTransform const& transf);

    /**
     * Constructor for root CS
     */
    CoordinateSystem();

    CoordinateSystem(CoordinateSystem&&) noexcept = default;
    CoordinateSystem& operator=(CoordinateSystem&&) noexcept = default;

    // default resource allocation
    CoordinateSystem(CoordinateSystem const&) = delete;
    CoordinateSystem& operator=(CoordinateSystem const& pCS) =
        delete; // avoid making copies
    ~CoordinateSystem() = default;

    /**
     * Checks if this is the root CS
     */
    bool isRoot() const { return !parent_; }

    /**
     * Get parent coordinate system
     */
    const CoordinateSystem* getParent() const;

    /**
     * Get transformation from this CS to parent
     */
    EigenTransform const& getTransformToParent() const;

    /**
     * Get cached transformation from this CS to root (lazy computed)
     */
    EigenTransform const& getTransformToRoot() const;

    /**
     * Get cached transformation from root to this CS (lazy computed)
     */
    EigenTransform const& getTransformFromRoot() const;

    bool operator==(CoordinateSystem const&) const;
    bool operator!=(CoordinateSystem const&) const;

  protected:
    /**
     * \name Friends
     * Manipulation and creation functions.
     * \{
     **/

    friend class CoordinateSystemRegistry;

    friend const CoordinateSystem& get_root_CoordinateSystem();

    friend const CoordinateSystem& make_translation(
        const CoordinateSystem& cs, QuantityVector<length_d> const& vector);
    template <typename TDim>
    friend const CoordinateSystem& make_rotationToZ(const CoordinateSystem& cs,
                                                    Vector<TDim> const& vVec);
    template <typename TDim>
    friend const CoordinateSystem& make_rotation(const CoordinateSystem& cs,
                                                 QuantityVector<TDim> const& axis,
                                                 double const angle);
    template <typename TDim>
    friend const CoordinateSystem& make_translationAndRotation(
        const CoordinateSystem& cs, QuantityVector<length_d> const& translation,
        QuantityVector<TDim> const& axis, double const angle);

    /** \} **/

  private:
    const CoordinateSystem* parent_;
    EigenTransform transformToParent_;

    // Lazy-computed cached transforms (thread-safe via std::call_once)
    mutable EigenTransform transformToRoot_;
    mutable EigenTransform transformFromRoot_;
    mutable std::once_flag cacheFlag_;
  };

  /**
   * Transformation matrix from one reference system to another.
   *
   * returns the transformation matrix necessary to transform primitives with coordinates
   * in \a pFrom to \a pTo, e.g.
   * \f$ \vec{v}^{\text{(to)}} = \mathcal{M} \vec{v}^{\text{(from)}} \f$
   * (\f$ \vec{v}^{(.)} \f$ denotes the coordinates/components of the component in
   * the indicated CoordinateSystem).
   *
   * \todo make this a protected member of CoordinateSystem
   */
  EigenTransform get_transformation(CoordinateSystem const& c1,
                                    CoordinateSystem const& c2);

} // namespace corsika

#include <corsika/detail/framework/geometry/CoordinateSystem.inl>