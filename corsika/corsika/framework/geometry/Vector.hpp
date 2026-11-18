/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/BaseVector.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>

namespace corsika {

  /**
   * A Vector represents a 3-vector in Euclidean space.
   *
   * Vectors are stored internally with components in the root coordinate system.
   * This design reduces memory footprint (no CS pointer stored) and optimizes
   * the hot path (direct component access without transformation).
   *
   * Physical dimensions are part of the type, preventing mixing incompatible
   * quantities (e.g., electric and magnetic fields).
   *
   * Vectors transform via rotation only (translation-invariant).
   */

  template <typename TDimension>
  class Vector : public BaseVector<TDimension> {
  public:
    using quantity_type = phys::units::quantity<TDimension, double>;
    using quantity_square_type =
        decltype(std::declval<quantity_type>() * std::declval<quantity_type>());

    /**
     * Construct Vector with components already in root CS (direct storage)
     */
    Vector(QuantityVector<TDimension> const& pQVector)
        : BaseVector<TDimension>(pQVector) {}

    /**
     * Construct Vector with components in root CS from individual values
     */
    Vector(quantity_type const x = quantity_type::zero(),
           quantity_type const y = quantity_type::zero(),
           quantity_type const z = quantity_type::zero())
        : BaseVector<TDimension>(QuantityVector<TDimension>(x, y, z)) {}

    /**
     * Construct Vector from components in arbitrary CS (transforms to root)
     */
    Vector(const CoordinateSystem& cs, QuantityVector<TDimension> const& pQVector);

    /**
     * Construct Vector from individual components in arbitrary CS (transforms to root)
     */
    Vector(const CoordinateSystem& cs, quantity_type const x, quantity_type const y,
           quantity_type const z);

    /**
     * Get components in root CS (zero overhead, const reference)
     */
    QuantityVector<TDimension> const& getComponents() const;

    /**
     * Get mutable components in root CS (zero overhead, non-const reference)
     */
    QuantityVector<TDimension>& getComponents();

    /**
     * Get components transformed to arbitrary CS (computed on-demand)
     */
    QuantityVector<TDimension> getComponents(const CoordinateSystem& pCS) const;

    /**
     * Get X, Y, Z components in root CS
     */
    quantity_type getX() const;
    quantity_type getY() const;
    quantity_type getZ() const;

    /**
     * Get X, Y, Z components in arbitrary CS
     *
     *
     * Note: If you need multiple components (X, Y, Z),
     * call getComponents(cs) once and access components individually:
     *
     *   // Instead of:
     *   auto x = vec.getX(cs);
     *   auto y = vec.getY(cs);
     *   auto z = vec.getZ(cs);
     *
     *   // Do:
     *   auto comp = vec.getComponents(cs);
     *   auto x = comp.getX();
     *   auto y = comp.getY();
     *   auto z = comp.getZ();
     *
     */
    quantity_type getX(const CoordinateSystem& pCS) const;
    quantity_type getY(const CoordinateSystem& pCS) const;
    quantity_type getZ(const CoordinateSystem& pCS) const;

    /*!
     * returns the norm/length of the Vector. Before using this method,
     * think about whether squaredNorm() might be cheaper for your computation.
     */
    quantity_type getNorm() const;

    /*!
     * returns the squared norm of the Vector. Before using this method,
     * think about whether norm() might be cheaper for your computation.
     */
    quantity_square_type getSquaredNorm() const;

    /**
     * Returns parallel projection of this vector onto another vector.
     */
    template <typename TDimension2>
    auto getParallelProjectionOnto(Vector<TDimension2> const& pVec) const;

    Vector operator+(Vector<TDimension> const& pVec) const;

    Vector operator-(Vector<TDimension> const& pVec) const;

    auto& operator*=(double const p);

    template <typename TScalarDim>
    auto operator*(phys::units::quantity<TScalarDim, double> const p) const;
    template <typename TScalarDim>
    auto operator/(phys::units::quantity<TScalarDim, double> const p) const;

    auto operator*(double const p) const;

    auto operator/(double const p) const;

    auto& operator+=(Vector<TDimension> const& pVec);

    auto& operator-=(Vector<TDimension> const& pVec);

    auto operator-() const;

    auto normalized() const;

    template <typename TDimension2>
    auto cross(Vector<TDimension2> const& pV) const;

    template <typename TDimension2>
    auto dot(Vector<TDimension2> const& pV) const;
  };

  /**
   * Free operator to allow commutative multiplications of quantities and Vector.
   *
   * @tparam TDimension
   * @tparam UDimension
   * @param n
   * @param vec
   * @return auto
   */
  template <typename TDimension, typename UDimension>
  Vector<phys::units::detail::product_d<TDimension, UDimension>> operator*(
      quantity<UDimension> const n, Vector<TDimension> const& vec);

  /**
   * Free operator to allow commutative multiplications of normal double with Vector.
   *
   * @tparam TDimension
   * @tparam UDimension
   * @param n
   * @param vec
   * @return auto
   */
  template <typename TDimension>
  Vector<TDimension> operator*(double const n, Vector<TDimension> const& vec);

} // namespace corsika

#include <corsika/detail/framework/geometry/Vector.inl>