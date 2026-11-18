/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <Eigen/Dense>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <ostream>
#include <utility>

namespace corsika {

  class CoordinateSystem; // fwd decl
  class Point;            // fwd decl
  template <typename T>
  class Vector; // fwd decl

  /*!
   * A QuantityVector is a three-component container based on Eigen::Vector3d
   * with a phys::units::si::dimension. Arithmethic operators are defined that
   * propagate the dimensions by dimensional analysis.
   *
   * \todo review QuantityVector: can this be a protected-only
   * inheritance to other objects? We don't want to expose access to
   * low-level Eigen objects anywhere.
   */

  template <typename TDimension>
  class QuantityVector {
  public:
    using quantity_type =
        phys::units::quantity<TDimension, double>; //< the phys::units::quantity
                                                   // corresponding to the dimension
    using quantity_square_type =
        decltype(std::declval<quantity_type>() * std::declval<quantity_type>());

    QuantityVector(Eigen::Vector3d const& pBareVector)
        : eigenVector_(pBareVector) {}

  public:
    typedef TDimension dimension_type; //!< should be a phys::units::dimension

    QuantityVector(quantity_type const a, quantity_type const b, quantity_type const c)
        : eigenVector_{a.magnitude(), b.magnitude(), c.magnitude()} {}

    QuantityVector(double const a, double const b, double const c)
        : eigenVector_{a, b, c} {
      static_assert(
          std::is_same_v<TDimension, phys::units::dimensionless_d>,
          "initialization of dimensionful QuantityVector with pure numbers not allowed!");
    }

    quantity_type operator[](size_t const index) const;
    quantity_type getX() const;
    quantity_type getY() const;
    quantity_type getZ() const;
    Eigen::Vector3d const& getEigenVector() const { return eigenVector_; }
    Eigen::Vector3d& getEigenVector() { return eigenVector_; }

    quantity_type getNorm() const;

    quantity_square_type getSquaredNorm() const;

    QuantityVector operator+(QuantityVector<TDimension> const& pQVec) const;

    QuantityVector operator-(QuantityVector<TDimension> const& pQVec) const;

    template <typename TScalarDim>
    auto operator*(phys::units::quantity<TScalarDim, double> const p) const;

    template <typename TScalarDim>
    auto operator/(phys::units::quantity<TScalarDim, double> const p) const;

    auto operator*(double const p) const;

    auto operator/(double const p) const;

    auto& operator/=(double const p);

    auto& operator*=(double const p);

    auto& operator+=(QuantityVector<TDimension> const& pQVec);

    auto& operator-=(QuantityVector<TDimension> const& pQVec);

    auto operator-() const;

    auto normalized() const;

    auto operator==(QuantityVector<TDimension> const& p) const;

    // friends:
    friend class CoordinateSystem;
    friend class Point;
    template <typename T>
    friend class corsika::Vector;
    template <typename TDim>
    friend std::ostream& operator<<(std::ostream& os, QuantityVector<TDim> const& qv);

  protected:
    Eigen::Vector3d
        eigenVector_; //!< the actual container where the raw numbers are stored
  };

  /**
   * streaming operator
   **/

  template <typename TDimension>
  std::ostream& operator<<(std::ostream& os,
                           corsika::QuantityVector<TDimension> const& qv);

} // namespace corsika

#include <corsika/detail/framework/geometry/QuantityVector.inl>
