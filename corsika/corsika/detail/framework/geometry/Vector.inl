/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/BaseVector.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  // Constructors that transform to root CS
  template <typename TDimension>
  inline Vector<TDimension>::Vector(const CoordinateSystem& cs,
                                    QuantityVector<TDimension> const& pQVector)
      : BaseVector<TDimension>(
            cs.isRoot() ? pQVector
                        : QuantityVector<TDimension>(cs.getTransformToRoot().linear() *
                                                     pQVector.eigenVector_)) {}

  template <typename TDimension>
  inline Vector<TDimension>::Vector(const CoordinateSystem& cs,
                                    typename Vector<TDimension>::quantity_type const x,
                                    typename Vector<TDimension>::quantity_type const y,
                                    typename Vector<TDimension>::quantity_type const z)
      : Vector(cs, QuantityVector<TDimension>(x, y, z)) {}

  // getComponents: direct access in root CS (zero overhead)
  template <typename TDimension>
  inline QuantityVector<TDimension> const& Vector<TDimension>::getComponents() const {
    return BaseVector<TDimension>::getQuantityVector();
  }

  template <typename TDimension>
  inline QuantityVector<TDimension>& Vector<TDimension>::getComponents() {
    return BaseVector<TDimension>::getQuantityVector();
  }

  // getComponents(cs): transform from root to target CS
  template <typename TDimension>
  inline QuantityVector<TDimension> Vector<TDimension>::getComponents(
      const CoordinateSystem& pCS) const {
    if (pCS.isRoot()) { return BaseVector<TDimension>::getQuantityVector(); }
    // Transform: root -> pCS (rotation only for vectors)
    return QuantityVector<TDimension>(
        pCS.getTransformFromRoot().linear() *
        BaseVector<TDimension>::getQuantityVector().eigenVector_);
  }

  // getX/Y/Z: direct access to root CS components
  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getX() const {
    return getComponents()[0];
  }

  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getY() const {
    return getComponents()[1];
  }

  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getZ() const {
    return getComponents()[2];
  }

  // getX/Y/Z: transform individual components to target CS
  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getX(
      const CoordinateSystem& pCS) const {
    return getComponents(pCS)[0];
  }

  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getY(
      const CoordinateSystem& pCS) const {
    return getComponents(pCS)[1];
  }

  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getZ(
      const CoordinateSystem& pCS) const {
    return getComponents(pCS)[2];
  }

  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_type Vector<TDimension>::getNorm() const {
    return BaseVector<TDimension>::getQuantityVector().getNorm();
  }

  template <typename TDimension>
  inline typename Vector<TDimension>::quantity_square_type
  Vector<TDimension>::getSquaredNorm() const {
    return BaseVector<TDimension>::getQuantityVector().getSquaredNorm();
  }

  template <typename TDimension>
  template <typename TDimension2>
  inline auto Vector<TDimension>::getParallelProjectionOnto(
      Vector<TDimension2> const& pVec) const {
    // Both vectors in root CS, no transformation needed
    auto const& a = getComponents().eigenVector_;
    auto const& b = pVec.getComponents().eigenVector_;

    return Vector<TDimension>(
        QuantityVector<TDimension>(b * ((a.dot(b)) / b.squaredNorm())));
  }

  template <typename TDimension>
  inline Vector<TDimension> Vector<TDimension>::operator+(
      Vector<TDimension> const& pVec) const {
    // Both in root CS, direct addition
    return Vector<TDimension>(getComponents() + pVec.getComponents());
  }

  template <typename TDimension>
  inline Vector<TDimension> Vector<TDimension>::operator-(
      Vector<TDimension> const& pVec) const {
    // Both in root CS, direct subtraction
    return Vector<TDimension>(getComponents() - pVec.getComponents());
  }

  template <typename TDimension>
  inline auto& Vector<TDimension>::operator*=(double const p) {
    BaseVector<TDimension>::getQuantityVector() *= p;
    return *this;
  }

  template <typename TDimension>
  template <typename TScalarDim>
  inline auto Vector<TDimension>::operator*(
      phys::units::quantity<TScalarDim, double> const p) const {
    using ProdDim = phys::units::detail::product_d<TDimension, TScalarDim>;

    return Vector<ProdDim>(BaseVector<TDimension>::getQuantityVector() * p);
  }

  template <typename TDimension>
  template <typename TScalarDim>
  inline auto Vector<TDimension>::operator/(
      phys::units::quantity<TScalarDim, double> const p) const {
    return (*this) * (1 / p);
  }

  template <typename TDimension>
  inline auto Vector<TDimension>::operator*(double const p) const {
    return Vector<TDimension>(BaseVector<TDimension>::getQuantityVector() * p);
  }

  template <typename TDimension>
  inline auto Vector<TDimension>::operator/(double const p) const {
    return Vector<TDimension>(BaseVector<TDimension>::getQuantityVector() / p);
  }

  template <typename TDimension>
  inline auto& Vector<TDimension>::operator+=(Vector<TDimension> const& pVec) {
    BaseVector<TDimension>::getQuantityVector() += pVec.getComponents();
    return *this;
  }

  template <typename TDimension>
  inline auto& Vector<TDimension>::operator-=(Vector<TDimension> const& pVec) {
    BaseVector<TDimension>::getQuantityVector() -= pVec.getComponents();
    return *this;
  }

  template <typename TDimension>
  inline auto Vector<TDimension>::operator-() const {
    return Vector<TDimension>(-BaseVector<TDimension>::getQuantityVector());
  }

  template <typename TDimension>
  inline auto Vector<TDimension>::normalized() const {
    return (*this) * (1 / getNorm());
  }

  template <typename TDimension>
  template <typename TDimension2>
  inline auto Vector<TDimension>::cross(Vector<TDimension2> const& pV) const {
    // Both in root CS, direct cross product
    auto const c1 = getComponents().eigenVector_;
    auto const c2 = pV.getComponents().eigenVector_;
    auto const bareResult = c1.cross(c2);

    using ProdDim = phys::units::detail::product_d<TDimension, TDimension2>;
    return Vector<ProdDim>(bareResult);
  }

  template <typename TDimension>
  template <typename TDimension2>
  inline auto Vector<TDimension>::dot(Vector<TDimension2> const& pV) const {
    // Both in root CS, direct dot product
    auto const c1 = getComponents().eigenVector_;
    auto const c2 = pV.getComponents().eigenVector_;
    auto const bareResult = c1.dot(c2);

    using ProdDim = phys::units::detail::product_d<TDimension, TDimension2>;

    return phys::units::quantity<ProdDim, double>(phys::units::detail::magnitude_tag,
                                                  bareResult);
  }

  template <typename TDimension>
  inline std::ostream& operator<<(std::ostream& os,
                                  corsika::Vector<TDimension> const& v) {
    auto const& qv = v.getComponents();
    os << qv << " (root CS)";
    return os;
  }

  /*
   * scalar * vector multiplication
   */
  template <typename TDimension, typename UDimension>
  inline Vector<phys::units::detail::product_d<TDimension, UDimension>> operator*(
      quantity<UDimension> const n, Vector<TDimension> const& vec) {
    return vec * n;
  }

  template <typename TDimension>
  inline Vector<TDimension> operator*(double const n, Vector<TDimension> const& vec) {
    return vec * n;
  }
} // namespace corsika