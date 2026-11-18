/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <utility>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  template <typename TDimension>
  inline typename QuantityVector<TDimension>::quantity_type QuantityVector<TDimension>::
  operator[](size_t const index) const {
    return quantity_type(phys::units::detail::magnitude_tag, eigenVector_[index]);
  }

  template <typename TDimension>
  inline typename QuantityVector<TDimension>::quantity_type
  QuantityVector<TDimension>::getX() const {
    return (*this)[0];
  }

  template <typename TDimension>
  inline typename QuantityVector<TDimension>::quantity_type
  QuantityVector<TDimension>::getY() const {
    return (*this)[1];
  }

  template <typename TDimension>
  inline typename QuantityVector<TDimension>::quantity_type
  QuantityVector<TDimension>::getZ() const {
    return (*this)[2];
  }

  template <typename TDimension>
  inline typename QuantityVector<TDimension>::quantity_type
  QuantityVector<TDimension>::getNorm() const {
    return quantity_type(phys::units::detail::magnitude_tag, eigenVector_.norm());
  }

  template <typename TDimension>
  inline typename QuantityVector<TDimension>::quantity_square_type
  QuantityVector<TDimension>::getSquaredNorm() const {
    using QuantitySquared =
        decltype(std::declval<quantity_type>() * std::declval<quantity_type>());
    return QuantitySquared(phys::units::detail::magnitude_tag,
                           eigenVector_.squaredNorm());
  }

  template <typename TDimension>
  inline QuantityVector<TDimension> QuantityVector<TDimension>::operator+(
      QuantityVector<TDimension> const& pQVec) const {
    return QuantityVector<TDimension>(eigenVector_ + pQVec.eigenVector_);
  }

  template <typename TDimension>
  inline QuantityVector<TDimension> QuantityVector<TDimension>::operator-(
      QuantityVector<TDimension> const& pQVec) const {
    return QuantityVector<TDimension>(eigenVector_ - pQVec.eigenVector_);
  }

  template <typename TDimension>
  template <typename TScalarDim>
  inline auto QuantityVector<TDimension>::operator*(
      phys::units::quantity<TScalarDim, double> const p) const {
    using ResQuantity =
        phys::units::detail::Product<TScalarDim, TDimension, double, double>;

    if constexpr (std::is_same<ResQuantity, double>::value) // result dimensionless, not
                                                            // a "quantity_type" anymore
    {
      return QuantityVector<phys::units::dimensionless_d>(eigenVector_ * p.magnitude());
    } else {
      return QuantityVector<typename ResQuantity::dimension_type>(eigenVector_ *
                                                                  p.magnitude());
    }
  }

  template <typename TDimension>
  template <typename TScalarDim>
  inline auto QuantityVector<TDimension>::operator/(
      phys::units::quantity<TScalarDim, double> const p) const {
    return (*this) * (1 / p);
  }

  template <typename TDimension>
  inline auto QuantityVector<TDimension>::operator*(double const p) const {
    return QuantityVector<TDimension>(eigenVector_ * p);
  }

  template <typename TDimension>
  inline auto QuantityVector<TDimension>::operator/(double const p) const {
    return QuantityVector<TDimension>(eigenVector_ / p);
  }

  template <typename TDimension>
  inline auto& QuantityVector<TDimension>::operator/=(double const p) {
    eigenVector_ /= p;
    return *this;
  }

  template <typename TDimension>
  inline auto& QuantityVector<TDimension>::operator*=(double const p) {
    eigenVector_ *= p;
    return *this;
  }

  template <typename TDimension>
  inline auto& QuantityVector<TDimension>::operator+=(
      QuantityVector<TDimension> const& pQVec) {
    eigenVector_ += pQVec.eigenVector_;
    return *this;
  }

  template <typename TDimension>
  inline auto& QuantityVector<TDimension>::operator-=(
      QuantityVector<TDimension> const& pQVec) {
    eigenVector_ -= pQVec.eigenVector_;
    return *this;
  }

  template <typename TDimension>
  inline auto QuantityVector<TDimension>::operator-() const {
    return QuantityVector<TDimension>(-eigenVector_);
  }

  template <typename TDimension>
  inline auto QuantityVector<TDimension>::normalized() const {
    return QuantityVector<TDimension>(eigenVector_.normalized());
  }

  template <typename TDimension>
  inline auto QuantityVector<TDimension>::operator==(
      QuantityVector<TDimension> const& p) const {
    return eigenVector_ == p.eigenVector_;
  }

  template <typename TDimension>
  inline std::ostream& operator<<(std::ostream& os,
                                  corsika::QuantityVector<TDimension> const& qv) {
    using quantity_type = phys::units::quantity<TDimension, double>;

    os << '(' << qv.eigenVector_(0) << ' ' << qv.eigenVector_(1) << ' '
       << qv.eigenVector_(2) << ") "
       << phys::units::to_unit_symbol<TDimension, double>(
              quantity_type(phys::units::detail::magnitude_tag, 1));
    return os;
  }

} // namespace corsika
