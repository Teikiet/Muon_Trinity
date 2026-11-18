/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <Eigen/Dense>

#include <stdexcept>
#include <mutex>

namespace corsika {

  inline const CoordinateSystem* CoordinateSystem::getParent() const { return parent_; }

  inline EigenTransform const& CoordinateSystem::getTransformToParent() const {
    return transformToParent_;
  }

  inline EigenTransform const& CoordinateSystem::getTransformToRoot() const {
    std::call_once(cacheFlag_, [this]() {
      if (parent_ == nullptr) {
        // Root CS: identity
        transformToRoot_ = EigenTransform::Identity();
        transformFromRoot_ = EigenTransform::Identity();
      } else {
        // Compose transformations up to root
        transformToRoot_ = parent_->getTransformToRoot() * transformToParent_;
        transformFromRoot_ = transformToRoot_.inverse(Eigen::TransformTraits::Isometry);
      }
    });
    return transformToRoot_;
  }

  inline EigenTransform const& CoordinateSystem::getTransformFromRoot() const {
    // Trigger cache computation if needed
    getTransformToRoot();
    return transformFromRoot_;
  }

  inline bool CoordinateSystem::operator==(CoordinateSystem const& cs) const {
    if (this == &cs) return true;
    return parent_ == cs.parent_ &&
           transformToParent_.matrix() == cs.transformToParent_.matrix();
  }

  inline bool CoordinateSystem::operator!=(CoordinateSystem const& cs) const {
    return !(cs == *this);
  }

  /**
   * Transformation matrix from one coordinate system to another.
   * Uses cached transformations to/from root for efficiency.
   */
  inline EigenTransform get_transformation(CoordinateSystem const& pFrom,
                                           CoordinateSystem const& pTo) {
    // Transform: from -> root -> to
    // T = T_to_from_root * T_from_to_root
    return pTo.getTransformFromRoot() * pFrom.getTransformToRoot();
  }

  inline const CoordinateSystem& make_translation(
      const CoordinateSystem& cs, QuantityVector<length_d> const& vector) {
    EigenTransform const translation{EigenTranslation(vector.getEigenVector())};
    return CoordinateSystemRegistry::getInstance().registerCS(&cs, translation);
  }

  template <typename TDim>
  inline const CoordinateSystem& make_rotationToZ(const CoordinateSystem& cs,
                                                  Vector<TDim> const& vVec) {
    auto const vVecComp = vVec.getComponents(cs);
    if (vVecComp.getX().magnitude() == 0 && vVecComp.getY().magnitude() == 0 &&
        vVecComp.getZ().magnitude() == 0) {
      return cs;
    }

    auto const a = vVecComp.normalized().getEigenVector();
    auto const a1 = a(0), a2 = a(1), a3 = a(2);

    Eigen::Matrix3d A, B;

    if (a3 > 0) {
      auto const c = 1 / (1 + a3);
      A << 1, 0, a1,                      // comment to prevent clang-format
          0, 1, a2,                       // .
          -a1, -a2, 1;                    // .
      B << -a1 * a1 * c, -a1 * a2 * c, 0, // .
          -a1 * a2 * c, -a2 * a2 * c, 0,  // .
          0, 0, -(a1 * a1 + a2 * a2) * c; // .

    } else {
      auto const c = 1 / (1 - a3);
      A << 1, 0, a1,                      // .
          0, -1, a2,                      // .
          a1, -a2, -1;                    // .
      B << -a1 * a1 * c, +a1 * a2 * c, 0, // .
          -a1 * a2 * c, +a2 * a2 * c, 0,  // .
          0, 0, (a1 * a1 + a2 * a2) * c;  // .
    }

    return CoordinateSystemRegistry::getInstance().registerCS(&cs, EigenTransform{A + B});
  }

  template <typename TDim>
  inline const CoordinateSystem& make_rotation(const CoordinateSystem& cs,
                                               QuantityVector<TDim> const& axis,
                                               double const angle) {
    if (axis.getEigenVector().isZero()) {
      throw std::runtime_error("null-vector given as axis parameter");
    }

    EigenTransform const rotation{
        Eigen::AngleAxisd(angle, axis.getEigenVector().normalized())};

    return CoordinateSystemRegistry::getInstance().registerCS(&cs, rotation);
  }

  template <typename TDim>
  inline const CoordinateSystem& make_translationAndRotation(
      const CoordinateSystem& cs, QuantityVector<length_d> const& translation,
      QuantityVector<TDim> const& axis, double const angle) {
    if (axis.getEigenVector().isZero()) {
      throw std::runtime_error("null-vector given as axis parameter");
    }

    EigenTransform const transf{
        Eigen::AngleAxisd(angle, axis.getEigenVector().normalized()) *
        EigenTranslation(translation.getEigenVector())};

    return CoordinateSystemRegistry::getInstance().registerCS(&cs, transf);
  }

} // namespace corsika
