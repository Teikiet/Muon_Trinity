/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <Eigen/Dense>

namespace corsika {

  /**
   * @defgroup Utilities
   *
   * Collection of classes and methods to perform recurring tasks.
   */

  /**
   * @class COMBoost
   * @ingroup Utilities
   *
   * This utility class handles Lorentz boost (in one spatial direction)
   * between different referenence frames, using FourVector.
   *
   * The class is initialized with projectile and optionally target
   * energy/momentum data. During initialization, a rotation matrix is
   * calculated to represent the projectile movement (and thus the
   * boost) along the z-axis. Also the inverse of this rotation is
   * calculated. The Lorentz boost matrix and its inverse are
   * determined as 2x2 matrices considering the energy and
   * pz-momentum.
   *
   * Different constructors are offered with different specialization
   * for the cases of collisions (projectile-target) or just decays
   * (projectile only).
   *
   * Two sets of constructors are provided:
   * 1. Primary constructors with explicit originalCS parameter (first argument)
   * 2. Convenience constructors assuming originalCS = rootCS
   *
   * Use explicit CS constructors when working with lab frames different from root:
   *   auto const& labCS = make_rotation(...);
   *   COMBoost boost(labCS, P4projectile, massTarget);
   *
   * Use convenience constructors for typical cases (originalCS = rootCS):
   *   COMBoost boost(P4projectile, massTarget);
   */

  class COMBoost {

  public:
    // ========================================================================
    // Constructors (assume root CS as originalCS)
    // ========================================================================

    /**
     * Construct a COMBoost given four-vector of projectile and mass of target (target at
     * rest).
     *
     * The FourMomentum and mass define the lab system.
     */
    COMBoost(FourMomentum const& P4projectile, HEPEnergyType const massTarget);

    /**
     * Construct a COMBoost to boost into the rest frame of a particle given its
     * 3-momentum and mass.
     */
    COMBoost(MomentumVector const& momentum, HEPEnergyType const mass);

    /**
     * Construct a COMBoost given two four-vectors of projectile target.
     *
     * The two FourMomentum can define an arbitrary system.
     */
    COMBoost(FourMomentum const& P4projectile, FourMomentum const& P4target);

    //! transforms a 4-momentum from lab frame to the center-of-mass frame
    template <typename FourVector>
    FourVector toCoM(FourVector const& p4) const;

    //! transforms a 4-momentum from the center-of-mass frame back to lab frame
    template <typename FourVector>
    FourVector fromCoM(FourVector const& p4) const;

    /**
     * Wrap raw generator components (in boosted frame, z = projectile direction)
     * into a MomentumVector in root CS.
     *
     * Use this instead of MomentumVector(getRotatedCS(), components) at all call sites.
     */
    MomentumVector makeCoMVector(
        QuantityVector<hepmomentum_d> const& rawComponents) const;

    //! returns the rotation matrix (rotated frame -> root CS)
    Eigen::Matrix3d getRotation() const;

    //! returns the inverse rotation matrix (root CS -> rotated frame)
    Eigen::Matrix3d getInverseRotation() const;

  protected:
    //! internal method
    void setBoost(double const coshEta, double const sinhEta);
    void setRotation(Eigen::Vector3d const& direction);

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Matrix2d boost_;
    Eigen::Matrix2d inverseBoost_;
    Eigen::Matrix3d rotation_; // rotated frame (z = projectile direction) -> root CS
    Eigen::Matrix3d inverseRotation_; // root CS -> rotated frame
  };
} // namespace corsika

#include <corsika/detail/framework/utility/COMBoost.inl>