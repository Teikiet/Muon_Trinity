/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <cmath>

#include <Eigen/Dense>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/Logging.hpp>

namespace corsika {

  inline void COMBoost::setRotation(Eigen::Vector3d const& direction) {
    double const norm = direction.norm();
    if (norm == 0.) {
      rotation_ = Eigen::Matrix3d::Identity();
      inverseRotation_ = Eigen::Matrix3d::Identity();
      return;
    }

    auto const a = direction / norm;
    double const a1 = a(0), a2 = a(1), a3 = a(2);

    Eigen::Matrix3d A, B;
    if (a3 > 0) {
      double const c = 1. / (1. + a3);
      A << 1, 0, a1, 0, 1, a2, -a1, -a2, 1;
      B << -a1 * a1 * c, -a1 * a2 * c, 0, -a1 * a2 * c, -a2 * a2 * c, 0, 0, 0,
          -(a1 * a1 + a2 * a2) * c;
    } else {
      double const c = 1. / (1. - a3);
      A << 1, 0, a1, 0, -1, a2, a1, -a2, -1;
      B << -a1 * a1 * c, +a1 * a2 * c, 0, -a1 * a2 * c, +a2 * a2 * c, 0, 0, 0,
          (a1 * a1 + a2 * a2) * c;
    }

    rotation_ = A + B;
    inverseRotation_ = rotation_.transpose();
  }

  // ==========================================================================
  // Primary constructors assuming originalCS = rootCS
  // ==========================================================================

  inline COMBoost::COMBoost(FourMomentum const& P4projectile,
                            HEPMassType const massTarget) {
    auto const& pProjectile = P4projectile.getSpaceLikeComponents();
    auto const pProjNormSquared = pProjectile.getSquaredNorm();
    auto const pProjNorm = sqrt(pProjNormSquared);

    setRotation(pProjectile.getComponents().getEigenVector());

    auto const eProjectile = P4projectile.getTimeLikeComponent();
    auto const massProjectileSquared = eProjectile * eProjectile - pProjNormSquared;
    auto const s =
        massTarget * massTarget + massProjectileSquared + 2 * eProjectile * massTarget;

    auto const sqrtS = sqrt(s);
    auto const sinhEta = -pProjNorm / sqrtS;
    auto const coshEta = sqrt(1 + pProjNormSquared / s);

    setBoost(coshEta, sinhEta);
    CORSIKA_LOG_TRACE("COMBoost (1-beta)={}, gamma={}, det={}", 1 - sinhEta / coshEta,
                      coshEta, boost_.determinant() - 1);
  }

  inline COMBoost::COMBoost(FourMomentum const& P4projectile,
                            FourMomentum const& P4target) {
    // this is the center-of-momentum CM frame
    auto const pCM =
        P4projectile.getSpaceLikeComponents() + P4target.getSpaceLikeComponents();
    auto const pCM2 = pCM.getSquaredNorm();
    auto const pCMnorm = sqrt(pCM2);

    setRotation(pCM.getComponents().getEigenVector());

    auto const s = (P4projectile + P4target).getNormSqr();
    auto const sqrtS = sqrt(s);
    auto const sinhEta = -pCMnorm / sqrtS;
    auto const coshEta = sqrt(1 + pCM2 / s);

    setBoost(coshEta, sinhEta);
    CORSIKA_LOG_TRACE("COMBoost (1-beta)={}, gamma={}, det={}", 1 - sinhEta / coshEta,
                      coshEta, boost_.determinant() - 1);
  }

  inline COMBoost::COMBoost(MomentumVector const& momentum, HEPEnergyType const mass) {
    auto const squaredNorm = momentum.getSquaredNorm();
    auto const norm = sqrt(squaredNorm);

    setRotation(momentum.getComponents().getEigenVector());

    auto const sinhEta = -norm / mass;
    auto const coshEta = sqrt(1 + squaredNorm / (mass * mass));
    setBoost(coshEta, sinhEta);
    CORSIKA_LOG_TRACE("COMBoost (1-beta)={}, gamma={}, det={}", 1 - sinhEta / coshEta,
                      coshEta, boost_.determinant() - 1);
  }

  // ==========================================================================
  // Transformation methods
  // ==========================================================================

  template <typename FourVector>
  inline FourVector COMBoost::toCoM(FourVector const& p4) const {
    // rotate momentum from root CS to rotated frame (z = projectile direction)
    Eigen::Vector3d eVecRotated =
        inverseRotation_ * p4.getSpaceLikeComponents().getComponents().getEigenVector();
    Eigen::Vector2d lab;

    lab << (p4.getTimeLikeComponent() * (1 / 1_GeV)),
        (eVecRotated(2) * (1 / 1_GeV).magnitude());

    auto const boostedZ = boost_ * lab;
    auto const E_CoM = boostedZ(0) * 1_GeV;

    eVecRotated(2) = boostedZ(1) * (1_GeV).magnitude();

    // result stays in rotated frame
    CORSIKA_LOG_TRACE("E0={}, p={}, E0'={}, p'={}", p4.getTimeLikeComponent() / 1_GeV,
                      eVecRotated(2) * (1 / 1_GeV).magnitude(), E_CoM / 1_GeV, boostedZ);

    return FourVector(E_CoM, MomentumVector(QuantityVector<hepmomentum_d>(eVecRotated)));
  }

  template <typename FourVector>
  inline FourVector COMBoost::fromCoM(FourVector const& p4) const {
    // p4 components are in rotated frame (constructed via makeCoMVector or toCoM)
    Eigen::Vector3d pVec = p4.getSpaceLikeComponents().getComponents().getEigenVector();

    auto const Ecm = p4.getTimeLikeComponent();

    Eigen::Vector2d com;
    com << (Ecm * (1 / 1_GeV)), (pVec(2) * (1 / 1_GeV).magnitude());

    CORSIKA_LOG_TRACE("Ecm={} GeV, pcm={} GeV (norm = {} GeV), invariant mass={} GeV",
                      Ecm / 1_GeV, p4.getSpaceLikeComponents().getComponents() / 1_GeV,
                      p4.getSpaceLikeComponents().getNorm() / 1_GeV,
                      p4.getNorm() / 1_GeV);

    auto const boostedZ = inverseBoost_ * com;
    auto const E_lab = boostedZ(0) * 1_GeV;

    pVec(2) = boostedZ(1) * (1_GeV).magnitude();

    // rotate from rotated frame back to root CS
    MomentumVector const pLab(QuantityVector<hepmomentum_d>(rotation_ * pVec));

    CORSIKA_LOG_TRACE(
        "COMBoost::fromCoM --> Elab={} GeV,"
        " plab={} GeV (norm={} GeV),"
        " invariant mass={}",
        E_lab / 1_GeV, pLab.getComponents() / 1_GeV, pLab.getNorm() / 1_GeV,
        FourVector{E_lab, pLab}.getNorm() / 1_GeV);

    return FourVector{E_lab, pLab};
  }

  inline MomentumVector COMBoost::makeCoMVector(
      QuantityVector<hepmomentum_d> const& rawComponents) const {
    return MomentumVector(rawComponents);
  }

  inline Eigen::Matrix3d COMBoost::getRotation() const { return rotation_; }

  inline Eigen::Matrix3d COMBoost::getInverseRotation() const { return inverseRotation_; }

  inline void COMBoost::setBoost(double const coshEta, double const sinhEta) {
    boost_ << coshEta, sinhEta, sinhEta, coshEta;
    inverseBoost_ << coshEta, -sinhEta, -sinhEta, coshEta;
  }

} // namespace corsika