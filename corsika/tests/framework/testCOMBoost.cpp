/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/utility/COMBoost.hpp>

using namespace corsika;
using Catch::Approx;

double constexpr absMargin = 1e-6;

/**
 * \todo such helper functions should be moved to the FourVector class:
 */
// helper function for energy-momentum
// relativistic energy
auto const energy = [](HEPMassType const m, MomentumVector const& p) {
  return sqrt(m * m + p.getSquaredNorm());
};

auto const momentum = [](HEPEnergyType const E, HEPMassType const m) {
  return sqrt(E * E - m * m);
};

// helper function for mandelstam-s
auto const s = [](HEPEnergyType const E, QuantityVector<hepmomentum_d> const& p) {
  return E * E - p.getSquaredNorm();
};

TEST_CASE("rotation") {

  logging::set_level(logging::level::info);

  // define projectile kinematics in lab frame
  HEPMassType const projectileMass = 1_GeV;
  HEPMassType const targetMass = 1.0e300_eV;
  MomentumVector pProjectileLab{0_GeV, 0_PeV, 1_GeV};
  HEPEnergyType const eProjectileLab = energy(projectileMass, pProjectileLab);
  FourVector const PprojLab(eProjectileLab, pProjectileLab);

  MomentumVector e1(1_GeV, 0_GeV, 0_GeV);
  MomentumVector e2(0_GeV, 1_GeV, 0_GeV);
  MomentumVector e3(0_GeV, 0_GeV, 1_GeV);

  auto const toRotatedFrame = [](COMBoost const& boost, MomentumVector const& v) {
    Eigen::Vector3d const vRot =
        boost.getInverseRotation() * v.getComponents().getEigenVector();
    return MomentumVector(QuantityVector<hepmomentum_d>(vRot));
  };

  // define boost to com frame
  SECTION("pos. z-axis") {
    COMBoost boost({eProjectileLab, {0_GeV, 0_GeV, 1_GeV}}, targetMass);

    MomentumVector e1rot = toRotatedFrame(boost, e1);
    MomentumVector e2rot = toRotatedFrame(boost, e2);
    MomentumVector e3rot = toRotatedFrame(boost, e3);

    // length of e1, e2 and e3 must all be 1_GeV in rotated frame (not boosted!)
    CHECK(e1rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e2rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e3rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));

    // z-axis is along z-boost
    CHECK(e3rot.getX() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e3rot.getY() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e3rot.getZ() / 1_GeV == Approx(1).margin(absMargin));
  }

  SECTION("y-axis in upper half") {
    COMBoost boost({eProjectileLab, {0_GeV, 1_GeV, 1_meV}}, targetMass);

    MomentumVector e1rot = toRotatedFrame(boost, e1);
    MomentumVector e2rot = toRotatedFrame(boost, e2);
    MomentumVector e3rot = toRotatedFrame(boost, e3);

    // length of e1, e2 and e3 must all be 1_GeV in rotated CS (not boosted!)
    CHECK(e1rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e2rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e3rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));

    // z-axis is along y-boost
    CHECK(e2rot.getX() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e2rot.getY() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e2rot.getZ() / 1_GeV == Approx(1).margin(absMargin));
  }

  SECTION("x-axis in upper half") {
    COMBoost boost({eProjectileLab, {1_GeV, 0_GeV, 1_meV}}, targetMass);

    MomentumVector e1rot = toRotatedFrame(boost, e1);
    MomentumVector e2rot = toRotatedFrame(boost, e2);
    MomentumVector e3rot = toRotatedFrame(boost, e3);

    // length of e1, e2 and e3 must all be 1_GeV in rotated CS (not boosted!)
    CHECK(e1rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e2rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e3rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));

    // z-axis is along x-boost
    CHECK(e1rot.getX() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e1rot.getY() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e1rot.getZ() / 1_GeV == Approx(1).margin(absMargin));
  }

  SECTION("neg. z-axis") {
    COMBoost boost({eProjectileLab, {0_GeV, 0_GeV, -1_GeV}}, targetMass);

    MomentumVector e1rot = toRotatedFrame(boost, e1);
    MomentumVector e2rot = toRotatedFrame(boost, e2);
    MomentumVector e3rot = toRotatedFrame(boost, e3);

    // length of e1, e2 and e3 must all be 1_GeV in rotated CS (not boosted!)
    CHECK(e1rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e2rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e3rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));

    // z-axis is along -z-boost
    CHECK(-e3rot.getX() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(-e3rot.getY() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(-e3rot.getZ() / 1_GeV == Approx(1).margin(absMargin));
  }

  SECTION("x-axis lower half") {
    COMBoost boost({eProjectileLab, {1_GeV, 0_GeV, -1_meV}}, targetMass);

    MomentumVector e1rot = toRotatedFrame(boost, e1);
    MomentumVector e2rot = toRotatedFrame(boost, e2);
    MomentumVector e3rot = toRotatedFrame(boost, e3);

    // length of e1, e2 and e3 must all be 1_GeV in rotated CS (not boosted!)
    CHECK(e1rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e2rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e3rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));

    // z-axis is along x-boost
    CHECK(e1rot.getX() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e1rot.getY() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e1rot.getZ() / 1_GeV == Approx(1).margin(absMargin));
  }

  SECTION("y-axis lower half") {
    COMBoost boost({eProjectileLab, {0_GeV, 1_GeV, -1_meV}}, targetMass);

    MomentumVector e1rot = toRotatedFrame(boost, e1);
    MomentumVector e2rot = toRotatedFrame(boost, e2);
    MomentumVector e3rot = toRotatedFrame(boost, e3);

    // length of e1, e2 and e3 must all be 1_GeV in rotated CS (not boosted!)
    CHECK(e1rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e2rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));
    CHECK(e3rot.getNorm() / 1_GeV == Approx(1).margin(absMargin));

    // z-axis is along y-boost
    CHECK(e2rot.getX() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e2rot.getY() / 1_GeV == Approx(0).margin(absMargin));
    CHECK(e2rot.getZ() / 1_GeV == Approx(1).margin(absMargin));
  }
}

TEST_CASE("boosts") {

  logging::set_level(logging::level::trace);

  // define target kinematics in lab frame
  HEPMassType const targetMass = 1_GeV;
  MomentumVector pTargetLab{0_eV, 0_eV, 0_eV};
  HEPEnergyType const eTargetLab = energy(targetMass, pTargetLab);

  /*
    General tests check the interface and basic operation
   */

  SECTION("General tests") {

    // define projectile kinematics in lab frame
    HEPMassType const projectileMass = 1_GeV;
    MomentumVector pProjectileLab{0_GeV, 1_PeV, 0_GeV};
    HEPEnergyType const eProjectileLab = energy(projectileMass, pProjectileLab);
    FourVector const PprojLab(eProjectileLab, pProjectileLab);

    // define boost to com frame
    COMBoost boost(PprojLab, targetMass);

    // boost projecticle
    auto const PprojCoM = boost.toCoM(PprojLab);

    // boost target
    auto const PtargCoM = boost.toCoM(FourVector(targetMass, pTargetLab));

    // sum of momenta in CoM, should be 0
    auto const sumPCoM =
        PprojCoM.getSpaceLikeComponents() + PtargCoM.getSpaceLikeComponents();
    CHECK(sumPCoM.getNorm() / 1_GeV == Approx(0).margin(absMargin));

    // mandelstam-s should be invariant under transformation
    CHECK(s(eProjectileLab + eTargetLab,
            pProjectileLab.getComponents() + pTargetLab.getComponents()) /
              1_GeV / 1_GeV ==
          Approx(s(PprojCoM.getTimeLikeComponent() + PtargCoM.getTimeLikeComponent(),
                   PprojCoM.getSpaceLikeComponents().getComponents() +
                       PtargCoM.getSpaceLikeComponents().getComponents()) /
                 1_GeV / 1_GeV));

    // boost back...
    auto const PprojBack = boost.fromCoM(PprojCoM);

    // ...should yield original values before the boosts
    CHECK(PprojBack.getTimeLikeComponent() / PprojLab.getTimeLikeComponent() ==
          Approx(1));
    CHECK((PprojBack.getSpaceLikeComponents() - PprojLab.getSpaceLikeComponents())
                  .getNorm() /
              PprojLab.getSpaceLikeComponents().getNorm() ==
          Approx(0).margin(absMargin));
  }

  /*
    special case: projectile along -z
   */

  SECTION("Test boost along z-axis") {

    // define projectile kinematics in lab frame
    HEPMassType const projectileMass = 1_GeV;
    MomentumVector pProjectileLab{0_GeV, 0_PeV, -1_PeV};
    HEPEnergyType const eProjectileLab = energy(projectileMass, pProjectileLab);
    FourVector const PprojLab(eProjectileLab, pProjectileLab);

    auto const sqrt_s_lab =
        sqrt(s(eProjectileLab + targetMass, pProjectileLab.getComponents()));

    // define boost to com frame
    COMBoost boost(PprojLab, targetMass);

    // boost projecticle
    auto const PprojCoM = boost.toCoM(PprojLab);

    // toCoM() result is stored in raw rotated-frame components
    auto const a = PprojCoM.getSpaceLikeComponents().getComponents();
    CHECK(a.getX() / 1_GeV == Approx(0));
    CHECK(a.getY() / 1_GeV == Approx(0));
    CHECK(a.getZ() / (momentum(sqrt_s_lab / 2, projectileMass)) == Approx(1));

    // boost target
    auto const PtargCoM = boost.toCoM(FourVector(targetMass, pTargetLab));
    CHECK(PtargCoM.getTimeLikeComponent() / sqrt_s_lab == Approx(.5));

    // sum of momenta in CoM, should be 0
    auto const sumPCoM =
        PprojCoM.getSpaceLikeComponents() + PtargCoM.getSpaceLikeComponents();
    CHECK(sumPCoM.getNorm() / 1_GeV == Approx(0).margin(absMargin));
  }

  /*
    special case: projectile with arbitrary direction
   */

  SECTION("Test boost along tilted axis") {

    HEPMomentumType const P0 = 1_PeV;
    double theta = 33.;
    double phi = -10.;
    auto momentumComponents = [](double theta, double phi, HEPMomentumType ptot) {
      return std::make_tuple(ptot * sin(theta) * cos(phi), ptot * sin(theta) * sin(phi),
                             -ptot * cos(theta));
    };
    auto const [px, py, pz] =
        momentumComponents(theta / 180. * M_PI, phi / 180. * M_PI, P0);

    // define projectile kinematics in lab frame
    HEPMassType const projectileMass = 1_GeV;
    MomentumVector pProjectileLab{px, py, pz};
    HEPEnergyType const eProjectileLab = energy(projectileMass, pProjectileLab);
    FourVector const PprojLab(eProjectileLab, pProjectileLab);

    // define boost to com frame
    COMBoost boost(PprojLab, targetMass);

    // boost projecticle
    auto const PprojCoM = boost.toCoM(PprojLab);

    // boost target
    auto const PtargCoM = boost.toCoM(FourVector(targetMass, pTargetLab));

    // sum of momenta in CoM, should be 0
    auto const sumPCoM =
        PprojCoM.getSpaceLikeComponents() + PtargCoM.getSpaceLikeComponents();
    CHECK(sumPCoM.getNorm() / 1_GeV == Approx(0).margin(absMargin));
  }

  /*
    test the ultra-high energy behaviour: E=ZeV
   */

  SECTION("High energy") {
    // define projectile kinematics in lab frame
    HEPMassType const projectileMass = 1_GeV;
    HEPMomentumType P0 = 1_ZeV;
    MomentumVector pProjectileLab{0_GeV, 0_PeV, -P0};
    HEPEnergyType const eProjectileLab = energy(projectileMass, pProjectileLab);
    FourVector const PprojLab(eProjectileLab, pProjectileLab);

    // define boost to com frame
    COMBoost boost(PprojLab, targetMass);

    // boost projecticle
    auto const PprojCoM = boost.toCoM(PprojLab);

    // boost target
    auto const PtargCoM = boost.toCoM(FourVector(targetMass, pTargetLab));

    // sum of momenta in CoM, should be 0
    auto const sumPCoM =
        PprojCoM.getSpaceLikeComponents() + PtargCoM.getSpaceLikeComponents();
    CHECK(sumPCoM.getNorm() / P0 == Approx(0).margin(absMargin)); // MAKE RELATIVE CHECK
  }

  SECTION("CoM system") {

    MomentumVector pCM{0_GeV, 0_GeV, 5_GeV};

    COMBoost boostCMS({energy(1_GeV, pCM), pCM}, {energy(1_GeV, pCM), -pCM});

    auto test1 =
        boostCMS.fromCoM(FourMomentum{0_GeV, MomentumVector{0_GeV, 0_GeV, 0_GeV}});
    CHECK(test1.getNorm() == 0_GeV);
    auto test2 = boostCMS.toCoM(FourMomentum{0_GeV, MomentumVector{0_GeV, 0_GeV, 0_GeV}});
    CHECK(test2.getNorm() == 0_GeV);

    HEPEnergyType const sqrtS =
        (FourMomentum{energy(1_GeV, pCM), pCM} + FourMomentum{energy(1_GeV, pCM), -pCM})
            .getNorm();
    HEPEnergyType const eLab =
        (static_pow<2>(sqrtS) - 2 * static_pow<2>(1_GeV)) / (2 * 1_GeV);
    COMBoost boostLab({eLab, MomentumVector{momentum(eLab, 1_GeV), 0_eV, 0_eV}},
                      {1_GeV, MomentumVector{0_eV, 0_eV, 0_eV}});

    FourMomentum p4lab_trans(10_GeV, MomentumVector{0_eV, momentum(10_GeV, 1_GeV), 0_eV});
    FourMomentum p4lab_long(10_GeV, MomentumVector{momentum(10_GeV, 1_GeV), 0_GeV, 0_eV});
    // boost of transverse momentum
    CHECK(boostLab.toCoM(p4lab_trans).getNorm() / 1_GeV == Approx(1));
    CHECK(boostLab.toCoM(p4lab_trans).getTimeLikeComponent() / 1_GeV == Approx(50.99));
    // boost of longitudinal momentum
    CHECK(boostLab.toCoM(p4lab_long).getNorm() / 1_GeV == Approx(1));
    CHECK(boostLab.toCoM(p4lab_long).getTimeLikeComponent() / 1_GeV ==
          Approx(1.24).margin(0.1));
  }
}

TEST_CASE("rest frame") {

  logging::set_level(logging::level::info);

  HEPMassType const projectileMass = 1_GeV;
  HEPMomentumType const P0 = 1_TeV;
  MomentumVector pProjectileLab{{0_GeV, P0, 0_GeV}};
  HEPEnergyType const eProjectileLab = energy(projectileMass, pProjectileLab);
  const FourVector PprojLab(eProjectileLab, pProjectileLab);

  COMBoost boostRest(pProjectileLab, projectileMass);
  FourVector const rest4Mom = boostRest.toCoM(PprojLab);

  CHECK(rest4Mom.getTimeLikeComponent() / 1_GeV == Approx(projectileMass / 1_GeV));
  CHECK(rest4Mom.getSpaceLikeComponents().getNorm() / 1_GeV ==
        Approx(0).margin(absMargin));

  FourVector const a{0_eV, MomentumVector{0_eV, 5_GeV, 0_eV}};

  FourVector const bRoot{0_eV, MomentumVector{3_GeV, 0_eV, 0_eV}};

  // Since fromCoM() now expects rotated-frame input, we
  // explicitly pre-rotate root -> rotated before calling fromCoM, using
  // getInverseRotation() (root -> rotated)
  FourVector const b{
      0_eV, MomentumVector(QuantityVector<hepmomentum_d>(
                boostRest.getInverseRotation() *
                bRoot.getSpaceLikeComponents().getComponents().getEigenVector()))};

  auto const aLab = boostRest.fromCoM(a);
  auto const bLab = boostRest.fromCoM(b);

  CHECK(aLab.getNorm() / a.getNorm() == Approx(1));

  // recover the y-component of aLab as expressed in the rotated frame
  // (root -> rotated via getInverseRotation()), matching the old
  // aLab.getComponents(csPrime)[1] check
  auto const aLabRotated = QuantityVector<hepmomentum_d>(
      boostRest.getInverseRotation() *
      aLab.getSpaceLikeComponents().getComponents().getEigenVector());
  CHECK(aLabRotated[1].magnitude() == Approx((5_GeV).magnitude()));

  // bLab is in root CS (fromCoM output convention)
  CHECK(bLab.getSpaceLikeComponents().getComponents()[0].magnitude() ==
        Approx((3_GeV).magnitude()));
}
