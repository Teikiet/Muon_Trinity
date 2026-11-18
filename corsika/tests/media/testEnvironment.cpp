/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <corsika/media/density/DensityFunction.hpp>
#include <corsika/media/density_and_composition/FlatExponential.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/refractivity/UniformRefractiveIndex.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/interfaces/IMediumPropertyModel.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>
#include <corsika/media/density_and_composition/InhomogeneousMedium.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/density/integrator/LinearApproximationIntegrator.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/density_and_composition/SlidingPlanarExponential.hpp>
#include <corsika/media/density_and_composition/SlidingPlanarTabular.hpp>
#include <corsika/media/VolumeTreeNode.hpp>

#include <SetupTestTrajectory.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

Point const gOrigin(0_m, 0_m, 0_m);

TEST_CASE("VolumeTree") {
  logging::set_level(logging::level::info);
  media::Environment<IEmpty> env;
  auto& universe = *(env.getUniverse());
  auto world =
      media::Environment<IEmpty>::createNode<Sphere>(Point{0_m, 0_m, 0_m}, 150_km);
  // volume cut partly by "world"
  auto vol1 =
      media::Environment<IEmpty>::createNode<Sphere>(Point{0_m, 0_m, 140_km}, 20_km);
  // partly overlap with "vol1"
  auto vol2 =
      media::Environment<IEmpty>::createNode<Sphere>(Point{0_m, 0_m, 120_km}, 5_km);

  vol1->excludeOverlapWith(vol2);
  world->addChild(std::move(vol1));
  world->addChild(std::move(vol2));

  // in world
  CHECK(dynamic_cast<Sphere const&>(
            world->getContainingNode(Point(0_m, 0_m, 0_m))->getVolume())
            .getRadius() == 150_km);
  // in vol1
  CHECK(dynamic_cast<Sphere const&>(
            world->getContainingNode(Point(0_m, 0_m, 149_km))->getVolume())
            .getRadius() == 20_km);
  // outside world, in universe
  CHECK(world->getContainingNode(Point(0_m, 151_km, 0_m)) == nullptr);
  // in vol2
  CHECK(dynamic_cast<Sphere const&>(
            world->getContainingNode(Point(0_m, 0_km, 119_km))->getVolume())
            .getRadius() == 5_km);
  CHECK(dynamic_cast<Sphere const&>(
            world->getContainingNode(Point(0_m, 0_km, 121_km))->getVolume())
            .getRadius() == 5_km);

  // contained in world
  auto nestingPoint = Point{0_m, 0_m, -50_km};
  auto r3 = 3_km;
  auto vol3 = media::Environment<IEmpty>::createNode<Sphere>(nestingPoint, r3);
  world->addChildToContainingNode(nestingPoint, std::move(vol3));
  // check that vol3 has been added correctly
  CHECK(dynamic_cast<Sphere const&>(world->getContainingNode(nestingPoint)->getVolume())
            .getRadius() == r3);
  // check nesting of vol3 inside world
  CHECK(dynamic_cast<Sphere const&>(
            world->getContainingNode(nestingPoint)->getParent()->getVolume())
            .getRadius() == 150_km);

  auto r4 = 1_km;
  auto vol4 = media::Environment<IEmpty>::createNode<Sphere>(nestingPoint, r4);
  world->addChildToContainingNode(nestingPoint, std::move(vol4));
  // check that vol4 has been added correctly
  CHECK(dynamic_cast<Sphere const&>(world->getContainingNode(nestingPoint)->getVolume())
            .getRadius() == r4);
  // check nesting of vol4 inside vol3
  CHECK(dynamic_cast<Sphere const&>(
            world->getContainingNode(nestingPoint)->getParent()->getVolume())
            .getRadius() == r3);

  // throws if you try to add child outside of all existing nodes
  CHECK_THROWS(world->addChildToContainingNode(Point(0_m, 151_km, 0_m), std::move(vol3)));
  CHECK(world->getContainingNode(Point(0_m, 151_km, 0_m)) == nullptr);

  universe.addChild(std::move(world));
}

TEST_CASE("HomogeneousMedium") {

  logging::set_level(logging::level::info);

  media::NuclearComposition const protonComposition(std::vector<Code>{Code::Proton},
                                                    {1.});
  HomogeneousMedium<media::IMediumModel> const medium(19.2_g / cube(1_cm),
                                                      protonComposition);

  CHECK_THROWS(NuclearComposition({Code::Proton}, {1.1}));
  CHECK_THROWS(NuclearComposition({Code::Proton}, {0.99}));
}

TEST_CASE("FlatExponential") {

  logging::set_level(logging::level::info);

  media::NuclearComposition const protonComposition({Code::Proton}, {1.});

  Vector const axis(QuantityVector<dimensionless_d>(0, 0, 1));
  LengthType const lambda = 3_m;
  auto const rho0 = 1_g / cube(1_cm);
  FlatExponential<media::IMediumModel> const medium(gOrigin, axis, rho0, lambda,
                                                    protonComposition);
  SpeedType const speed = 20_m / second;
  LengthType const length = 2_m;
  TimeType const tEnd = length / speed;

  CHECK((medium.getMassDensity(gOrigin)) == rho0);
  CHECK(medium.getNuclearComposition().getFractions() == std::vector<double>{1.});
  CHECK(medium.getNuclearComposition().getComponents() ==
        std::vector<Code>{Code::Proton});

  SECTION("horizontal") {
    // Check that not moving along axis does not change density
    CHECK(medium.getMassDensity(Point(1_m, 0_m, 0_m)) == rho0);
    CHECK(medium.getMassDensity(Point(1_m, 0_m, 0_m)) == rho0);
    CHECK(medium.getMassDensity(Point(-1_m, 0_m, 0_m)) == rho0);
    CHECK(medium.getMassDensity(Point(0_m, 1_m, 0_m)) == rho0);
    CHECK(medium.getMassDensity(Point(0_m, -1_m, 0_m)) == rho0);

    Line const line(gOrigin,
                    Vector<SpeedType::dimension_type>(speed, 0_m / second, 0_m / second));
    setup::Trajectory const trajectory =
        setup::testing::make_track<setup::Trajectory>(line, tEnd);

    CHECK((medium.getIntegratedGrammage(trajectory) / (rho0 * length)) == Approx(1));
    CHECK((medium.getArclengthFromGrammage(trajectory, rho0 * length) / length) ==
          Approx(1));
  }

  SECTION("vertical") {
    // Moving along axis does change density
    CHECK(medium.getMassDensity(Point(0_m, 0_m, 1_m)) > rho0);
    CHECK(medium.getMassDensity(Point(0_m, 0_m, -1_m)) < rho0);

    Line const line(gOrigin,
                    Vector<SpeedType::dimension_type>(0_m / second, 0_m / second, speed));
    setup::Trajectory const trajectory =
        setup::testing::make_track<setup::Trajectory>(line, tEnd);
    GrammageType const exact = rho0 * lambda * (exp(length / lambda) - 1);

    CHECK((medium.getIntegratedGrammage(trajectory) / exact) == Approx(1));
    CHECK((medium.getArclengthFromGrammage(trajectory, exact) / length) == Approx(1));
  }

  SECTION("escape grammage") {
    Line const line(gOrigin, Vector<SpeedType::dimension_type>(
                                 SpeedType::zero(), SpeedType::zero(), -speed));
    setup::Trajectory const trajectory =
        setup::testing::make_track<setup::Trajectory>(line, tEnd);
    GrammageType const escapeGrammage = rho0 * lambda;

    CHECK(trajectory.getDirection(0).dot(axis).magnitude() < 0);
    CHECK(medium.getArclengthFromGrammage(trajectory, 1.2 * escapeGrammage) ==
          std::numeric_limits<typename GrammageType::value_type>::infinity() * 1_m);
  }

  SECTION("inclined") {
    Line const line(gOrigin, Vector<SpeedType::dimension_type>(
                                 0_m / second, speed / sqrt(2.), speed / sqrt(2.)));
    setup::Trajectory const trajectory =
        setup::testing::make_track<setup::Trajectory>(line, tEnd);
    double const cosTheta = M_SQRT1_2;
    GrammageType const exact =
        rho0 * lambda * (exp(cosTheta * length / lambda) - 1) / cosTheta;
    CHECK((medium.getIntegratedGrammage(trajectory) / exact) == Approx(1));
    CHECK((medium.getArclengthFromGrammage(trajectory, exact) / length) == Approx(1));
  }
}

TEST_CASE("SlidingPlanarExponential") {

  logging::set_level(logging::level::info);

  media::NuclearComposition const protonComposition(std::vector<Code>{Code::Proton},
                                                    {1.});

  LengthType const lambda = 3_m;
  auto const rho0 = 1_g / static_pow<3>(1_cm);
  auto const tEnd = 5_s;

  SlidingPlanarExponential<media::IMediumModel> const medium(gOrigin, rho0, lambda,
                                                             protonComposition);

  SECTION("density") {
    CHECK(medium.getMassDensity({0_m, 0_m, 3_m}) /
              medium.getMassDensity({0_m, 3_m, 0_m}) ==
          Approx(1));
  }

  SECTION("vertical") {
    Vector const axis(QuantityVector<dimensionless_d>(0, 0, 1));
    FlatExponential<media::IMediumModel> const flat(gOrigin, axis, rho0, lambda,
                                                    protonComposition);
    Line const line({0_m, 0_m, 1_m}, Vector<SpeedType::dimension_type>(
                                         0_m / second, 0_m / second, 5_m / second));
    setup::Trajectory const trajectory =
        setup::testing::make_track<setup::Trajectory>(line, tEnd);

    CHECK(medium.getMassDensity({0_mm, 0_m, 3_m}).magnitude() ==
          flat.getMassDensity({0_mm, 0_m, 3_m}).magnitude());
    CHECK(medium.getIntegratedGrammage(trajectory).magnitude() ==
          flat.getIntegratedGrammage(trajectory).magnitude());
    CHECK(medium.getArclengthFromGrammage(trajectory, rho0 * 5_m).magnitude() ==
          flat.getArclengthFromGrammage(trajectory, rho0 * 5_m).magnitude());
  }
}

struct RhoFuncConst {
  MassDensityType operator()(LengthType) const { return 1_g / cube(1_cm); }
  static GrammageType integrate(LengthType dL) { return dL * 1_g / cube(1_cm); }
};

struct RhoFuncExp {
  MassDensityType operator()(LengthType height) const {
    return 1_g / cube(1_cm) * exp(-height / 1000_m);
  }
  static GrammageType integrate(BaseTrajectory const& traj, Point const& origin,
                                LengthType const& refH) {
    LengthType height1 = (traj.getPosition(0) - origin).getNorm() - refH;
    LengthType height2 = (traj.getPosition(1) - origin).getNorm() - refH;
    if (height1 > height2) { std::swap(height1, height2); }

    DirectionVector const axis(
        (traj.getPosition(0) - origin).normalized()); // to gravity center
    double const cosTheta = axis.dot(traj.getDirection(0));

    CORSIKA_LOG_INFO("h1={} h2={} cT={} rho1={}, rho2={}", height1, height2, cosTheta,
                     1_g / cube(1_cm) * exp(-height1 / 1000_m),
                     1_g / cube(1_cm) * exp(-height2 / 1000_m));
    return (1_km * 1_g / cube(1_cm) * exp(-height1 / 1000_m) -
            1_km * 1_g / cube(1_cm) * exp(-height2 / 1000_m)) /
           cosTheta;
  }
  static GrammageType integrate(BaseTrajectory const& traj, LengthType const& length,
                                Point const& origin, LengthType const& refH) {
    LengthType height1 = (traj.getPosition(0) - origin).getNorm() - refH;
    LengthType height2 =
        (traj.getPosition(0) + traj.getDirection(0) * length - origin).getNorm() - refH;
    if (height1 > height2) { std::swap(height1, height2); }

    DirectionVector const axis(
        (traj.getPosition(0) - origin).normalized()); // to gravity center
    double const cosTheta = axis.dot(traj.getDirection(0));

    CORSIKA_LOG_INFO("h1={} h2={} cT={}", height1, height2, cosTheta);
    return (1_km * 1_g / cube(1_cm) * exp(-height1 / 1000_m) -
            1_km * 1_g / cube(1_cm) * exp(-height2 / 1000_m)) /
           cosTheta;
  }
};

TEST_CASE("SlidingPlanarTabular") {

  logging::set_level(logging::level::info);

  media::NuclearComposition const protonComposition(std::vector<Code>{Code::Proton},
                                                    {1.});

  RhoFuncConst rhoFunc;
  SlidingPlanarTabular<media::IMediumModel> const medium(gOrigin, rhoFunc, 1000, 10_m,
                                                         protonComposition);

  SECTION("not possible") {
    CHECK_THROWS(medium.getMassDensity({0_m, 1e10_m, 0_m}));

    SpeedType const speed = 5_m / second;
    TimeType const tEnd = 1e10_s;
    Line const line({0_m, 0_m, 1_m},
                    Vector<SpeedType::dimension_type>(0_m / second, 0_m / second, speed));
    setup::Trajectory const trajectory =
        setup::testing::make_track<setup::Trajectory>(line, tEnd);
    CHECK_THROWS(medium.getIntegratedGrammage(trajectory));

    Line const line2({0_m, 0_m, 1e9_m}, Vector<SpeedType::dimension_type>(
                                            0_m / second, 0_m / second, speed));
    setup::Trajectory const trajectory2 =
        setup::testing::make_track<setup::Trajectory>(line2, tEnd);
    CHECK_THROWS(medium.getArclengthFromGrammage(trajectory2, 1e3_g / square(1_cm)));
  }

  SECTION("density") {
    CHECK(medium.getMassDensity({0_m, 0_m, 3_m}) /
              medium.getMassDensity({0_m, 3_m, 0_m}) ==
          Approx(1));
    CHECK(medium.getMassDensity({0_mm, 0_m, 3_m}) == 1_g / cube(1_cm));
    CHECK(medium.getMassDensity({0_mm, 0_m, 300_m}) == 1_g / cube(1_cm));
  }

  SECTION("vertical") {
    SpeedType const speed = 5_m / second;
    TimeType const tEnd1 = 1_s;
    LengthType const length1 = speed * tEnd1;
    TimeType const tEnd2 = 300_s;
    LengthType const length2 = speed * tEnd2;
    Line const line({0_m, 0_m, 1_m},
                    Vector<SpeedType::dimension_type>(0_m / second, 0_m / second, speed));
    setup::Trajectory const trajectory1 =
        setup::testing::make_track<setup::Trajectory>(line, tEnd1);
    Line const line1Reverse(
        trajectory1.getPosition(1),
        Vector<SpeedType::dimension_type>(0_m / second, 0_m / second, -speed));
    setup::Trajectory const trajectory1Reverse =
        setup::testing::make_track<setup::Trajectory>(line1Reverse, tEnd1);

    setup::Trajectory const trajectory2 =
        setup::testing::make_track<setup::Trajectory>(line, tEnd2);
    Line const line2Reverse(
        trajectory2.getPosition(0),
        Vector<SpeedType::dimension_type>(0_m / second, 0_m / second, -speed));

    // failures
    CHECK_THROWS(medium.getArclengthFromGrammage(trajectory1, -1_kg / square(1_cm)));

    MassDensityType const rho0 = 1_g / cube(1_cm);

    // short track
    CHECK(medium.getIntegratedGrammage(trajectory1) == length1 * rho0);
    LengthType const testD1 = length1 / 200; // within bin
    CHECK(medium.getArclengthFromGrammage(trajectory1, rho0 * testD1) / testD1 ==
          Approx(1));
    // short track, reverse
    CHECK(medium.getIntegratedGrammage(trajectory1Reverse) == length1 * rho0);
    CHECK(medium.getArclengthFromGrammage(trajectory1Reverse, rho0 * testD1) / testD1 ==
          Approx(1));

    // long track
    CHECK(medium.getIntegratedGrammage(trajectory2) == length2 * 1_g / cube(1_cm));
    LengthType const testD2 = length2 / 25; // multi bin
    CHECK(medium.getArclengthFromGrammage(trajectory2, rho0 * testD2) == testD2);
  }

  SECTION("inclined") {
    SpeedType const speed = 5_m / second;
    TimeType const tEnd1 = 1_s;
    LengthType const length1 = speed * tEnd1;
    TimeType const tEnd2 = 300_s;
    LengthType const length2 = speed * tEnd2;
    Line const line({0_m, 0_m, 1_m},
                    Vector<SpeedType::dimension_type>(speed / sqrt(2.), 0_m / second,
                                                      speed / sqrt(2.)));
    setup::Trajectory const trajectory1 =
        setup::testing::make_track<setup::Trajectory>(line, tEnd1);
    Line const line1Reverse(trajectory1.getPosition(1),
                            Vector<SpeedType::dimension_type>(
                                -speed / sqrt(2.), 0_m / second, -speed / sqrt(2.)));
    setup::Trajectory const trajectory1Reverse =
        setup::testing::make_track<setup::Trajectory>(line1Reverse, tEnd1);

    setup::Trajectory const trajectory2 =
        setup::testing::make_track<setup::Trajectory>(line, tEnd2);
    Line const line2Reverse(trajectory2.getPosition(1),
                            Vector<SpeedType::dimension_type>(
                                -speed / sqrt(2.), 0_m / second, -speed / sqrt(2.)));
    setup::Trajectory const trajectory2Reverse =
        setup::testing::make_track<setup::Trajectory>(line2Reverse, tEnd2);

    MassDensityType const rho0 = 1_g / cube(1_cm);

    // short track
    CHECK(medium.getIntegratedGrammage(trajectory1) / (length1 * rho0) == Approx(1));
    LengthType const testD1 = length1 / 200; // within bin
    CHECK(medium.getArclengthFromGrammage(trajectory1, RhoFuncConst::integrate(testD1)) /
              testD1 ==
          Approx(1));
    // short track, reverse
    CHECK(medium.getIntegratedGrammage(trajectory1Reverse) / (length1 * rho0) ==
          Approx(1));
    CHECK(medium.getArclengthFromGrammage(trajectory1Reverse, rho0 * testD1) / testD1 ==
          Approx(1));

    // long track
    CHECK(medium.getIntegratedGrammage(trajectory2) / (length2 * rho0) ==
          Approx(1).epsilon(0.01));
    LengthType const testD2 = length2 / 25; // multi bin
    CHECK(medium.getArclengthFromGrammage(trajectory2, rho0 * testD2) / testD2 ==
          Approx(1).epsilon(0.01));
    // long track reverse
    CORSIKA_LOG_INFO("length2={}", length2);
    CHECK(medium.getIntegratedGrammage(trajectory2Reverse) / (length2 * rho0) ==
          Approx(1).epsilon(0.01));
    CHECK(medium.getArclengthFromGrammage(trajectory2Reverse, rho0 * testD2) / testD2 ==
          Approx(1).epsilon(0.01));
  }

  /*The exponential test is taken over phase-space where the exponential is not so steep
   * and is samples in sufficient substeps. An reference-height offset of 1000_km is used.
   * Thus, density is given from 1000 to 1010 km. And curvature effects are small.
   */

  RhoFuncExp rhoFuncExp;
  SlidingPlanarTabular<media::IMediumModel> const mediumExp(
      gOrigin, rhoFuncExp, 1000, 10_m, protonComposition, 1000_km);

  SECTION("exponential") {

    SpeedType const speed = 5_m / second;
    TimeType const tEnd1 = 1_s;
    LengthType const length1 = speed * tEnd1;
    TimeType const tEnd2 = 300_s;
    LengthType const length2 = speed * tEnd2;
    Line const line({0_m, 0_m, 1000.005_km},
                    Vector<SpeedType::dimension_type>(speed / sqrt(2.), 0_m / second,
                                                      speed / sqrt(2.)));
    setup::Trajectory const trajectory1 =
        setup::testing::make_track<setup::Trajectory>(line, tEnd1);
    Line const line1Reverse(trajectory1.getPosition(1),
                            Vector<SpeedType::dimension_type>(
                                -speed / sqrt(2.), 0_m / second, -speed / sqrt(2.)));
    setup::Trajectory const trajectory1Reverse =
        setup::testing::make_track<setup::Trajectory>(line1Reverse, tEnd1);

    setup::Trajectory const trajectory2 =
        setup::testing::make_track<setup::Trajectory>(line, tEnd2);

    CORSIKA_LOG_INFO("{} {}", RhoFuncExp::integrate(trajectory1, gOrigin, 1000_km),
                     length1);

    // short track
    GrammageType const testShortX = RhoFuncExp::integrate(trajectory1, gOrigin, 1000_km);
    CHECK(mediumExp.getIntegratedGrammage(trajectory1) / testShortX ==
          Approx(1).epsilon(0.01));
    LengthType const testD1 = length1 / 200; // within bin
    GrammageType const testD1X =
        RhoFuncExp::integrate(trajectory1, testD1, gOrigin, 1000_km);
    CHECK(mediumExp.getArclengthFromGrammage(trajectory1, testD1X) / testD1 ==
          Approx(1).epsilon(0.01));
    // short track, reverse
    CHECK(mediumExp.getIntegratedGrammage(trajectory1Reverse) / testShortX ==
          Approx(1).epsilon(0.01));
    CHECK(mediumExp.getArclengthFromGrammage(trajectory1Reverse, testD1X) / testD1 ==
          Approx(1).epsilon(0.01));

    // long track
    GrammageType const testLongX = RhoFuncExp::integrate(trajectory2, gOrigin, 1000_km);
    CORSIKA_LOG_INFO("testLongX={}", testLongX);
    CHECK(mediumExp.getIntegratedGrammage(trajectory2) / testLongX ==
          Approx(1).epsilon(0.01));
    LengthType const testD2 = length2 / 25; // multi bin
    GrammageType const testD2X =
        RhoFuncExp::integrate(trajectory2, testD2, gOrigin, 1000_km);
    CHECK(mediumExp.getArclengthFromGrammage(trajectory2, testD2X) / testD2 ==
          Approx(1).epsilon(0.01));
    // long track, reverse

    // first full trajectory2 reverse
    Line line2Reverse(trajectory2.getPosition(1),
                      Vector<SpeedType::dimension_type>(-speed / sqrt(2.), 0_m / second,
                                                        -speed / sqrt(2.)));
    setup::Trajectory trajectory2Reverse =
        setup::testing::make_track<setup::Trajectory>(line2Reverse, tEnd2);

    CHECK(mediumExp.getIntegratedGrammage(trajectory2Reverse) / testLongX ==
          Approx(1).epsilon(0.01));

    // but now shorter trajectory2 reversed to correspond 100% to testD2

    line2Reverse = Line(trajectory2.getPosition(0) + trajectory2.getDirection(0) * testD2,
                        Vector<SpeedType::dimension_type>(-speed / sqrt(2.), 0_m / second,
                                                          -speed / sqrt(2.)));
    auto const trajectory2ReverseShort =
        setup::testing::make_track<setup::Trajectory>(line2Reverse, testD2 / speed);

    CORSIKA_LOG_INFO("here {} {} {}", trajectory2ReverseShort.getLength(), testD2,
                     testD2X / 1_g * square(1_cm));
    CHECK(mediumExp.getArclengthFromGrammage(trajectory2ReverseShort, testD2X) / testD2 ==
          Approx(1).epsilon(0.01));
  }
}

MassDensityType constexpr rho0 = 1_kg / 1_m / 1_m / 1_m;

struct ExponentialTest {
  auto operator()(Point const& p) const {
    return exp(p.getCoordinates()[0] / 1_m) * rho0;
  }

  template <int N>
  auto getDerivative(Point const& p, DirectionVector const& v) const {
    return v.getComponents()[0] * (*this)(p) / static_pow<N>(1_m);
  }

  auto getFirstDerivative(Point const& p, DirectionVector const& v) const {
    return getDerivative<1>(p, v);
  }

  auto getSecondDerivative(Point const& p, DirectionVector const& v) const {
    return getDerivative<2>(p, v);
  }
};

TEST_CASE("InhomogeneousMedium") {

  logging::set_level(logging::level::info);

  Vector direction(QuantityVector<dimensionless_d>(1, 0, 0));

  SpeedType const speed = 20_m / second;
  Line line(gOrigin, Vector<SpeedType::dimension_type>(speed, SpeedType::zero(),
                                                       SpeedType::zero()));

  // the tested LinearApproximationIntegrator really does a single step only. It is very
  // poor for exponentials with a bit larger step-width.
  TimeType const tEnd = 0.001_s;
  setup::Trajectory const trajectory =
      setup::testing::make_track<setup::Trajectory>(line, tEnd);

  ExponentialTest const expTest;
  DensityFunction<ExponentialTest, LinearApproximationIntegrator> const rho(expTest);

  SECTION("DensityFunction") {
    CHECK(expTest.getDerivative<1>(gOrigin, direction) / (1_kg / 1_m / 1_m / 1_m / 1_m) ==
          Approx(1));
    CHECK(rho.evaluateAt(gOrigin) == expTest(gOrigin));
  }

  auto const exactGrammage = [](auto l) { return 1_m * rho0 * (exp(l / 1_m) - 1); };
  auto const exactLength = [](auto X) { return 1_m * log(1 + X / (rho0 * 1_m)); };

  LengthType const length = tEnd * speed;

  media::NuclearComposition const composition{{Code::Proton}, {1.}};
  InhomogeneousMedium<IMediumModel, decltype(rho)> const inhMedium(composition, rho);

  CORSIKA_LOG_INFO("test={} l={} {} {}", rho.getIntegrateGrammage(trajectory), length,
                   exactGrammage(length), 1_m * rho0 * (exp(length / 1_m) - 1));

  SECTION("Integration") {
    CORSIKA_LOG_INFO("test={} {} {}", rho.getIntegrateGrammage(trajectory),
                     exactGrammage(length),
                     rho.getIntegrateGrammage(trajectory) / exactGrammage(length));
    CHECK(rho.getIntegrateGrammage(trajectory) / exactGrammage(length) ==
          Approx(1).epsilon(1e-2));
    CHECK(rho.getArclengthFromGrammage(trajectory, exactGrammage(length)) /
              exactLength(exactGrammage(length)) ==
          Approx(1).epsilon(1e-2));
    CHECK(rho.getMaximumLength(trajectory, 1e-2) >
          length); // todo: write reasonable test when implementation is working

    CHECK(rho.getIntegrateGrammage(trajectory) ==
          inhMedium.getIntegratedGrammage(trajectory));
    CHECK(rho.getArclengthFromGrammage(trajectory, 20_g / (1_cm * 1_cm)) ==
          inhMedium.getArclengthFromGrammage(trajectory, 20_g / (1_cm * 1_cm)));
    CHECK(inhMedium.getNuclearComposition() == composition);
    CHECK(inhMedium.getMassDensity({0_m, 0_m, 0_m}) == 1_kg / static_pow<3>(1_m));
  }
}

TEST_CASE("LayeredSphericalAtmosphereBuilder") {

  logging::set_level(logging::level::info);

  LayeredSphericalAtmosphereBuilder builder =
      make_layered_spherical_atmosphere_builder<>::create(gOrigin,
                                                          constants::EarthRadius::Mean);

  builder.setNuclearComposition({{{Code::Nitrogen, Code::Oxygen}}, {{.6, .4}}});

  builder.addLinearLayer(1_g / (1_cm * 1_cm), 1_km, 10_km);
  builder.addLinearLayer(1_g / (1_cm * 1_cm), 2_km, 20_km);
  builder.addExponentialLayer(540.1778_g / (1_cm * 1_cm), 772170.16_cm, 30_km);

  CHECK_THROWS(builder.addLinearLayer(1_g / (1_cm * 1_cm), 0.5_km, 5_km));

  CHECK(builder.getSize() == 3);

  auto const builtEnv = builder.assemble();
  auto const& univ = builtEnv.getUniverse();

  CHECK(builder.getSize() == 0);

  auto const R = builder.getPlanetRadius();

  CHECK(univ->getChildNodes().size() == 1);

  CHECK(univ->getContainingNode(Point(0_m, 0_m, R + 35_km)) == univ.get());
  CHECK(dynamic_cast<Sphere const&>(
            univ->getContainingNode(Point(0_m, 0_m, R + 8_km))->getVolume())
            .getRadius() == R + 10_km);
  CHECK(dynamic_cast<Sphere const&>(
            univ->getContainingNode(Point(0_m, 0_m, R + 12_km))->getVolume())
            .getRadius() == R + 20_km);
  CHECK(dynamic_cast<Sphere const&>(
            univ->getContainingNode(Point(0_m, 0_m, R + 24_km))->getVolume())
            .getRadius() == R + 30_km);

  CHECK(corsika::media::get_all_elements_in_universe(builtEnv).size() > 0);
}

TEST_CASE("LayeredSphericalAtmosphereBuilder w/ magnetic field") {

  logging::set_level(logging::level::info);

  // setup our interface types
  using ModelInterface = media::IMagneticFieldModel<media::IMediumModel>;

  // the composition we use for the homogenous medium
  media::NuclearComposition const protonComposition(std::vector<Code>{Code::Proton},
                                                    {1.});

  // create magnetic field vectors
  Vector B0(0_T, 0_T, 1_T);

  LayeredSphericalAtmosphereBuilder builder =
      media::make_layered_spherical_atmosphere_builder<
          ModelInterface, UniformMagneticField>::create(gOrigin,
                                                        constants::EarthRadius::Mean, B0);

  builder.setNuclearComposition({{{Code::Nitrogen, Code::Oxygen}}, {{.6, .4}}});
  builder.addLinearLayer(1_g / (1_cm * 1_cm), 1_km, 10_km);
  builder.addExponentialLayer(1222.6562_g / (1_cm * 1_cm), 994186.38_cm, 20_km);

  CHECK(builder.getSize() == 2);

  auto const builtEnv = builder.assemble();
  auto const& univ = builtEnv.getUniverse();

  CHECK(builder.getSize() == 0);
  CHECK(univ->getChildNodes().size() == 1);
  auto const R = builder.getPlanetRadius();

  // check magnetic field at several locations
  const Point pTest(-10_m, 4_m, R + 35_m);
  CHECK(B0.getComponents() == univ->getContainingNode(pTest)
                                  ->getModelProperties()
                                  .getMagneticField(pTest)
                                  .getComponents());
  const Point pTest2(10_m, -4_m, R + 15_km);
  CHECK(B0.getComponents() == univ->getContainingNode(pTest2)
                                  ->getModelProperties()
                                  .getMagneticField(pTest2)
                                  .getComponents());
}

TEST_CASE("media", "LayeredSphericalAtmosphereBuilder USStd") {
  // setup environment, geometry
  Point const center{0_m, 0_m, 0_m};

  // setup our interface types
  auto builder = media::make_layered_spherical_atmosphere_builder<>::create(
      center, constants::EarthRadius::Mean);

  builder.setNuclearComposition(
      {{Code::Nitrogen, Code::Oxygen},
       {0.7847f, 1.f - 0.7847f}}); // values taken from AIRES manual, Ar removed for now

  builder.addExponentialLayer(1222.6562_g / (1_cm * 1_cm), 994186.38_cm, 4_km);
  builder.addExponentialLayer(1144.9069_g / (1_cm * 1_cm), 878153.55_cm, 10_km);
  builder.addExponentialLayer(1305.5948_g / (1_cm * 1_cm), 636143.04_cm, 40_km);
  builder.addExponentialLayer(540.1778_g / (1_cm * 1_cm), 772170.16_cm, 100_km);
  builder.addLinearLayer(1_g / (1_cm * 1_cm), 1e9_cm, 112.8_km);

  media::Environment<media::IMediumModel> env;
  builder.assemble(env);

  typedef
      typename media::Environment<media::IMediumModel>::BaseNodeType::VTN_type node_type;
  node_type const* universe = env.getUniverse().get();

  // far out there is the universe
  CHECK(universe->getContainingNode(Point(10000_km, 0_m, 0_m)) == universe);
  CHECK(universe->getContainingNode(Point(0_m, 10000_km, 0_m)) == universe);

  // at 112.8km there is transition to atmosphere
  CHECK(universe->getContainingNode(
            Point(constants::EarthRadius::Mean + 112.8_km + 1_cm, 0_m, 0_m)) == universe);
  CHECK(universe->getContainingNode(
            Point(0_m, constants::EarthRadius::Mean + 112.8_km + 1_cm, 0_m)) == universe);

  // check layer transition at 112.8km

  node_type const* layer1_not_yet = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 112.8_km + 1_cm, 0_m, 0_m));
  node_type const* layer1 = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 112.8_km - 1_cm, 0_m, 0_m));
  node_type const* layer1_also = universe->getContainingNode(
      Point(0_m, constants::EarthRadius::Mean + 112.8_km - 1_cm, 0_m));

  CHECK(layer1_not_yet == universe);
  CHECK(layer1 != universe);
  CHECK(layer1 == layer1_also);

  // check layer transition at 100km

  node_type const* layer2_not_yet = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 100_km + 1_cm, 0_m, 0_m));
  node_type const* layer2 = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 100_km - 1_cm, 0_m, 0_m));
  node_type const* layer2_also = universe->getContainingNode(
      Point(0_m, constants::EarthRadius::Mean + 100_km - 1_cm, 0_m));

  CHECK(layer2_not_yet == layer1);
  CHECK(layer2 != layer1);
  CHECK(layer2 == layer2_also);

  // check layer transition at 40km

  node_type const* layer3_not_yet = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 40_km + 1_cm, 0_m, 0_m));
  node_type const* layer3 = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 40_km - 1_cm, 0_m, 0_m));
  node_type const* layer3_also = universe->getContainingNode(
      Point(0_m, constants::EarthRadius::Mean + 40_km - 1_cm, 0_m));

  CHECK(layer3_not_yet == layer2);
  CHECK(layer3 != layer2);
  CHECK(layer3 == layer3_also);

  // check layer transition at 10km

  node_type const* layer4_not_yet = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 10_km + 1_cm, 0_m, 0_m));
  node_type const* layer4 = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 10_km - 1_cm, 0_m, 0_m));
  node_type const* layer4_also = universe->getContainingNode(
      Point(0_m, constants::EarthRadius::Mean + 10_km - 1_cm, 0_m));

  CHECK(layer4_not_yet == layer3);
  CHECK(layer4 != layer3);
  CHECK(layer4 == layer4_also);

  // check layer transition at 4km

  node_type const* layer5_not_yet = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 4_km + 1_cm, 0_m, 0_m));
  node_type const* layer5 = universe->getContainingNode(
      Point(constants::EarthRadius::Mean + 4_km - 1_cm, 0_m, 0_m));
  node_type const* layer5_also = universe->getContainingNode(
      Point(0_m, constants::EarthRadius::Mean + 4_km - 1_cm, 0_m));

  CHECK(layer5_not_yet == layer4);
  CHECK(layer5 != layer4);
  CHECK(layer5 == layer5_also);
}
