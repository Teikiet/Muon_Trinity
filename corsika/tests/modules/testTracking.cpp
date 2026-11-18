/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/Tracking.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>
#include <SetupTestTrajectory.hpp>

#include <catch2/catch_all.hpp>

#include <boost/type_index.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

using DummyEnvironmentInterface = IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>;
using DummyEnvironment = Environment<DummyEnvironmentInterface>;

struct NonExistingDummyObject : public IVolume {
  NonExistingDummyObject const& getVolume() const { return *this; }
  bool contains(Point const&) const { return false; }
};

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/**
 * @file testTracking.cpp
 *
 * This is the unified and common unit test for all Tracking algorithms:
 *
 * - tracking_leapfrog_curved::Tracking
 * - tracking_leapfrog_straight::Tracking
 * - tracking_line::Tracking
 *
 *
 * The main part of tests are to inject particles at 10GeV momentum at
 * (-Rg,0,0) in +x direction into a sphere of radius Rg, where Rg is
 * the gyroradius (or 10m for neutral particles). Then it is checked
 * where the particles leaves the sphere for different charges
 * (-1,0,+1) and field strength (-50uT, 0T, +50uT).
 *
 * Each test is perfromed once, with the particles starting logically
 * outside of the Rg sphere (thus it first has to enter insides) and a
 * second time with the particle already logically inside the sphere.
 *
 * There is a second smaller, internal sphere at +z displacement. Only
 * neutral particles are allowed and expected to hit this.
 *
 * All those tests are parameterized, thus, they can be easily extended
 * or applied to new algorithms.
 */

TEMPLATE_TEST_CASE("Tracking", "tracking", tracking_leapfrog_curved::Tracking,
                   tracking_leapfrog_straight::Tracking, tracking_line::Tracking) {

  logging::set_level(logging::level::info);

  const HEPEnergyType P0 = 10_GeV;

  auto PID = GENERATE(as<Code>{}, Code::MuPlus, Code::MuPlus, Code::Photon);
  // test also special case: movement parallel to field (along x)
  auto isParallel = GENERATE(as<bool>{}, true, false);
  // for algorithms that know magnetic deflections choose: +-50uT, 0uT
  // otherwise just 0uT
  auto Bfield = GENERATE_COPY(filter(
      []([[maybe_unused]] MagneticFluxType v) {
        if constexpr (std::is_same_v<TestType, tracking_line::Tracking>)
          return v == 0_uT;
        else
          return true;
      },
      values<MagneticFluxType>({50_uT, 0_uT, (isParallel ? 0_uT : -50_uT)})));
  // particle --> (world) --> | --> (target)
  // true: start inside "world" volume
  // false: start inside "target" volume
  auto outer = GENERATE(as<bool>{}, true, false);

  SECTION(fmt::format("Tracking PID={}, Bfield={} uT, isParallel={}, from outside={}",
                      PID, Bfield / 1_uT, isParallel, outer)) {

    CORSIKA_LOG_DEBUG(
        "********************\n                          TEST algo={} section PID={}, "
        "Bfield={} "
        "uT, field is parallel={}, start_outside={}",
        boost::typeindex::type_id<TestType>().pretty_name(), PID, Bfield / 1_uT,
        isParallel, outer);

    const int chargeNumber = get_charge_number(PID);
    LengthType radius = 10_m;
    int deflect = 0;
    if (chargeNumber != 0 && Bfield != 0_T && !isParallel) {
      deflect = -sgn(chargeNumber) * sgn(Bfield / 1_T); // direction of deflection
      LengthType const gyroradius = (convert_HEP_to_SI<MassType::dimension_type>(P0) *
                                     constants::c / (abs(get_charge(PID)) * abs(Bfield)));
      CORSIKA_LOG_DEBUG("Rg={} deflect={}", gyroradius, deflect);
      radius = gyroradius;
    }

    auto [env, csPtr, worldPtr] =
        corsika::setup::testing::setup_environment(Code::Oxygen, Bfield);
    { [[maybe_unused]] const auto& env_dummy = env; }
    auto const& cs = *csPtr;

    TestType tracking;
    Point const center(cs, {0_m, 0_m, 0_m});
    auto target = DummyEnvironment::createNode<Sphere>(center, radius);

    // every particle should hit target_2
    // it is very close to injection and not so small
    auto target_2 = DummyEnvironment::createNode<Sphere>(
        Point(cs, {-radius * 3 / 4, 0_m, 0_m}), radius * 0.2);

    // only neutral particles hit_target_neutral
    // this is far from injection and really small
    auto target_neutral = DummyEnvironment::createNode<Sphere>(
        Point(cs, {radius / 2, 0_m, 0_m}), radius * 0.1);

    // target to be overlapped entirely by target_2
    auto target_2_behind = DummyEnvironment::createNode<Sphere>(
        Point(cs, {-radius * 3 / 4, 0_m, 0_m}), radius * 0.1);

    // target to be overlapped partly by target_2
    auto target_2_partly_behind = DummyEnvironment::createNode<Sphere>(
        Point(cs, {-radius * 3 / 4 + radius * 0.1, 0_m, 0_m}), radius * 0.2);

    using MyHomogeneousModel = media::MediumPropertyModel<
        UniformMagneticField<media::HomogeneousMedium<DummyEnvironmentInterface>>>;

    MagneticFieldVector magneticfield(cs, 0_T, 0_T, Bfield);
    target->setModelProperties<MyHomogeneousModel>(
        Medium::AirDry1Atm, magneticfield, 1_g / (1_m * 1_m * 1_m),
        media::NuclearComposition(std::vector<Code>{Code::Oxygen},
                                  std::vector<double>{1.}));
    target_neutral->setModelProperties<MyHomogeneousModel>(
        Medium::AirDry1Atm, magneticfield, 1_g / (1_m * 1_m * 1_m),
        media::NuclearComposition(std::vector<Code>{Code::Oxygen},
                                  std::vector<double>{1.}));
    target_2->setModelProperties<MyHomogeneousModel>(
        Medium::AirDry1Atm, magneticfield, 1_g / (1_m * 1_m * 1_m),
        media::NuclearComposition(std::vector<Code>{Code::Oxygen},
                                  std::vector<double>{1.}));
    target_2_behind->setModelProperties<MyHomogeneousModel>(
        Medium::AirDry1Atm, magneticfield, 1_g / (1_m * 1_m * 1_m),
        media::NuclearComposition(std::vector<Code>{Code::Oxygen},
                                  std::vector<double>{1.}));
    target_2_partly_behind->setModelProperties<MyHomogeneousModel>(
        Medium::AirDry1Atm, magneticfield, 1_g / (1_m * 1_m * 1_m),
        media::NuclearComposition(std::vector<Code>{Code::Oxygen},
                                  std::vector<double>{1.}));
    auto* targetPtr = target.get();
    auto* targetPtr_2 = target_2.get();
    auto* targetPtr_neutral = target_neutral.get();
    auto* targetPtr_2_behind = target_2_behind.get();
    auto* targetPtr_2_partly_behind = target_2_partly_behind.get();
    target_2_behind->excludeOverlapWith(target_2);
    target_2_partly_behind->excludeOverlapWith(target_2);
    worldPtr->addChild(std::move(target));
    targetPtr->addChild(std::move(target_2));
    targetPtr->addChild(std::move(target_neutral));
    targetPtr->addChild(std::move(target_2_behind));
    targetPtr->addChild(std::move(target_2_partly_behind));

    auto [stack, viewPtr] = setup::testing::setup_stack(PID, P0, targetPtr, cs);
    { [[maybe_unused]] auto& viewPtr_dum = viewPtr; }
    auto particle = stack->first();
    // Note: momentum in X-direction
    //       magnetic field in Z-direction
    //       put particle on x_start, 0, 0
    //       expect intersections somewere in +-y_start

    if (outer) {
      particle.setNode(worldPtr); // set particle inside "target" volume
    } else {
      particle.setNode(targetPtr); // set particle outside "target" volume
    }
    particle.setPosition(Point(cs, -radius, 0_m, 0_m));

    auto [traj, nextVol] = tracking.getTrack(particle);
    particle.setNode(nextVol);
    particle.setPosition(traj.getPosition(1));
    particle.setDirection(traj.getDirection(1));
    SpeedType const speed_0 = particle.getVelocity().getNorm();
    if (outer) {
      // now we know we are in target volume, depending on "outer"
      CHECK(traj.getLength(1) / 1_m == Approx(0).margin(1e-3));
      CHECK(nextVol == targetPtr);
    }
    // move forward, until we leave target volume
    bool hit_neutral = false;
    bool hit_2nd = false;
    bool hit_2nd_behind = false;
    [[maybe_unused]] bool hit_2nd_partly_behind = false;
    while (nextVol != worldPtr) {
      if (nextVol == targetPtr_neutral) hit_neutral = true;
      if (nextVol == targetPtr_2) hit_2nd = true;
      if (nextVol == targetPtr_2_behind) hit_2nd_behind = true;
      if (nextVol == targetPtr_2_partly_behind) hit_2nd_partly_behind = true;
      const auto [traj2, nextVol2] = tracking.getTrack(particle);
      nextVol = nextVol2;
      particle.setNode(nextVol);
      particle.setPosition(traj2.getPosition(1));
      particle.setDirection(traj2.getDirection(1));
      CORSIKA_LOG_TRACE("pos={}, p={}, |p|={} |v|={}, delta-l={}, delta-t={}",
                        particle.getPosition(), particle.getMomentum(),
                        particle.getMomentum().getNorm(),
                        particle.getVelocity().getNorm(), traj2.getLength(1),
                        traj2.getLength(1) / particle.getVelocity().getNorm());
      CHECK(speed_0 / traj2.getVelocity(1).getNorm() == Approx(1));
    }
    CHECK_FALSE(hit_2nd_behind); // this can never happen
    // the next line is maybe an actual BUG: this should be investigated and eventually
    // fixed:
    // CHECK(hit_2nd == hit_2nd_partly_behind); // if one is hit, the other also must be
    CHECK(nextVol == worldPtr);
    CHECK(hit_2nd == true);
    CHECK(hit_neutral == (deflect == 0 ? true : false));

    Point pointCheck(cs, (deflect == 0 ? radius : 0_m), (deflect * radius), 0_m);

    CORSIKA_LOG_DEBUG(
        "testTrackingLineStack: deflect={}, momentum={}, pos={}, pos_check={}", deflect,
        particle.getMomentum().getComponents(), particle.getPosition().getCoordinates(),
        pointCheck.getCoordinates());

    CHECK((particle.getPosition() - pointCheck).getNorm() / radius ==
          Approx(0).margin(1e-1));
  }
}

/** specifc test for curved leap-frog algorithm. */

TEST_CASE("TrackingLeapFrogCurved") {

  logging::set_level(logging::level::info);

  const HEPEnergyType P0 = 10_GeV;

  corsika::Code PID = Code::MuPlus;

  using MyHomogeneousModel = media::MediumPropertyModel<
      UniformMagneticField<media::HomogeneousMedium<DummyEnvironmentInterface>>>;

  SECTION("infinite sphere / universe") {

    auto [env, csPtr, worldPtr] =
        corsika::setup::testing::setup_environment(Code::Oxygen, 100_uT);
    { [[maybe_unused]] const auto& env_dummy = env; }
    auto const& cs = *csPtr;

    tracking_leapfrog_curved::Tracking tracking;
    Point const center(cs, {0_m, 0_m, 0_m});

    auto [stack, viewPtr] = setup::testing::setup_stack(PID, P0, worldPtr, cs);
    { [[maybe_unused]] auto& viewPtr_dum = viewPtr; }
    auto particle = stack->first();
    // Note: momentum in X-direction
    //       magnetic field in X-direction

    particle.setPosition(Point(cs, 0_m, 0_m, 0_m));

    Sphere sphere(center, 1_km * std::numeric_limits<double>::infinity());

    auto intersections = tracking.intersect(particle, sphere);
    // this must be a "linear trajectory" with no curvature

    CHECK(intersections.hasIntersections() == false);
  }

  SECTION("momentum along field") {

    auto [env, csPtr, worldPtr] =
        corsika::setup::testing::setup_environment(Code::Oxygen, 100_uT);
    { [[maybe_unused]] const auto& env_dummy = env; }
    auto const& cs = *csPtr;

    tracking_leapfrog_curved::Tracking tracking;
    Point const center(cs, {0_m, 0_m, 0_m});
    auto target = DummyEnvironment::createNode<Sphere>(center, 10_km);

    MagneticFieldVector magneticfield(cs, 100_T, 0_T, 0_uT);
    target->setModelProperties<MyHomogeneousModel>(
        Medium::AirDry1Atm, magneticfield, 1_g / (1_m * 1_m * 1_m),
        media::NuclearComposition(std::vector<Code>{Code::Oxygen},
                                  std::vector<double>{1.}));
    auto* targetPtr = target.get();
    worldPtr->addChild(std::move(target));

    auto [stack, viewPtr] = setup::testing::setup_stack(PID, P0, targetPtr, cs);
    { [[maybe_unused]] auto& viewPtr_dum = viewPtr; } // prevent warning
    auto particle = stack->first();
    // Note: momentum in X-direction
    //       magnetic field in X-direction

    particle.setNode(targetPtr); // set particle outside "target" volume
    particle.setPosition(Point(cs, 0_m, 0_m, 0_m));

    auto [traj, nextVol] = tracking.getTrack(particle);
    { [[maybe_unused]] auto const& dummy = nextVol; } // prevent warning

    // this must be a "linear trajectory" with no curvature
    CHECK(traj.getDirection(0).getComponents() == traj.getDirection(1).getComponents());
  }
}

TEMPLATE_TEST_CASE("TrackingFail", "doesntwork", tracking_leapfrog_curved::Tracking,
                   tracking_leapfrog_straight::Tracking, tracking_line::Tracking) {

  logging::set_level(logging::level::info);

  const HEPEnergyType P0 = 10_GeV;

  auto [env, csPtr, worldPtr] =
      corsika::setup::testing::setup_environment(Code::Oxygen, 1_mT);
  { [[maybe_unused]] const auto& env_dummy = env; }
  auto const& cs = *csPtr;

  TestType tracking;
  Point const center(cs, {0_m, 0_m, 0_m});

  auto [stack, viewPtr] = setup::testing::setup_stack(Code::Proton, P0, worldPtr, cs);
  { [[maybe_unused]] auto& viewPtr_dum = viewPtr; }
  auto particle = stack->first();
  NonExistingDummyObject const dummy;
  CHECK_THROWS(tracking.intersect(particle, dummy));
}

TEMPLATE_TEST_CASE("TrackingSeparationPlane", "separation_plane",
                   tracking_leapfrog_curved::Tracking,
                   tracking_leapfrog_straight::Tracking, tracking_line::Tracking) {

  logging::set_level(logging::level::info);

  const HEPEnergyType P0 = 10_GeV;

  auto PID = GENERATE(as<Code>{}, Code::MuPlus, Code::MuPlus, Code::Photon);
  // for algorithms that know magnetic deflections choose: +-50uT, 0uT
  // otherwise just 0uT
  auto Bfield = GENERATE(filter(
      []([[maybe_unused]] MagneticFluxType v) {
        if constexpr (std::is_same_v<TestType, tracking_line::Tracking>)
          return v == 0_uT;
        else
          return true;
      },
      values<MagneticFluxType>({50_uT, 0_uT, -50_uT})));

  SECTION(fmt::format("Tracking PID={}, Bfield={} uT", PID, Bfield / 1_uT)) {

    CORSIKA_LOG_DEBUG(
        "********************\n                          TEST algo={} section PID={}, "
        "Bfield={}uT ",
        boost::typeindex::type_id<TestType>().pretty_name(), PID, Bfield / 1_uT);

    const int chargeNumber = get_charge_number(PID);
    LengthType radius = 10_m;
    int deflect = 0;
    if (chargeNumber != 0 and Bfield != 0_T) {
      deflect = 1;
      LengthType const gyroradius = (convert_HEP_to_SI<MassType::dimension_type>(P0) *
                                     constants::c / (abs(get_charge(PID)) * abs(Bfield)));
      CORSIKA_LOG_DEBUG("Rg={} deflect={}", gyroradius, deflect);
      radius = gyroradius;
    }

    auto [env, csPtr, worldPtr] =
        corsika::setup::testing::setup_environment(Code::Oxygen, Bfield);
    { [[maybe_unused]] const auto& env_dummy = env; }
    auto const& cs = *csPtr;

    TestType tracking;
    Point const center(cs, {0_m, 0_m, 0_m});

    auto [stack, viewPtr] = setup::testing::setup_stack(PID, P0, worldPtr, cs);
    { [[maybe_unused]] auto& viewPtr_dum = viewPtr; }
    auto particle = stack->first();
    // Note: momentum in X-direction
    //       magnetic field in Z-direction

    particle.setPosition(Point(cs, -radius, 0_m, 0_m));

    // put plane where we expect deflection by radius/2
    Plane const plane(Point(cs, radius * (1 - sqrt(3. / 4.)), 0_m, 0_m),
                      DirectionVector(cs, {-1., 0., 0.}));
    SeparationPlane const sepPlane{plane};

    Intersections const hit = tracking.intersect(particle, sepPlane);

    CORSIKA_LOG_DEBUG("entry={} exit={}", hit.getEntry(), hit.getExit());

    CHECK(hit.hasIntersections());
    CHECK(hit.getEntry() / 1_s == Approx(0.00275 * deflect).margin(0.0003));
    CHECK(hit.getExit() == 1_s * std::numeric_limits<double>::infinity());
  }
}

TEMPLATE_TEST_CASE("TrackingPlane", "plane", tracking_leapfrog_curved::Tracking,
                   tracking_leapfrog_straight::Tracking, tracking_line::Tracking) {

  logging::set_level(logging::level::info);

  const HEPEnergyType P0 = 10_GeV;

  auto PID = GENERATE(as<Code>{}, Code::MuPlus, Code::MuPlus, Code::Photon);
  // for algorithms that know magnetic deflections choose: +-50uT, 0uT
  // otherwise just 0uT
  auto Bfield = GENERATE(filter(
      []([[maybe_unused]] MagneticFluxType v) {
        if constexpr (std::is_same_v<TestType, tracking_line::Tracking>)
          return v == 0_uT;
        else
          return true;
      },
      values<MagneticFluxType>({50_uT, 0_uT, -50_uT})));

  SECTION(fmt::format("Tracking PID={}, Bfield={} uT", PID, Bfield / 1_uT)) {

    CORSIKA_LOG_DEBUG(
        "********************\n                          TEST algo={} section PID={}, "
        "Bfield={}uT ",
        boost::typeindex::type_id<TestType>().pretty_name(), PID, Bfield / 1_uT);

    const int chargeNumber = get_charge_number(PID);
    LengthType radius = 10_m;
    int deflect = 0;
    if (chargeNumber != 0 and Bfield != 0_T) {
      deflect = 1;
      LengthType const gyroradius = (convert_HEP_to_SI<MassType::dimension_type>(P0) *
                                     constants::c / (abs(get_charge(PID)) * abs(Bfield)));
      CORSIKA_LOG_DEBUG("Rg={} deflect={}", gyroradius, deflect);
      radius = gyroradius;
    }

    auto [env, csPtr, worldPtr] =
        corsika::setup::testing::setup_environment(Code::Oxygen, Bfield);
    { [[maybe_unused]] const auto& env_dummy = env; }
    auto const& cs = *csPtr;

    TestType tracking;
    Point const center(cs, {0_m, 0_m, 0_m});

    auto [stack, viewPtr] = setup::testing::setup_stack(PID, P0, worldPtr, cs);
    { [[maybe_unused]] auto& viewPtr_dum = viewPtr; }
    auto particle = stack->first();
    // Note: momentum in X-direction
    //       magnetic field in Z-direction

    particle.setPosition(Point(cs, -radius, 0_m, 0_m));

    // put plane where we expect deflection by radius/2
    Plane const plane(Point(cs, radius * (1 - sqrt(3. / 4.)), 0_m, 0_m),
                      DirectionVector(cs, {-1., 0., 0.}));
    Intersections const hit = tracking.intersect(particle, plane);

    CORSIKA_LOG_DEBUG("entry={} exit={}", hit.getEntry(), hit.getExit());

    CHECK(hit.hasIntersections());
    CHECK(hit.getEntry() / 1_s == Approx(0.00275 * deflect).margin(0.0003));
    CHECK(hit.getExit() == 1_s * std::numeric_limits<double>::infinity());
  }
}

TEST_CASE("Intersections") {
  Intersections test;
  CHECK(test.getEntry() ==
        std::numeric_limits<TimeType::value_type>::infinity() * second);
  CHECK(test.getExit() == std::numeric_limits<TimeType::value_type>::infinity() * second);
}
