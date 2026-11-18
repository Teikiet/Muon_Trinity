/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/output/OutputManager.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/refractivity/UniformRefractiveIndex.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <corsika/modules/radio/RadioProcess.hpp>
#include <corsika/modules/radio/CoREAS.hpp>
#include <corsika/modules/radio/ZHS.hpp>
#include <corsika/modules/radio/observers/TimeDomainObserver.hpp>
#include <corsika/modules/radio/detectors/ObserverCollection.hpp>
#include <corsika/modules/radio/propagators/DummyTestPropagator.hpp>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace corsika;
using namespace std;

//
// A simple example to get the electric field trace of an electron using manual tracking
//
int main() {

  // create a suitable environment
  using IModelInterface = media::IRefractiveIndexModel<
      media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>>;
  using AtmModel = media::UniformRefractiveIndex<media::MediumPropertyModel<
      media::UniformMagneticField<media::HomogeneousMedium<IModelInterface>>>>;
  using EnvType = media::Environment<AtmModel>;
  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  // a refractive index for the vacuum
  const double ri_{1};
  // the constant density
  const auto density{19.2_g / cube(1_cm)};
  // the composition we use for the homogeneous medium
  media::NuclearComposition const Composition({Code::Nitrogen}, {1.});
  // create magnetic field vector
  Vector B1(rootCS, 0_T, 0_T, 0.3809_T);
  // create a Sphere for the medium
  auto Medium =
      EnvType::createNode<Sphere>(center, 1_km * std::numeric_limits<double>::infinity());
  // set the environment properties
  auto const props = Medium->setModelProperties<AtmModel>(ri_, media::Medium::AirDry1Atm,
                                                          B1, density, Composition);
  // bind things together
  env.getUniverse()->addChild(std::move(Medium));
  auto const& node_ = env.getUniverse()->getChildNodes().front();

  // the observers location
  const auto point1{Point(rootCS, 30000_m, 0_m, 0_m)};

  // create times for the observer
  const TimeType start{0.994e-4_s};
  const TimeType duration{1.07e-4_s - 0.994e-4_s};
  // 3 km observer
  //    const TimeType start{0.994e-5_s};
  //    const TimeType duration{1.7e-5_s - 0.994e-5_s};
  const InverseTimeType sampleRate_{5e+11_Hz};

  std::cout << "number of points in time: " << duration * sampleRate_ << std::endl;

  // create 2 observers
  TimeDomainObserver obs1("CoREAS observer", point1, rootCS, start, duration, sampleRate_,
                          start);
  TimeDomainObserver obs2("ZHS observer", point1, rootCS, start, duration, sampleRate_,
                          start);

  // construct a radio detector instance to store our observers
  ObserverCollection<TimeDomainObserver> detectorCoREAS;
  ObserverCollection<TimeDomainObserver> detectorZHS;

  // add the observers to the detector
  detectorCoREAS.addObserver(obs1);
  detectorZHS.addObserver(obs2);

  // create a new stack for each trial
  setup::Stack<EnvType> stack;
  stack.clear();

  const Code particle{Code::Electron};
  const HEPMassType pmass{get_mass(particle)};

  // construct an energy // move in the for loop
  const HEPEnergyType E0{11.4_MeV};

  // construct the output manager
  OutputManager outputs("synchrotron_radiation_manual_tracking-output");

  // the radio signal propagator
  auto SP = make_dummy_test_radio_propagator(env);

  // create a radio process instance using CoREAS
  RadioProcess<ObserverCollection<TimeDomainObserver>,
               CoREAS<ObserverCollection<TimeDomainObserver>, decltype(SP)>, decltype(SP)>
      coreas(detectorCoREAS, SP);
  // register CoREAS to the output manager
  outputs.add("CoREAS", coreas);

  // create a radio process instance using ZHS
  RadioProcess<ObserverCollection<TimeDomainObserver>,
               ZHS<ObserverCollection<TimeDomainObserver>, decltype(SP)>, decltype(SP)>
      zhs(detectorZHS, SP);
  // register ZHS to the output manager
  outputs.add("ZHS", zhs);

  // trigger the start of the library and the first event
  outputs.startOfLibrary();
  outputs.startOfShower();

  // the number of points that make up the circle
  int const n_points{100000};
  LengthType const radius{100_m};
  TimeType timeCounter{0._s};

  // loop over the circular trajectory (this produces 1 pulse)
  for (size_t i = 0; i <= (n_points); i++) {
    Point const point_1(rootCS, {radius * cos(M_PI * 2 * i / n_points),
                                 radius * sin(M_PI * 2 * i / n_points), 0_m});
    Point const point_2(rootCS, {radius * cos(M_PI * 2 * (i + 1) / n_points),
                                 radius * sin(M_PI * 2 * (i + 1) / n_points), 0_m});
    TimeType t{(point_2 - point_1).getNorm() / (0.999 * constants::c)};
    timeCounter = timeCounter + t;
    VelocityVector v{(point_2 - point_1) / t};
    auto beta{v / constants::c};
    auto gamma{E0 / pmass};
    auto plab{beta * pmass * gamma};
    Line l{point_1, v};
    StraightTrajectory track{l, t};
    auto particle1{stack.addParticle(std::make_tuple(
        particle, calculate_kinetic_energy(plab.getNorm(), get_mass(particle)),
        plab.normalized(), point_1, timeCounter))};
    particle1.setNode(node_.get());
    Step step(particle1, track);
    coreas.doContinuous(step, true);
    zhs.doContinuous(step, true);
    stack.clear();
  }

  // trigger the manager to write the data to disk
  outputs.endOfShower();
  outputs.endOfLibrary();
}
