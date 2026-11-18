/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <array>
#include <iomanip>
#include <iostream>

#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Vector.hpp>

using namespace corsika;
using namespace std;

const int nSteps = 10;

class Process1 : public ContinuousProcess<Process1> {
public:
  Process1() {}
  template <typename D>
  ProcessReturn doContinuous(Step<D>& d, bool const) const {
    LengthVector displacement_{get_root_CoordinateSystem(), 1_m, 0_m, 0_m};
    DirectionVector dU_{get_root_CoordinateSystem(), {1, 0, 0}};
    for (int i = 0; i < nSteps; ++i) {
      d.add_displacement(displacement_);
      d.add_dU(dU_);
    }
    return ProcessReturn::Ok;
  }
};

class Process2 : public ContinuousProcess<Process2> {
public:
  Process2() {}

  template <typename D>
  inline ProcessReturn doContinuous(Step<D>& d, bool const) const {
    LengthVector displacement_{get_root_CoordinateSystem(), 0.1_m, 0_m, 0_m};
    for (int i = 0; i < nSteps; ++i) { d.add_displacement(-displacement_); }
    return ProcessReturn::Ok;
  }
};

class Process3 : public ContinuousProcess<Process3> {
public:
  Process3() {}

  template <typename D>
  inline ProcessReturn doContinuous(Step<D>&, bool const) const {
    return ProcessReturn::Ok;
  }
};

class Process4 : public ContinuousProcess<Process4> {
public:
  Process4(const double v)
      : fV(v) {}
  template <typename D>
  inline ProcessReturn doContinuous(Step<D>& d, bool const) const {
    for (int i = 0; i < nSteps; ++i) { d.add_dEkin(i * fV * 1_eV); }
    return ProcessReturn::Ok;
  }

private:
  double fV;
};

struct DummyData {
  Point getPosition() const { return Point(get_root_CoordinateSystem(), 0_m, 0_m, 0_m); }
  DirectionVector getDirection() const {
    return DirectionVector{get_root_CoordinateSystem(), {0, 0, 0}};
  }
};
struct DummyStack {
  void clear() {}
};
struct DummyTrajectory {
  TimeType getDuration(int u) const { return 0_s; }
  Point getPosition(int u) const {
    return Point(get_root_CoordinateSystem(), 0_m, 0_m, 0_m);
  }
  DirectionVector getDirection(int u) const {
    return DirectionVector{get_root_CoordinateSystem(), {0, 0, 0}};
  }
};

void modular() {

  //              = 0
  Process1 m1;      // + 1.0m
  Process2 m2;      // - 0.1m
  Process3 m3;      // * 1.0
  Process4 m4(1.5); // * i * 1.5_eV

  auto sequence = make_sequence(m1, m2, m3, m4);

  DummyData particle;
  DummyTrajectory track;
  Step step(particle, track);

  sequence.doContinuous(step, ContinuousProcessIndex(0));
  double X = 0;
  double energy = 0;
  double direction = 0;
  for (int i = 0; i < nSteps; ++i) {
    X += 1. - 0.1;
    energy += i * 1.5;
    direction += 1.;
  }

  assert(static_cast<int>(step.getDisplacement().getX(get_root_CoordinateSystem()) /
                          1_m) == static_cast<int>(X));
  assert((step.getDiffEkin() / 1_eV) == energy);
  assert((step.getDiffDirection().getX(get_root_CoordinateSystem())) == direction);

  cout << " done (checking...) " << endl;
}

int main() {

  logging::set_level(logging::level::info);

  std::cout << "staticsequence_example" << std::endl;

  modular();
  return 0;
}
