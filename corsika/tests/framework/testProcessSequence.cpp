/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#define CORSIKA_UNIT_TESTING

#include <corsika/framework/process/ContinuousProcessStepLength.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/ProcessTraits.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Step.hpp>

#include <corsika/framework/utility/COMBoost.hpp>

#include <corsika/media/composition/NuclearComposition.hpp>

#include <catch2/catch_all.hpp>

#include <array>
#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <boost/type_index.hpp>

/*
  Unit test for testing all Process types and their arrangement in
  containers ProcessSequence and SwitchProcessSequence
 */

using namespace corsika;
using namespace media;
using namespace std;
using Catch::Approx;

struct DummyRNG {
  static constexpr int max() { return 10; }
  static constexpr int min() { return 0; }
  double operator()() const { return 0.5; }
};

static int const nData = 10;

// DummyNode is only needed for BoundaryCrossingProcess
struct DummyNode {
  DummyNode(int v)
      : data_(v) {}
  int data_ = 0;
};

// our data object (particle) is a simple arrary of doubles
struct DummyData {
  double data_[nData] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  typedef DummyNode node_type; // for BoundaryCrossingProcess
  Code getPID() const { return Code::Proton; }
  MomentumVector getMomentum() const {
    // only need the coordinate system
    return MomentumVector{0_eV, 0_eV, 0_eV};
  }
  HEPEnergyType getEnergy() const { return 10_GeV; }
  Point getPosition() const { return Point(0_m, 0_m, 0_m); }
  DirectionVector getDirection() const { return DirectionVector{{0, 0, 0}}; }
};

// The stack is non-existent for this example
struct DummyStack {};

// there is no real trajectory/track
struct DummyTrajectory {
  TimeType getDuration([[maybe_unused]] int u) const { return 0_s; }
  Point getPosition([[maybe_unused]] int u) const { return Point(0_m, 0_m, 0_m); }
  DirectionVector getDirection([[maybe_unused]] int u) const {
    return DirectionVector{{0, 0, 0}};
  }
};

// since there is no stack, there is also no view. This is a simplistic dummy
// object sufficient here.
struct DummyView {
  DummyView(DummyData& p)
      : p_(p) {}
  DummyData& p_;
  DummyData& parent() { return p_; }
  // this is only needed because of PROPOSAL interface right now:
};

int globalCount = 0; // simple counter

int checkDecay = 0;       // use this as a bit field
int checkInteract = 0;    // use this as a bit field
int checkSec = 0;         // use this as a bit field
int checkCont = 0;        // use this as a bit field
int checkSecondaries = 0; // use this as a bit field

class ContinuousProcess1 : public ContinuousProcess<ContinuousProcess1> {
public:
  ContinuousProcess1(int const v, LengthType const step)
      : v_(v)
      , step_(step) {

    CORSIKA_LOG_DEBUG(
        "globalCount: {} "
        ", v_: {} ",
        globalCount, v_);
    globalCount++;
  }

  void setStep(LengthType const v) { step_ = v; }

  template <typename D>
  ProcessReturn doContinuous(Step<D>& d, bool flag) const {
    flag_ = flag;
    CORSIKA_LOG_TRACE("ContinuousProcess1::DoContinuous");
    checkCont |= 1;
    LengthVector displacement_{get_root_CoordinateSystem(), 1_m, 0_m, 0_m};
    DirectionVector dU_{get_root_CoordinateSystem(), {1, 0, 0}};
    d.add_dt(1_s);
    d.add_displacement(displacement_);
    d.add_dU(dU_);
    d.add_dEkin(1_eV);
    return ProcessReturn::Ok;
  }

  template <typename TParticle, typename TTrack>
  LengthType getMaxStepLength(TParticle&, TTrack&) {
    return step_;
  }

  bool getFlag() const { return flag_; }
  void resetFlag() { flag_ = false; }

private:
  int v_ = 0;
  LengthType step_ = 0_m;
  mutable bool flag_ = false;
};

class ContinuousProcess2 : public ContinuousProcess<ContinuousProcess2> {
public:
  ContinuousProcess2([[maybe_unused]] int const v, LengthType const step)
      : step_(step) {
    CORSIKA_LOG_DEBUG("globalCount: {}", globalCount);
    globalCount++;
  }

  void setStep(LengthType const v) { step_ = v; }

  template <typename D>
  ProcessReturn doContinuous(Step<D>& d, bool const flag) const {
    flag_ = flag;
    CORSIKA_LOG_DEBUG("ContinuousProcess2::DoContinuous");
    checkCont |= 2;
    d.add_dt(10_s);
    return ProcessReturn::Ok;
  }

  template <typename TParticle, typename TTrack>
  LengthType getMaxStepLength(TParticle&, TTrack&) {
    return step_;
  }

  bool getFlag() const { return flag_; }
  void resetFlag() { flag_ = false; }

private:
  LengthType step_ = 0_m;
  mutable bool flag_ = false;
};

class ContinuousProcess3 : public ContinuousProcess<ContinuousProcess3> {
public:
  ContinuousProcess3([[maybe_unused]] int const v, LengthType const step)
      : step_(step) {
    CORSIKA_LOG_DEBUG("globalCount: {}", globalCount);
    globalCount++;
  }

  void setStep(LengthType const v) { step_ = v; }

  template <typename D>
  ProcessReturn doContinuous(Step<D>& d, bool const flag) const {
    flag_ = flag;
    CORSIKA_LOG_DEBUG("ContinuousProcess3::DoContinuous");
    checkCont |= 4;
    for (int i = 0; i < nData; ++i) d.add_dEkin(1_eV);
    return ProcessReturn::Ok;
  }

  template <typename TParticle, typename TTrack>
  LengthType getMaxStepLength(TParticle&, TTrack&) {
    return step_;
  }

  bool getFlag() const { return flag_; }
  void resetFlag() { flag_ = false; }

private:
  LengthType step_ = 0_m;
  mutable bool flag_ = false;
};

class Process1 : public InteractionProcess<Process1> {
public:
  Process1(int const v)
      : v_(v) {
    CORSIKA_LOG_DEBUG(
        "globalCount: {}"
        ", v_: {}",
        globalCount, v_);
    ;
    globalCount++;
  }

  template <typename TView>
  void doInteraction(TView& v, Code const, Code const, FourMomentum const&,
                     FourMomentum const&) const {
    checkInteract |= 1;
    for (int i = 0; i < nData; ++i) v.parent().data_[i] += 1 + i;
  }

  CrossSectionType getCrossSection(Code const, Code const, FourMomentum const&,
                                   FourMomentum const&) const {
    return 10_mb;
  }

private:
  int v_;
};

class Process2 : public InteractionProcess<Process2> {
public:
  Process2(int const v)
      : v_(v) {
    CORSIKA_LOG_DEBUG(
        "globalCount: {}"
        ", v_: {}",
        globalCount, v_);
    globalCount++;
  }

  template <typename TView>
  void doInteraction(TView& v, Code const, Code const, FourMomentum const&,
                     FourMomentum const&) const {
    checkInteract |= 2;
    for (int i = 0; i < nData; ++i) v.parent().data_[i] /= 1.1;
    CORSIKA_LOG_DEBUG("Process2::doInteraction");
  }

  CrossSectionType getCrossSection(Code const, Code const, FourMomentum const&,
                                   FourMomentum const&) const {
    CORSIKA_LOG_DEBUG("Process2::getCrossSection");
    return 20_mb;
  }

private:
  int v_ = 0;
};

class Process3 : public InteractionProcess<Process3> {
public:
  Process3(int const v)
      : v_(v) {
    CORSIKA_LOG_DEBUG(
        "globalCount: {}"
        ", v_: {}",
        globalCount, v_);
    globalCount++;
  }

  template <typename TView>
  void doInteraction(TView& v, Code const, Code const, FourMomentum const&,
                     FourMomentum const&) const {
    checkInteract |= 4;
    for (int i = 0; i < nData; ++i) v.parent().data_[i] *= 1.01;
    CORSIKA_LOG_DEBUG("Process3::doInteraction");
  }

  CrossSectionType getCrossSection(Code const, Code const, FourMomentum const&,
                                   FourMomentum const&) const {
    CORSIKA_LOG_DEBUG("Process3::getCrossSection");
    return 30_mb;
  }

private:
  int v_ = 0;
};

class Process4 : public BaseProcess<Process4> {
public:
  Process4(int const v)
      : v_(v) {
    CORSIKA_LOG_DEBUG(
        "globalCount: {}"
        ", v_: {}",
        globalCount, v_);
    globalCount++;
  }

  template <typename D>
  ProcessReturn doContinuous(Step<D>& d, bool const) const {
    CORSIKA_LOG_DEBUG("Base::doContinuous");
    checkCont |= 8;
    for (int i = 0; i < nData; ++i) { d.add_dEkin(1_eV); }
    return ProcessReturn::Ok;
  }
  template <typename TView>
  void doInteraction(TView&, Code const, Code const, FourMomentum const&,
                     FourMomentum const&) const {
    checkInteract |= 8;
  }

private:
  int v_ = 0;
};

class Decay1 : public DecayProcess<Decay1> {
public:
  Decay1(int const) {
    CORSIKA_LOG_DEBUG("Decay1()");
    globalCount++;
  }

  template <typename Particle>
  TimeType getLifetime(Particle&) const {
    return 1_s;
  }
  template <typename TView>
  void doDecay(TView&) const {
    checkDecay |= 1;
  }
};

class Decay2 : public DecayProcess<Decay2> {
public:
  Decay2(int const) {
    CORSIKA_LOG_DEBUG("Decay2()");
    globalCount++;
  }

  template <typename Particle>
  TimeType getLifetime(Particle&) const {
    return 2_s;
  }
  void doDecay(DummyView&) const { checkDecay |= 2; }
};

class Stack1 : public StackProcess<Stack1> {
public:
  Stack1(int const n)
      : StackProcess(n) {}
  template <typename TStack>
  void doStack(TStack const&) {
    count_++;
  }
  int getCount() const { return count_; }

private:
  int count_ = 0;
};

class Secondaries1 : public SecondariesProcess<Secondaries1> {
public:
  template <typename TView>
  void doSecondaries(TView const&) {
    checkSecondaries |= 1;
  }
};

class Secondaries2 : public SecondariesProcess<Secondaries2> {
public:
  template <typename TView>
  void doSecondaries(TView const&) {
    checkSecondaries |= 2;
  }
};

class Boundary1 : public BoundaryCrossingProcess<Boundary1> {
public:
  Boundary1(double const v = 1.0)
      : v_(v) {}

  template <typename Particle>
  ProcessReturn doBoundaryCrossing(Particle& p, typename Particle::node_type const& from,
                                   typename Particle::node_type const& to) {

    for (int i = 0; i < nData; ++i) { p.data_[i] += v_ * (from.data_ - to.data_); }
    return ProcessReturn::Ok;
  }

private:
  double v_ = 0.0;
};

TEST_CASE("ProcessSequence General", "ProcessSequence") {

  logging::set_level(logging::level::info);

  SECTION("BaseProcess") {

    Process1 m1(0);
    const Process4 m4(3);

    CHECK(is_process_v<Process1>);
    CHECK(is_process_v<decltype(m4)>);
    CHECK(is_process_v<decltype(Decay1(1))>);
    CHECK(is_process_v<decltype(ContinuousProcess3{3, 3_m})>);
  }

  SECTION("Check construction") {
    globalCount = 0;
    Process1 m1(0);
    CHECK(globalCount == 1);
    Process2 m2(1);
    CHECK(globalCount == 2);
    Process3 m3(2);
    CHECK(globalCount == 3);
    Process4 m4(3);
    CHECK(globalCount == 4);

    auto sequence1 = make_sequence(m1, m2, m3, m4);
    CHECK(is_process_v<decltype(sequence1)>);
    CHECK(is_process_v<decltype(m2)>);
    CHECK(decltype(sequence1)::is_process_sequence);
    CHECK_FALSE(decltype(m2)::is_process_sequence);
    CHECK_FALSE(decltype(sequence1)::is_switch_process_sequence);
    CHECK_FALSE(decltype(m2)::is_switch_process_sequence);

    CHECK_FALSE(decltype(Decay1(7))::is_process_sequence);
    CHECK_FALSE(decltype(Decay1(7))::is_switch_process_sequence);

    auto sequence2 = make_sequence(m1, m2, m3);
    CHECK(decltype(sequence2)::is_process_sequence == true);

    auto sequence3 = make_sequence(m4);
    CHECK(decltype(sequence3)::is_process_sequence == true);

    CHECK(std::is_reference_v<decltype(sequence3.getProcess1())>);  // Process4&
    CHECK(!std::is_reference_v<decltype(sequence3.getProcess2())>); // NullModel

    CHECK(std::is_reference_v<decltype(sequence2.getProcess1())>);  // Process1&
    CHECK(!std::is_reference_v<decltype(sequence2.getProcess2())>); // ProcessSequence
    CHECK(std::is_reference_v<decltype(
              sequence2.getProcess2().getProcess1())>); // Process2&
    CHECK(std::is_reference_v<decltype(
              sequence2.getProcess2().getProcess2())>); // Process3&

    // and now with rvalue initialization

    auto sequence2_rv = make_sequence(Process1(0), m2, Process3(0));
    CHECK(!std::is_reference_v<decltype(sequence2_rv.getProcess1())>); // Process1
    CHECK(!std::is_reference_v<decltype(sequence2_rv.getProcess2())>); // ProcessSequence
    CHECK(std::is_reference_v<decltype(
              sequence2_rv.getProcess2().getProcess1())>); // Process2&
    CHECK(!std::is_reference_v<decltype(
              sequence2_rv.getProcess2().getProcess2())>); // Process3
  }

  SECTION("lifetime") {
    globalCount = 0;
    ContinuousProcess1 cp1(0, 1_m);
    Process2 m2(1);
    Process3 m3(2);
    Decay1 d3(3);

    DummyData particle;

    auto sequence2 = make_sequence(cp1, m2, m3, d3);
    TimeType const tot = sequence2.getLifetime(particle);
    InverseTimeType const tot_inv = sequence2.getInverseLifetime(particle);
    CORSIKA_LOG_DEBUG(
        "lambda_tot={}"
        "; lambda_tot_inv={}",
        tot, tot_inv);

    CHECK(tot / 1_s == 1);
    CHECK(tot_inv * 1_s == 1.);
    globalCount = 0;
  }

  SECTION("ContinuousProcess") {
    globalCount = 0;
    ContinuousProcess1 cp1(0, 1_m);   // += 0.933
    ContinuousProcess2 cp2(1, 1.1_m); // += 0.111
    Process2 m2(2);                   //  /= 1.1
    Process3 m3(3);                   //  *= 1.01

    auto sequence2 = make_sequence(cp1, m2, m3, cp2);

    std::cout << boost::typeindex::type_id<decltype(sequence2)>().pretty_name()
              << std::endl;

    DummyData particle;
    DummyTrajectory track;
    Step step(particle, track);

    cp1.resetFlag();
    cp2.resetFlag();

    ContinuousProcessStepLength const step1 = sequence2.getMaxStepLength(particle, track);
    CHECK(LengthType(step1) == 1_m);
    sequence2.doContinuous(step, step1);
    CHECK(cp1.getFlag());
    CHECK_FALSE(cp2.getFlag());
    CORSIKA_LOG_INFO("step1, l={}, i={}", LengthType(step1),
                     ContinuousProcessIndex(step1).getIndex());

    cp1.resetFlag();
    cp2.resetFlag();

    cp1.setStep(10_m);
    ContinuousProcessStepLength const step2 = sequence2.getMaxStepLength(particle, track);
    CHECK(LengthType(step2) == 1.1_m);
    CHECK(ContinuousProcessIndex(step1) != ContinuousProcessIndex(step2));
    sequence2.doContinuous(step, step2);
    CHECK_FALSE(cp1.getFlag());
    CHECK(cp2.getFlag());
    CORSIKA_LOG_INFO("step2, l={}, i={}", LengthType(step2),
                     ContinuousProcessIndex(step2).getIndex());

    CORSIKA_LOG_DEBUG("-->init sequence2");
    globalCount = 0;
    CORSIKA_LOG_DEBUG("-->docont");

    // reset
    particle = DummyData();
    track = DummyTrajectory();

    int const nLoop = 5;
    CORSIKA_LOG_DEBUG("Running loop with n={}", nLoop);
    for (int iLoop = 0; iLoop < nLoop; ++iLoop) {
      sequence2.doContinuous(step, ContinuousProcessIndex(&cp1));
    }
    CHECK(step.getDiffT() / 1_s == Approx(77));
    CHECK(step.getDiffEkin() / 1_eV == Approx(7));
    CHECK(step.getDisplacement().getX(get_root_CoordinateSystem()) / 1_m == Approx(7));
    CHECK(step.getDisplacement().getY(get_root_CoordinateSystem()) / 1_m == Approx(0));
    CHECK(step.getDisplacement().getZ(get_root_CoordinateSystem()) / 1_m == Approx(0));
    CHECK(step.getDiffDirection().getX(get_root_CoordinateSystem()) == Approx(7));
    CHECK(step.getDiffDirection().getY(get_root_CoordinateSystem()) == Approx(0));
    CHECK(step.getDiffDirection().getZ(get_root_CoordinateSystem()) == Approx(0));
    CORSIKA_LOG_DEBUG("done");
  }

  SECTION("StackProcess") {

    globalCount = 0;
    Stack1 s1(1);
    Stack1 s2(2);

    auto sequence1 = make_sequence(s1, s2);

    std::cout << boost::typeindex::type_id<decltype(sequence1)>().pretty_name()
              << std::endl;

    DummyStack stack;

    int const nLoop = 20;
    for (int i = 0; i < nLoop; ++i) { sequence1.doStack(stack); }

    CHECK(s1.getCount() == 20);
    CHECK(s2.getCount() == 10);

    ContinuousProcess2 cp2(1, 2_m); // += 0.111
    Process2 m2(2);                 //  /= 1.1
    auto sequence2 = make_sequence(cp2, m2);
    auto sequence3 = make_sequence(cp2, m2, s1);

    CHECK(decltype(sequence2)::is_process_sequence == true);
    CHECK(decltype(sequence3)::is_process_sequence == true);
    CHECK(contains_stack_process_v<decltype(sequence2)> == false);
    CHECK(contains_stack_process_v<decltype(sequence3)> == true);
  }

  SECTION("BoundaryCrossingProcess") {

    globalCount = 0;
    Boundary1 b1;

    auto sequence1 = make_sequence(b1);

    DummyData particle;
    DummyNode node_from(5);
    DummyNode node_to(4);

    int const nLoop = 20;
    for (int i = 0; i < nLoop; ++i) {
      sequence1.doBoundaryCrossing(particle, node_from, node_to);
    }

    for (int i = 0; i < nData; i++) {
      CORSIKA_LOG_DEBUG("data_[{}]={}", i, particle.data_[i]);
      CHECK(particle.data_[i] == Approx(nLoop).margin(1e-9));
    }

    CHECK(decltype(sequence1)::is_process_sequence == true);
    CHECK(contains_stack_process_v<decltype(sequence1)> == false);
    CHECK(count_processes<decltype(sequence1)>::count == 1);
  }
}

TEST_CASE("SwitchProcessSequence", "ProcessSequence") {

  logging::set_level(logging::level::info);

  /**
   * In this example switching is done only by "data_[0]>0", where
   * data in an arrray of doubles, DummyData.
   */

  struct SwitchSelect {
    bool operator()(DummyData const& p) const { return (p.data_[0] > 0); }
  };
  SwitchSelect select1;

  auto cp1 = ContinuousProcess1(0, 1_m);
  auto cp2 = ContinuousProcess2(0, 2_m);
  auto cp3 = ContinuousProcess3(0, 3_m);

  auto sec1 = Secondaries1();
  auto sec2 = Secondaries2();

  auto sequence1 =
      make_sequence(Process1(0), cp2, Decay1(0), sec1, Boundary1(1.0)); // 10 mb
  auto sequence2 = make_sequence(cp3, Process2(0), Boundary1(-1.0), Decay2(0),
                                 sec2); // 20 mb

  auto sequence3 = make_sequence(cp1, Process3(0), // 30 mb
                                 SwitchProcessSequence(select1, sequence1, sequence2));

  // it is even more typical to have just one sub-process inside the branches of
  // SwitchProcessSequence
  auto sequence3_short =
      make_sequence(cp1, Process3(0), // 30 mb
                    SwitchProcessSequence(select1, Process1(0), Process2(0)));

  auto sequence4 =
      make_sequence(cp1, Boundary1(2.0), Process3(0),
                    SwitchProcessSequence(select1, sequence1, Boundary1(-1.0)));

  SECTION("Check construction") {

    auto switch_seq = make_select(select1, sequence1, sequence2);
    CHECK(decltype(switch_seq)::is_process_sequence);
    CHECK(decltype(switch_seq)::is_switch_process_sequence);
    CHECK(decltype(SwitchProcessSequence(select1, sequence1,
                                         sequence2))::is_switch_process_sequence);

    CHECK(decltype(sequence3)::is_process_sequence);
    CHECK_FALSE(decltype(sequence3)::is_switch_process_sequence);

    auto sps1 = SwitchProcessSequence(select1, sequence1, sequence2);
    CHECK(decltype(sps1)::is_process_sequence);
    CHECK(decltype(sps1)::is_switch_process_sequence);

    std::cout << boost::typeindex::type_id<decltype(sequence3)>().pretty_name()
              << std::endl;

    CHECK(decltype(sequence3)::is_process_sequence);
    auto sps2 = SwitchProcessSequence(select1, sequence1, sequence2);
    CHECK(decltype(sps2)::is_process_sequence);

    CHECK(std::is_reference_v<decltype(switch_seq.getCondition())>);   //
    CHECK(std::is_reference_v<decltype(switch_seq.getSequence())>);    //
    CHECK(std::is_reference_v<decltype(switch_seq.getAltSequence())>); //

    // check with rvalue init
    auto switch_seq_rv =
        make_select(SwitchSelect(), make_sequence(Process1(0)), Process3(0));
    CHECK(!std::is_reference_v<decltype(switch_seq_rv.getCondition())>);
    CHECK(!std::is_reference_v<decltype(switch_seq_rv.getSequence())>);
    CHECK(!std::is_reference_v<decltype(switch_seq_rv.getAltSequence())>);
  }

  SECTION("Check interfaces") {

    DummyData particle;
    DummyTrajectory track;
    DummyView view(particle);
    Step step(particle, track);

    checkDecay = 0;
    checkInteract = 0;
    checkSec = 0;
    checkCont = 0;
    particle.data_[0] = 100; // data positive --> sequence1
    sequence3.doContinuous(step, ContinuousProcessIndex(&cp1));
    CHECK(checkInteract == 0);
    CHECK(checkDecay == 0);
    CHECK(checkCont == 0b011);
    CHECK(checkSec == 0);

    checkDecay = 0;
    checkInteract = 0;
    checkSec = 0;
    checkCont = 0;
    particle.data_[0] = -100; // data negative  --> sequence2
    sequence3.doContinuous(step, ContinuousProcessIndex(&cp1));
    CHECK(checkInteract == 0);
    CHECK(checkDecay == 0);
    CHECK(checkCont == 0b101);
    CHECK(checkSec == 0);

    // 30_mb is Process3
    CrossSectionType cx_select = .9 * 30_mb;
    InverseTimeType time_select = 0.1 / second; // for decay

    checkDecay = 0;
    checkInteract = 0;
    checkSec = 0;
    checkCont = 0;
    particle.data_[0] = 100; // data positive   --> sequence1

    DummyRNG rng;
    FourMomentum const projectileP4{10_GeV, {0_eV, 0_eV, 0_eV}};
    media::NuclearComposition const noComposition({Code::Nitrogen}, {1});
    sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
    sequence3.selectDecay(view, time_select);
    CHECK(checkInteract == 0b100); // this is Process3
    CHECK(checkDecay == 0b001);    // this is Decay1
    CHECK(checkCont == 0);
    CHECK(checkSec == 0);
    cx_select = 1.01 * 30_mb;
    checkInteract = 0;
    sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
    CHECK(checkInteract == 0b001); // this is Process1

    checkDecay = 0;
    checkInteract = 0;
    checkSec = 0;
    checkCont = 0;
    particle.data_[0] = -100; // data negative   --> sequence2
    sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
    sequence3.selectDecay(view, time_select);
    CHECK(checkInteract == 0b010); // this is Process2
    CHECK(checkDecay == 0b010);    // this is Decay2
    CHECK(checkCont == 0);
    CHECK(checkSec == 0);

    checkDecay = 0;
    checkInteract = 0;
    checkSec = 0;
    checkCont = 0;
    particle.data_[0] = -100; // data negative  --> sequence2
    sequence3.doSecondaries(view);
    Stack1 stack(0);
    sequence3.doStack(stack);
    CHECK(checkInteract == 0);
    CHECK(checkDecay == 0);
    CHECK(checkCont == 0);
    CHECK(checkSec == 0);

    // check the SwitchProcessSequence where no process is selected in
    // selected branch (fallthrough)

    checkDecay = 0;
    checkInteract = 0;
    checkSec = 0;
    checkCont = 0;
    particle.data_[0] = -100; // data negative --> sequence1
    sequence4.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
    sequence4.doSecondaries(view);
    sequence4.selectDecay(view, time_select);
    sequence4.doSecondaries(view);
    CHECK(checkInteract == 0);
    CHECK(checkDecay == 0);
    CHECK(checkCont == 0);
    CHECK(checkSec == 0);

    // now check sequence3, which contains a SwitchProcessSequence that contains
    // two longer sequences in each branch.
    {
      // check that large "select" value will correctly ignore the call
      cx_select = 1e5_mb;
      time_select = 1e5 / second;
      checkDecay = 0;
      checkInteract = 0;
      sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
      sequence3.selectDecay(view, time_select);
      CHECK(checkInteract == 0);
      CHECK(checkDecay == 0);

      // for a small cx_select selection must be successful
      cx_select = 28_mb; // -> Process3
      checkInteract = 0;
      particle.data_[0] = -100; // data negative --> sequence2
      CHECK(sequence3.getCrossSection(particle, Code::Oxygen,
                                      {Oxygen::mass, {0_eV, 0_eV, 0_eV}}) /
                1_mb ==
            Approx(50.));
      sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
      CHECK(checkInteract == 4); // 2^3

      particle.data_[0] = 100; // data positive --> sequence1
      checkInteract = 0;
      CHECK(sequence3.getCrossSection(particle, Code::Oxygen,
                                      {Oxygen::mass, {0_eV, 0_eV, 0_eV}}) /
                1_mb ==
            Approx(40.));
      sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
      CHECK(checkInteract == 4); // 2^3

      cx_select = 32_mb; // -> Process2 or Process1
      checkInteract = 0;
      particle.data_[0] = -100; // data negative --> Process2
      sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
      CHECK(checkInteract == 2); // 2^2

      particle.data_[0] = 100; // data positive --> Process1
      checkInteract = 0;
      sequence3.selectInteraction(view, projectileP4, noComposition, rng, cx_select);
      CHECK(checkInteract == 1); // 2^1
    }

    // now check sequence3, which contains a SwitchProcessSequence that contains
    // just two bare InteractionProcess-es in each branch.
    {
      // check that large "select" value will correctly ignore the call
      cx_select = 1e5_mb;
      checkInteract = 0;
      sequence3_short.selectInteraction(view, projectileP4, noComposition, rng,
                                        cx_select);
      CHECK(checkInteract == 0);

      // for a small cx_select selection must be sucessful
      cx_select = 28_mb; // -> Process3
      checkInteract = 0;
      particle.data_[0] = -100; // data negative --> sequence2
      CHECK(sequence3_short.getCrossSection(particle, Code::Oxygen,
                                            {Oxygen::mass, {0_eV, 0_eV, 0_eV}}) /
                1_mb ==
            Approx(50.));
      sequence3_short.selectInteraction(view, projectileP4, noComposition, rng,
                                        cx_select);
      CHECK(checkInteract == 4); // 2^3

      particle.data_[0] = 100; // data positive --> sequence1
      checkInteract = 0;
      CHECK(sequence3_short.getCrossSection(particle, Code::Oxygen,
                                            {Oxygen::mass, {0_eV, 0_eV, 0_eV}}) /
                1_mb ==
            Approx(40.));
      sequence3_short.selectInteraction(view, projectileP4, noComposition, rng,
                                        cx_select);
      CHECK(checkInteract == 4); // 2^3

      cx_select = 32_mb; // -> Process2 or Process1
      checkInteract = 0;
      particle.data_[0] = -100; // data negative --> Process2
      sequence3_short.selectInteraction(view, projectileP4, noComposition, rng,
                                        cx_select);
      CHECK(checkInteract == 2); // 2^2

      particle.data_[0] = 100; // data positive --> Process1
      checkInteract = 0;
      sequence3_short.selectInteraction(view, projectileP4, noComposition, rng,
                                        cx_select);
      CHECK(checkInteract == 1); // 2^1
    }
  }

  SECTION("Check SecondariesProcesses in SwitchProcessSequence") {

    DummyData particle;
    DummyView view(particle);

    checkSecondaries = 0;
    particle.data_[0] = 100; // data positive  --> sequence1
    sequence3.doSecondaries(view);
    CHECK(checkSecondaries == 1);

    checkSecondaries = 0;
    particle.data_[0] = -100; // data positive  --> sequence1
    sequence3.doSecondaries(view);
    CHECK(checkSecondaries == 2);
  }

  SECTION("Check ContinuousProcesses in SwitchProcessSequence") {

    DummyData particle;
    DummyTrajectory track;
    Step step(particle, track);

    particle.data_[0] = 100; // data positive, selects particular branch on
                             // SwitchProcessSequence

    cp1.setStep(10_m);
    cp2.setStep(15_m);
    cp3.setStep(100_m);

    cp1.resetFlag();
    cp2.resetFlag();
    cp3.resetFlag();

    ContinuousProcessStepLength const step1 = sequence3.getMaxStepLength(particle, track);
    CHECK(LengthType(step1) == 10_m);
    sequence3.doContinuous(step, step1);
    CHECK(cp1.getFlag());
    CHECK_FALSE(cp2.getFlag());
    CHECK_FALSE(cp3.getFlag());
    CORSIKA_LOG_INFO("step1, l={}, i={}", LengthType(step1),
                     ContinuousProcessIndex(step1).getIndex());

    particle.data_[0] = 100; // data positive, selects particular branch on
                             // SwitchProcessSequence

    cp1.setStep(50_m);
    cp2.setStep(15_m);
    cp3.setStep(100_m);

    cp1.resetFlag();
    cp2.resetFlag();
    cp3.resetFlag();

    ContinuousProcessStepLength const step2 = sequence3.getMaxStepLength(particle, track);
    CHECK(LengthType(step2) == 15_m);
    sequence3.doContinuous(step, step2);
    CHECK_FALSE(cp1.getFlag());
    CHECK(cp2.getFlag());
    CHECK_FALSE(cp3.getFlag());
    CORSIKA_LOG_INFO("step2, len_cont={}, indexLimit={} type={}", LengthType(step2),
                     ContinuousProcessIndex(step2).getIndex(),
                     boost::typeindex::type_id<decltype(sequence3)>().pretty_name());

    particle.data_[0] = -100; // data positive, selects particular branch on
                              // SwitchProcessSequence

    cp1.setStep(11_m);
    cp2.setStep(15_m);
    cp3.setStep(100_m);

    cp1.resetFlag();
    cp2.resetFlag();
    cp3.resetFlag();

    ContinuousProcessStepLength const step3 = sequence3.getMaxStepLength(particle, track);
    CHECK(LengthType(step3) == 11_m);
    sequence3.doContinuous(step, step3);
    CHECK(cp1.getFlag());
    CHECK_FALSE(cp2.getFlag());
    CHECK_FALSE(cp3.getFlag());
    CORSIKA_LOG_INFO("step3, len_cont={}, indexLimit={} type={}", LengthType(step3),
                     ContinuousProcessIndex(step3).getIndex(),
                     boost::typeindex::type_id<decltype(sequence3)>().pretty_name());

    particle.data_[0] = -100; // data positive, selects particular branch on
                              // SwitchProcessSequence

    cp1.setStep(11_m);
    cp2.setStep(15_m);
    cp3.setStep(2_m);

    cp1.resetFlag();
    cp2.resetFlag();
    cp3.resetFlag();

    ContinuousProcessStepLength const step4 = sequence3.getMaxStepLength(particle, track);
    CHECK(LengthType(step4) == 2_m);
    sequence3.doContinuous(step, step4);
    CHECK_FALSE(cp1.getFlag());
    CHECK_FALSE(cp2.getFlag());
    CHECK(cp3.getFlag());
    CORSIKA_LOG_INFO("step4, len_cont={}, indexLimit={} type={}", LengthType(step4),
                     ContinuousProcessIndex(step4).getIndex(),
                     boost::typeindex::type_id<decltype(sequence3)>().pretty_name());
  }

  SECTION("Check BoundaryCrossingProcess in SwitchProcessSequence") {

    DummyData particle;
    DummyNode node_from(1);
    DummyNode node_to(2);

    particle.data_[0] = 100; // data positive, selects particular branch on
                             // SwitchProcessSequence

    sequence4.doBoundaryCrossing(particle, node_from, node_to);

    CHECK(particle.data_[0] == 97); // 100 - 2*1 - 1*1

    particle.data_[0] = -100; // data positive, selects particular branch on
                              // SwitchProcessSequence

    sequence4.doBoundaryCrossing(particle, node_from, node_to);
    CHECK(particle.data_[0] == -101); // -100 - 2*1 + 1*1
  }
}

TEST_CASE("ProcessSequence Indexing", "ProcessSequence") {

  logging::set_level(logging::level::info);

  SECTION("Indexing") {

    int const n0 = count_processes<Decay2>::count;
    int const n1 = count_processes<ContinuousProcess3>::count;
    int const n2 = count_processes<ContinuousProcess2,
                                   count_processes<ContinuousProcess3>::count>::count;
    int const n1_b =
        count_processes<Process2, count_processes<ContinuousProcess3>::count>::count;
    int const n1_c =
        count_processes<ContinuousProcess3, count_processes<Process2>::count>::count;
    int const n12 =
        count_processes<ContinuousProcess2,
                        count_processes<ContinuousProcess3, 10>::count>::count;
    int const n11_b =
        count_processes<Process1, count_processes<ContinuousProcess3, 10>::count>::count;
    int const n11_c =
        count_processes<ContinuousProcess3, count_processes<Process1, 10>::count>::count;

    CHECK(n0 == 1);
    CHECK(n1 == 1);
    CHECK(n1_b == 2);
    CHECK(n1_c == 2);
    CHECK(n2 == 2);
    CHECK(n11_b == 12);
    CHECK(n11_c == 12);
    CHECK(n12 == 12);

    std::cout << count_processes<ContinuousProcess3>::count << std::endl;
    std::cout << count_processes<Process3>::count << std::endl;

    struct SwitchSelect {
      bool operator()(DummyData const& p) const {
        std::cout << "SwitchSelect data=" << p.data_[0] << std::endl;
        return (p.data_[0] > 0);
      }
    };

    auto sequence1 = make_sequence(Process1(0), ContinuousProcess2(0, 2_m), Decay1(0));
    auto sequence2 = make_sequence(ContinuousProcess3(0, 3_m), Process2(0), Decay2(0),
                                   ContinuousProcess1(0, 1_m));

    SwitchSelect select1;
    auto switch_seq = SwitchProcessSequence(select1, sequence1, sequence2);

    auto sequence3 = make_sequence(ContinuousProcess1(0, 1_m), Process3(0), switch_seq);
    auto sequence4 = make_sequence(ContinuousProcess1(0, 1_m), Process3(0),
                                   SwitchProcessSequence(select1, sequence1, sequence2));

    int const switch_seq_n = count_processes<decltype(switch_seq)>::count;
    int const sequence3_n = count_processes<decltype(sequence3)>::count;

    CHECK(decltype(sequence1)::getNumberOfProcesses() == 3);
    CHECK(count_processes<decltype(sequence1)>::count == 3);
    CHECK(count_processes<decltype(sequence2)>::count == 4);
    CHECK(switch_seq_n == 7);
    CHECK(sequence3_n == 9);
    CHECK(count_processes<decltype(sequence4)>::count == 9);

    std::cout << "switch_seq "
              << boost::typeindex::type_id<decltype(switch_seq)>().pretty_name()
              << std::endl;

    std::cout << "sequence3 "
              << boost::typeindex::type_id<decltype(sequence3)>().pretty_name()
              << std::endl;
  }
}

class ProcessZero : public InteractionProcess<ProcessZero> {
public:
  ProcessZero() = default;

  template <typename TView>
  void doInteraction(TView&, Code const, Code const, FourMomentum const&,
                     FourMomentum const&) const {
    FAIL("ProcessZero::doInteraction has been called!");
  }

  CrossSectionType getCrossSection(Code const, Code const, FourMomentum const&,
                                   FourMomentum const&) const {
    CORSIKA_LOG_DEBUG("ProcessZero::getCrossSection");
    return 0_mb;
  }
};

TEST_CASE("SelectInteractionZeroCrossSection", "ProcessSequence") {
  logging::set_level(logging::level::info);

  auto sequence = make_sequence(ProcessZero());

  FourMomentum const projectileP4{10_GeV, {0_eV, 0_eV, 0_eV}};

  DummyData particle;
  DummyView view(particle);

  CrossSectionType cx_select =
      sequence.getCrossSection(particle, Code::Nitrogen, projectileP4);
  CHECK(cx_select == 0_mb); // should be zero
  media::NuclearComposition const noComposition({Code::Nitrogen}, {1});
  DummyRNG rng;

  auto retValue =
      sequence.selectInteraction(view, projectileP4, noComposition, rng, 0_mb);
  CHECK(!isInteracted(retValue)); // cross section of process sequence is zero,
                                  // no process should cause an interaction
}

TEST_CASE("SwitchProcessSequence Indexing", "ProcessSequence") {
  // see issue https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/issues/573
  // and issue https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/issues/586

  logging::set_level(logging::level::info);

  struct SwitchSelect {
    bool operator()(DummyData const& p) const { return (p.data_[0] > 0); }
  };
  SwitchSelect select1;

  auto cp1 = ContinuousProcess1(0, 0_m);
  auto cp2 = ContinuousProcess1(0, 0_m);
  auto cp3 = ContinuousProcess1(0, 0_m);
  auto cp4 = ContinuousProcess1(0, 0_m);

  auto switch_sequence = make_select(select1, cp1, cp2);
  auto sequence = make_sequence(switch_sequence, cp3, cp4);

  CHECK(sequence.getNumberOfProcesses() == 4);
  CHECK(switch_sequence.getNumberOfProcesses() == 2);

  DummyData particle;
  DummyTrajectory track;
  DummyView view(particle);
  Step step(particle, track);

  // cp1 selected
  cp1.setStep(1_m);
  cp2.setStep(100_m);
  cp3.setStep(100_m);
  cp4.setStep(100_m);
  particle.data_[0] = 1; // positive so that cp1 is selected
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  ContinuousProcessStepLength continuousMaxStep =
      sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  cp1.resetFlag();
  cp2.resetFlag();
  cp3.resetFlag();
  cp4.resetFlag();

  // cp2 selected
  cp1.setStep(100_m);
  cp2.setStep(1_m);
  cp3.setStep(100_m);
  cp4.setStep(100_m);
  particle.data_[0] = -1; // negative so that cp2 is selected
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  continuousMaxStep = sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK_FALSE(cp1.getFlag());
  CHECK(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  cp1.resetFlag();
  cp2.resetFlag();
  cp3.resetFlag();
  cp4.resetFlag();

  // cp3 selected
  cp1.setStep(100_m);
  cp2.setStep(100_m);
  cp3.setStep(1_m);
  cp4.setStep(100_m);
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  continuousMaxStep = sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  cp1.resetFlag();
  cp2.resetFlag();
  cp3.resetFlag();
  cp4.resetFlag();

  // cp4 selected
  cp1.setStep(100_m);
  cp2.setStep(100_m);
  cp3.setStep(100_m);
  cp4.setStep(1_m);
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK_FALSE(cp4.getFlag());
  continuousMaxStep = sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  CHECK(cp4.getFlag());
  cp1.resetFlag();
  cp2.resetFlag();
  cp3.resetFlag();
  cp4.resetFlag();
}

TEST_CASE("Nested ProcessSequence", "ProcessSequence") {
  //

  logging::set_level(logging::level::info);

  auto cp1 = ContinuousProcess1(0, 0_m);
  auto cp2 = ContinuousProcess1(0, 0_m);
  auto cp3 = ContinuousProcess1(0, 0_m);

  auto sequence = make_sequence(make_sequence(cp1, cp2), cp3);

  CHECK(sequence.getNumberOfProcesses() == 3);

  DummyData particle;
  DummyTrajectory track;
  DummyView view(particle);
  Step step(particle, track);

  // cp1 selected
  cp1.setStep(1_m);
  cp2.setStep(100_m);
  cp3.setStep(100_m);
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  ContinuousProcessStepLength continuousMaxStep =
      sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  cp1.resetFlag();
  cp2.resetFlag();
  cp3.resetFlag();

  // cp2 selected
  cp1.setStep(100_m);
  cp2.setStep(1_m);
  cp3.setStep(100_m);
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  continuousMaxStep = sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK_FALSE(cp1.getFlag());
  CHECK(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  cp1.resetFlag();
  cp2.resetFlag();
  cp3.resetFlag();

  // cp3 selected
  cp1.setStep(100_m);
  cp2.setStep(100_m);
  cp3.setStep(1_m);
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK_FALSE(cp3.getFlag());
  continuousMaxStep = sequence.getMaxStepLength(particle, track);
  sequence.doContinuous(step, continuousMaxStep.getIndex());
  CHECK_FALSE(cp1.getFlag());
  CHECK_FALSE(cp2.getFlag());
  CHECK(cp3.getFlag());
}
