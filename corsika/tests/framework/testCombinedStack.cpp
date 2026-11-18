/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#define protected public // to also test the internal state of objects

#include <corsika/framework/stack/CombinedStack.hpp>
#include <corsika/framework/stack/SecondaryView.hpp>
#include <corsika/framework/stack/Stack.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <testTestStack.hpp> // for testing: simple stack. This is a
// test-build, and inluce file is obtained from CMAKE_CURRENT_SOURCE_DIR

#include <iomanip>
#include <vector>

#include <catch2/catch_all.hpp>
using Catch::Approx;

using namespace corsika;
using namespace std;

////////////////////////////////////////////////////////////
// first level test: combine two stacks:
//                   StackTest = (TestStackData + TestStackData2)

// definition of stack-data object
class TestStackData2 {

public:
  // these functions are needed for the Stack interface
  void clear() { data2_.clear(); }
  unsigned int getSize() const { return data2_.size(); }
  unsigned int getCapacity() const { return data2_.size(); }
  void copy(const int i1, const int i2) { data2_[i2] = data2_[i1]; }
  void swap(const int i1, const int i2) {
    double tmp0 = data2_[i1];
    data2_[i1] = data2_[i2];
    data2_[i2] = tmp0;
  }

  // custom data access function
  void setData2(const int i, const double v) { data2_[i] = v; }
  double getData2(const int i) const { return data2_[i]; }

  // these functions are also needed by the Stack interface
  void incrementSize() { data2_.push_back(0.); }
  void decrementSize() {
    if (data2_.size() > 0) { data2_.pop_back(); }
  }

  // custom private data section
private:
  std::vector<double> data2_;
};

// defintion of a stack-readout object, the iteractor dereference
// operator will deliver access to these function
template <typename T>
class TestParticleInterface2 : public T {

public:
  using T::getIndex;
  using T::getStackData;
  using T::setParticleData;

  // default version for particle-creation from input data
  void setParticleData(std::tuple<double> const v = {0.}) { setData2(std::get<0>(v)); }
  void setParticleData(TestParticleInterface2<T>& parent,
                       std::tuple<double> const v = {0.}) {
    setData2(parent.getData2() + std::get<0>(v));
  }
  void setData2(const double v) { getStackData().setData2(getIndex(), v); }
  double getData2() const { return getStackData().getData2(getIndex()); }
};

// combined stack: StackTest = (TestStackData + TestStackData2)
template <typename TStackIter>
using CombinedTestInterfaceType =
    corsika::CombinedParticleInterface<TestParticleInterface, TestParticleInterface2,
                                       TStackIter>;

using StackTest = CombinedStack<TestStackData, TestStackData2, CombinedTestInterfaceType>;

TEST_CASE("Combined Stack", "[stack]") {

  logging::set_level(logging::level::info);

  // helper function for sum over stack data
  auto sum = [](const StackTest& stack) {
    double v = 0;
    for (const auto& p : stack) v += p.getData();
    return v;
  };
  auto sum2 = [](const StackTest& stack) {
    double v = 0;
    for (const auto& p : stack) v += p.getData2();
    return v;
  };

  SECTION("StackInterface") {

    // construct a valid Stack object
    StackTest s;
    s.clear();
    s.addParticle(std::tuple{0.});
    s.copy(s.cbegin(), s.begin());
    s.swap(s.begin(), s.begin());
    CHECK(s.getSize() == 1);
  }

  SECTION("construct") {

    // construct a valid, empty Stack object
    StackTest s;
  }

  SECTION("write and read") {

    StackTest s;
    s.addParticle(std::tuple{9.9});
    CHECK(sum2(s) == 0.);
    CHECK(sum(s) == 9.9);
  }

  SECTION("delete from stack") {

    StackTest s;
    CHECK(s.getSize() == 0);
    StackTest::stack_iterator_type p =
        s.addParticle(std::tuple{0.}); // valid way to access particle data
    p.setData(8.9);
    p.setData2(3.);
    CHECK(sum2(s) == 3.);
    CHECK(sum(s) == 8.9);
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 1);
    s.erase(p);
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 0);
  }

  SECTION("delete particle") {

    StackTest s;
    CHECK(s.getSize() == 0);
    auto p = s.addParticle(
        std::tuple{9.9}); // also valid way to access particle data, identical to above
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 1);
    p.erase();
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 0);
  }

  SECTION("create secondaries") {
    StackTest s;
    CHECK(s.getSize() == 0);
    auto iter = s.addParticle(std::tuple{9.9});
    iter.setData2(2);
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 1);
    iter.addSecondary(std::tuple{4.4});
    CHECK(s.getSize() == 2);
    CHECK(s.getEntries() == 2);
    // p.addSecondary(3.3, 2.2, 1.);
    // CHECK(s.getSize() == 3);
    double v = 0;
    for (const auto& i : s) {
      v += i.getData();
      CHECK(i.getData2() == 2);
    }
    CHECK(v == 9.9 + 4.4);
  }

  SECTION("get next particle") {
    StackTest s;
    CHECK(s.getSize() == 0);
    CHECK(s.getEntries() == 0);
    CHECK(s.isEmpty());

    auto p1 = s.addParticle(std::tuple{9.9});
    auto p2 = s.addParticle(std::tuple{8.8});
    p1.setData2(20.2);
    p2.setData2(20.3);
    CHECK(s.getSize() == 2);
    CHECK(s.getEntries() == 2);
    CHECK(!s.isEmpty());

    auto particle = s.getNextParticle(); // first particle
    CHECK(particle.getData() == 8.8);
    CHECK(particle.getData2() == 20.3);

    particle.erase(); // only marks (last) particle as "deleted"
    CHECK(s.getSize() == 2);
    CHECK(s.getEntries() == 1);
    CHECK(!s.isEmpty());

    /*
      This following call to GetNextParticle will realize that the
      current last particle on the stack was marked "deleted" and will
      purge it: stack size is reduced by one.
     */
    auto particle2 = s.getNextParticle(); // first particle
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 1);
    CHECK(!s.isEmpty());
    CHECK(particle2.getData() == 9.9);
    CHECK(particle2.getData2() == 20.2);

    particle2.erase(); // also mark this particle as "deleted"
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 0);
    CHECK(s.isEmpty());
  }

  SECTION("exceptions") {
    StackTest s;
    auto p1 = s.addParticle(std::tuple{9.9});
    auto p2 = s.addParticle(std::tuple{9.9});
    ++p2;
    CHECK_THROWS(s.copy(p1, p2));
    CHECK_THROWS(s.swap(p1, p2));
    CHECK(s.getSize() == 2);
  }
}

////////////////////////////////////////////////////////////
// next level: combine three stacks:
// combined stack: StackTest2 = ((TestStackData + TestStackData2) + TestStackData3)

// definition of stack-data object
class TestStackData3 {

public:
  // these functions are needed for the Stack interface
  void clear() { data3_.clear(); }
  unsigned int getSize() const { return data3_.size(); }
  unsigned int getCapacity() const { return data3_.size(); }
  void copy(const int i1, const int i2) { data3_[i2] = data3_[i1]; }
  void swap(const int i1, const int i2) {
    double tmp0 = data3_[i1];
    data3_[i1] = data3_[i2];
    data3_[i2] = tmp0;
  }

  // custom data access function
  void setData3(const int i, const double v) { data3_[i] = v; }
  double getData3(const int i) const { return data3_[i]; }

  // these functions are also needed by the Stack interface
  void incrementSize() { data3_.push_back(0.); }
  void decrementSize() {
    if (data3_.size() > 0) { data3_.pop_back(); }
  }

  // custom private data section
private:
  std::vector<double> data3_;
};

// ---------------------------------------
// defintion of a stack-readout object, the iteractor dereference
// operator will deliver access to these function
template <typename T>
class TestParticleInterface3 : public T {

public:
  using T::getIndex;
  using T::getStackData;
  using T::setParticleData;

  // default version for particle-creation from input data
  void setParticleData(std::tuple<double> const v = {0.}) { setData3(std::get<0>(v)); }
  void setParticleData(TestParticleInterface3<T>& parent,
                       std::tuple<double> const v = {0.}) {
    setData3(parent.getData3() + std::get<0>(v));
  }
  void setData3(const double v) { getStackData().setData3(getIndex(), v); }
  double getData3() const { return getStackData().getData3(getIndex()); }
};

// double combined stack:
// combined stack
template <typename TStackIter>
using CombinedTestInterfaceType2 =
    corsika::CombinedParticleInterface<StackTest::pi_type, TestParticleInterface3,
                                       TStackIter>;

using StackTest2 = CombinedStack<typename StackTest::stack_data_type, TestStackData3,
                                 CombinedTestInterfaceType2>;

TEST_CASE("Combined Stack - multi", "[stack]") {

  logging::set_level(logging::level::info);

  SECTION("create secondaries") {

    StackTest2 s;
    CHECK(s.getSize() == 0);
    CHECK(s.isEmpty()); // size = entries = 0

    // add new particle, only provide tuple data for StackTest
    auto p1 = s.addParticle(std::tuple{9.9});
    // add new particle, provide tuple data for both StackTest and TestStackData3
    auto p2 = s.addParticle(std::tuple{8.8}, std::tuple{0.1});

    CHECK(s.getSize() == 2);
    CHECK(!s.isEmpty()); // size = entries = 2

    // examples to explicitly change data on stack
    p2.setData2(0.1); // not clear why this is needed, need to check
                      // SetParticleData workflow for more complicated
                      // settings
    p1.setData3(20.2);
    p2.setData3(10.3);

    CHECK(p1.getData() == 9.9);
    CHECK(p1.getData2() == 0.);
    p1.setData2(10.2);
    CHECK(p1.getData2() == 10.2);
    CHECK(p1.getData3() == 20.2);

    CHECK(p2.getData() == 8.8);
    CHECK(p2.getData2() == 0.1);
    CHECK(p2.getData3() == 10.3);

    auto particle = s.getNextParticle(); // first particle
    CHECK(particle.getData() == 8.8);
    CHECK(particle.getData2() == 0.1);
    CHECK(particle.getData3() == 10.3);

    auto sec = particle.addSecondary(std::tuple{4.4});
    CHECK(s.getSize() == 3);
    CHECK(s.getEntries() == 3);
    CHECK(sec.getData() == 4.4);
    CHECK(sec.getData2() == 0.1);
    CHECK(sec.getData3() == 10.3);

    sec.erase(); // mark for deletion: size=3, entries=2
    CHECK(s.getSize() == 3);
    CHECK(s.getEntries() == 2);
    CHECK(!s.isEmpty());

    s.last().erase(); // mark for deletion: size=3, entries=1
    CHECK(s.getSize() == 3);
    CHECK(s.getEntries() == 1);
    CHECK(!s.isEmpty());

    /*
       GetNextParticle will find two entries marked as "deleted" and
       will purge this from the end of the stack: size = 1
    */
    s.getNextParticle().erase(); // mark for deletion: size=3, entries=0
    CHECK(s.getSize() == 1);
    CHECK(s.getEntries() == 0);
    CHECK(s.isEmpty());
  }
}

////////////////////////////////////////////////////////////

// final level test, create SecondaryView on StackTest2

/*
  See Issue 161

  unfortunately clang does not support this in the same way (yet) as
  gcc, so we have to distinguish here. If clang cataches up, we could
  remove the clang branch here and also in corsika::Cascade. The gcc
  code is much more generic and universal.
 */
template <typename TStackIter>
using CombinedTestInterfaceType2 =
    corsika::CombinedParticleInterface<StackTest::pi_type, TestParticleInterface3,
                                       TStackIter>;

using StackTest2 = CombinedStack<typename StackTest::stack_data_type, TestStackData3,
                                 CombinedTestInterfaceType2>;

#if defined(__clang__)
using StackTestView =
    SecondaryView<typename StackTest2::stack_data_type, CombinedTestInterfaceType2>;
#elif defined(__GNUC__) || defined(__GNUG__)
using StackTestView = corsika::MakeView<StackTest2>::type;
#endif

using Particle2 = typename StackTest2::particle_type;

TEST_CASE("Combined Stack - secondary view") {

  logging::set_level(logging::level::info);

  SECTION("create secondaries via secondaryview") {

    StackTest2 stack;
    auto particle = stack.addParticle(std::tuple{9.9});
    StackTestView view(particle);

    auto projectile = view.getProjectile();
    projectile.addSecondary(std::tuple{8.8});

    CHECK(stack.getSize() == 2);
  }
}
