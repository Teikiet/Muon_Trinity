/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#define protected public // to also test the internal state of objects

#include <corsika/framework/stack/Stack.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <testTestStack.hpp> // from tests/common

#include <iomanip>
#include <tuple>
#include <vector>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace std;

typedef Stack<TestStackData, TestParticleInterface> StackTest;

TEST_CASE("Stack", "[Stack]") {

  logging::set_level(logging::level::info);

  // helper function for sum over stack data
  auto sum = [](StackTest const& stack) {
    double v = 0;
    for (const auto& p : stack) v += p.getData();
    return v;
  };

  SECTION("StackInterface") {

    // construct a valid Stack object
    StackTest stack;
    stack.clear();
    CHECK(stack.getSize() == 0);
    CHECK(stack.isEmpty());                          // stack empty here
    auto pTest0 = stack.addParticle(std::tuple{0.}); // [0]
    CHECK(stack.getSize() == 1);
    CHECK(!stack.isEmpty());
    auto pTest1 = stack.addParticle(std::tuple{1.}); // [0,1]
    CHECK(stack.getSize() == 2);
    CHECK(pTest1.getData() == 1.);
    auto pTestAt = stack.at(1); // -> 1
    CHECK(pTestAt == pTest1);
    auto pTestFirst = stack.first(); // -> 0
    CHECK(pTestFirst == pTest0);
  }

  SECTION("construct") {

    // construct a valid, empty Stack object
    StackTest s;
  }

  SECTION("write and read") {

    StackTest stack;
    stack.addParticle(std::tuple{9.9});
    const double v = sum(stack);
    CHECK(v == 9.9);
    auto const particle_const = stack.begin();
    CHECK(particle_const.getData() == 9.9);
  }

  SECTION("delete from stack") {

    StackTest stack;

    CHECK_THROWS(stack.erase(stack.begin())); // nothing to delete

    CHECK(stack.getSize() == 0);
    StackTest::stack_iterator_type p =
        stack.addParticle(std::tuple{0.}); // valid way to access particle data
    p.setData(9.9);
    CHECK(stack.getSize() == 1);
    CHECK(stack.getEntries() == 1);
    CHECK(!p.isErased());
    stack.erase(p);
    CHECK(p.isErased());
    CHECK(stack.getSize() == 1);
    CHECK(stack.getEntries() == 0);
  }

  SECTION("delete particle") {

    StackTest stack;
    CHECK(stack.getSize() == 0);
    stack.addParticle(std::tuple{8.9});
    stack.addParticle(std::tuple{7.9});
    auto p = stack.addParticle(
        std::tuple{9.9}); // also valid way to access particle data, identical to above

    CHECK(stack.getSize() == 3);
    CHECK(stack.getEntries() == 3);
    CHECK(!stack.isEmpty());

    p.erase();               // mark for deletion: size=3, entries=2
    CHECK_THROWS(p.erase()); // already deleted
    CHECK(stack.getSize() == 3);
    CHECK(stack.getEntries() == 2);
    CHECK(!stack.isEmpty());

    stack.last().erase(); // mark for deletion: size=3, entries=1
    CHECK(stack.getSize() == 3);
    CHECK(stack.getEntries() == 1);
    CHECK(!stack.isEmpty());

    /*
       GetNextParticle will find two entries marked as "deleted" and
       will purge this from the end of the stack: size = 1
    */
    stack.getNextParticle().erase(); // mark for deletion: size=3, entries=0
    CHECK(stack.getSize() == 1);
    CHECK(stack.getEntries() == 0);
    CHECK(stack.isEmpty());
  }

  SECTION("create secondaries") {

    StackTest stack;
    CHECK(stack.getSize() == 0);
    auto iter = stack.addParticle(std::tuple{9.9});
    StackTest::particle_interface_type& p =
        *iter; // also this is valid to access particle data
    CHECK(stack.getSize() == 1);
    p.addSecondary(std::tuple{4.4});
    CHECK(stack.getSize() == 2);
  }

  SECTION("get next particle") {
    StackTest stack;
    CHECK(stack.getSize() == 0);
    CHECK(stack.getEntries() == 0);
    CHECK(stack.isEmpty());

    stack.addParticle(std::tuple{9.9});
    stack.addParticle(std::tuple{8.8});
    CHECK(stack.getSize() == 2);
    CHECK(stack.getEntries() == 2);
    CHECK(!stack.isEmpty());

    auto particle = stack.getNextParticle(); // first particle
    CHECK(particle.getData() == 8.8);

    particle.erase(); // only marks (last) particle as deleted
    CHECK(stack.getSize() == 2);
    CHECK(stack.getEntries() == 1);
    CHECK(!stack.isEmpty());

    /*
      This following call to GetNextParticle will realize that the
      current last particle on the stack was marked "deleted" and will
      purge it: stack size is reduced by one.
     */
    auto particle2 = stack.getNextParticle(); // first particle
    CHECK(particle2.getData() == 9.9);
    CHECK(stack.getSize() == 1);
    CHECK(stack.getEntries() == 1);
    CHECK(!stack.isEmpty());

    particle2.erase(); // also mark this particle as deleted

    CHECK(stack.getSize() == 1);
    CHECK(stack.getEntries() == 0);
    CHECK(stack.isEmpty());
  }

  SECTION("swap particle") {
    StackTest stack;
    CHECK(stack.getSize() == 0);
    CHECK(stack.getEntries() == 0);
    CHECK(stack.isEmpty());

    stack.addParticle(std::tuple{9.888});
    stack.addParticle(std::tuple{8.999});
    CHECK(stack.getSize() == 2);
    CHECK(stack.getEntries() == 2);
    CHECK(!stack.isEmpty());

    auto p1 = stack.begin();
    auto p2 = p1 + 1;

    CHECK(p1.getData() == 9.888);
    CHECK(p2.getData() == 8.999);

    stack.swap(p1, p2);

    CHECK(p1.getData() == 8.999);
    CHECK(p2.getData() == 9.888);
  }

  SECTION("copy particle") {
    StackTest stack;
    CHECK(stack.getSize() == 0);
    CHECK(stack.getEntries() == 0);
    CHECK(stack.isEmpty());

    stack.addParticle(std::tuple{9.888});
    stack.addParticle(std::tuple{8.999});
    CHECK(stack.getSize() == 2);
    CHECK(stack.getEntries() == 2);
    CHECK(!stack.isEmpty());

    auto p1 = stack.begin();
    auto p2 = p1 + 1;

    CHECK(p1.getData() == 9.888);
    CHECK(p2.getData() == 8.999);

    stack.copy(p1, p2);

    CHECK(p1.getData() == 9.888);
    CHECK(p2.getData() == 9.888);
  }
}
