/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#define protected public // to also test the internal state of objects

#include <corsika/framework/stack/SecondaryView.hpp>
#include <corsika/framework/stack/Stack.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <testTestStack.hpp> // for testing: simple stack. This is a
// test-build, and inluce file is obtained from CMAKE_CURRENT_SOURCE_DIR

#include <iomanip>
#include <vector>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace std;

typedef Stack<TestStackData, TestParticleInterface> StackTest;

/*
  See Issue 161

  unfortunately clang does not support this in the same way (yet) as
  gcc, so we have to distinguish here. If clang cataches up, we could
  remove the clang branch here and also in corsika::Cascade. The gcc
  code is much more generic and universal.
 */
#if defined(__clang__)
using StackTestView = SecondaryView<TestStackData, TestParticleInterface>;
#elif defined(__GNUC__) || defined(__GNUG__)
using StackTestView = MakeView<StackTest>::type;
#endif

using Particle = typename StackTest::particle_type;

TEST_CASE("SecondaryStack", "[stack]") {

  logging::set_level(logging::level::info);

  // helper function for sum over stack data
  auto sum = [](const StackTest& stack) {
    double v = 0;
    for (const auto& p : stack) { v += p.getData(); }
    return v;
  };

  auto sumView = [](const StackTestView& stack) {
    double value = 0;
    for (const auto& p : stack) { value += p.getData(); }
    return value;
  };

  SECTION("secondary view") {
    StackTest stack;
    CHECK(stack.getSize() == 0);
    CHECK(stack.isEmpty());

    stack.addParticle(std::tuple{9.9});
    stack.addParticle(std::tuple{8.8});
    const double sumS = 9.9 + 8.8; // helper, see below
    CHECK(stack.getSize() == 2);
    CHECK(stack.getEntries() == 2);
    CHECK(!stack.isEmpty());

    auto particle = stack.getNextParticle();

    StackTestView view(particle);
    CHECK(view.getSize() == 0);
    CHECK(view.getEntries() == 0);
    CHECK(view.isEmpty());
    CHECK_THROWS(view.erase(view.begin())); // nothing to delete yet

    {
      auto proj = view.getProjectile();
      CHECK(proj.getData() == particle.getData());
      proj.addSecondary(std::tuple{4.4});
    }
    CHECK(view.getSize() == 1);
    CHECK(view.getEntries() == 1);
    CHECK(!view.isEmpty());
    CHECK(stack.getSize() == 3);
    CHECK(stack.getEntries() == 3);
    CHECK(!stack.isEmpty());

    view.addSecondary(std::tuple{4.5});
    view.addSecondary(std::tuple{4.6});
    CHECK(view.getSize() == 3);
    CHECK(view.getEntries() == 3);
    CHECK(!view.isEmpty());
    CHECK(stack.getSize() == 5);
    CHECK(stack.getEntries() == 5);
    CHECK(!stack.isEmpty());

    CHECK(sum(stack) == sumS + 4.4 + 4.5 + 4.6);
    CHECK(sumView(view) == 4.4 + 4.5 + 4.6);

    view.last().erase();
    CHECK(view.getSize() == 3);
    CHECK(view.getEntries() == 2);
    CHECK(stack.getSize() == 5);
    CHECK(stack.getEntries() == 4);

    CHECK(sum(stack) == sumS + 4.4 + 4.5);
    CHECK(sumView(view) == 4.4 + 4.5);

    auto pDel = view.getNextParticle();
    view.erase(pDel);
    CHECK_THROWS(view.erase(pDel)); // already erased
    CHECK(view.getSize() == 2);
    CHECK(stack.getSize() == 4);

    CHECK(sum(stack) == sumS + 4.4 + 4.5 - pDel.getData());
    CHECK(sumView(view) == 4.4 + 4.5 - pDel.getData());

    view.erase(view.getNextParticle());
    CHECK(sum(stack) == sumS);
    CHECK(sumView(view) == 0);
    CHECK(view.isEmpty());

    {
      auto proj = view.getProjectile();
      CHECK(proj.getData() == particle.getData());
      CHECK(particle == view.parent());
    }

    // at() and first()
    auto pTestAt = view.at(0);
    auto pTestFirst = view.first();
    CHECK(pTestFirst == pTestAt);
  }

  SECTION("secondary view, construct from ParticleType") {
    StackTest stack;
    CHECK(stack.getSize() == 0);
    stack.addParticle(std::tuple{9.9});
    stack.addParticle(std::tuple{8.8});

    auto iterator = stack.getNextParticle();
    typename StackTest::particle_type& particle = iterator; // as in corsika::Cascade

    StackTestView view(particle);
    CHECK(view.getSize() == 0);

    view.addSecondary(std::tuple{4.4});

    CHECK(view.getSize() == 1);
  }

  SECTION("deletion") {
    StackTest stack;
    stack.addParticle(std::tuple{-99.});
    stack.addParticle(std::tuple{0.});

    {
      auto particle = stack.getNextParticle();
      StackTestView view(particle);

      auto proj = view.getProjectile();
      proj.addSecondary(std::tuple{-2.});
      proj.addSecondary(std::tuple{-1.});
      proj.addSecondary(std::tuple{1.});
      proj.addSecondary(std::tuple{2.});

      CHECK(stack.getSize() == 6); // -99, 0, -2, -1, 1, 2
      CHECK(view.getSize() == 4);  // -2, -1, 1, 2

      // now delete all negative entries, i.e. -1 and -2
      auto p = view.begin();
      while (p != view.end()) {
        auto data = p.getData();
        if (data < 0) { p.erase(); }
        ++p;
      }
      CHECK(stack.getSize() == 6);
      CHECK(stack.getEntries() == 4); // -99, 0, 2, 1 (order changes during deletion)
      CHECK(view.getSize() == 4);
      CHECK(view.getEntries() == 2); // 2, 1
    }

    // repeat

    {
      auto particle = stack.getNextParticle();
      StackTestView view(particle);

      // put -2,...,+2 on stack
      auto proj = view.getProjectile();
      proj.addSecondary(std::tuple{-2.});
      proj.addSecondary(std::tuple{-1.});
      proj.addSecondary(std::tuple{1.});
      proj.addSecondary(std::tuple{2.});
      // stack should contain -99, 0, 2, 1, [-2, -1, 1, 2]

      auto p = view.begin();
      while (p != view.end()) {
        auto data = p.getData();
        if (data < 0) { p.erase(); }
        ++p;
      }

      // stack should contain -99, 0, 2, 1, [2, 1]
      // view should contain 1, 2

      CHECK(stack.getEntries() == 6);
      CHECK(stack.getSize() == 10);
    }
  }

  SECTION("swap particle") {
    StackTest stack;
    stack.addParticle(std::tuple{-99.});

    StackTestView view(stack.first());
    view.addSecondary(std::tuple{-2.});
    view.addSecondary(std::tuple{-1.});
    view.addSecondary(std::tuple{1.});

    auto p1 = view.begin();
    auto p2 = p1 + 1;

    CHECK(p1.getData() == -2.);
    CHECK(p2.getData() == -1.);

    view.swap(p1, p2);

    CHECK(p1.getData() == -1);
    CHECK(p2.getData() == -2);
  }

  SECTION("copy particle") {
    StackTest stack;
    stack.addParticle(std::tuple{-99.});

    StackTestView view(stack.first());
    view.addSecondary(std::tuple{-2.});
    view.addSecondary(std::tuple{-1.});
    view.addSecondary(std::tuple{1.});

    auto p1 = view.begin();
    auto p2 = p1 + 1;

    CHECK(p1.getData() == -2.);
    CHECK(p2.getData() == -1.);

    view.copy(p1, p2);

    CHECK(p1.getData() == -2);
    CHECK(p2.getData() == -2);
  }
}
