/*
 * (c) Copyright 2019 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/modules/ExecTime.hpp>

#include <corsika/process/example_processors/DummyBoundaryCrossingProcess.hpp>
#include <corsika/process/example_processors/DummyContinuousProcess.hpp>
#include <corsika/process/example_processors/DummyDecayProcess.hpp>
#include <corsika/process/example_processors/DummyInteractionProcess.hpp>
#include <corsika/process/example_processors/DummySecondariesProcess.hpp>

#include <cmath>
#include <random>
#include <vector>

using namespace corsika::process;
using namespace corsika::process::analytic_processors;
using namespace corsika::process::example_processors;

TEST_CASE("Timing process", "[proccesses][analytic_processors ExecTime]") {

  logging::set_level(logging::level::info);

  int tmp = 0;

  SECTION("BoundaryCrossing") {
    ExecTime<DummyBoundaryCrossingProcess<10>> execTime;
    auto start = std::chrono::steady_clock::now();
    CHECK(execTime.DoBoundaryCrossing(tmp, 0, 0) == EProcessReturn::eOk);
    auto end = std::chrono::steady_clock::now();
    CHECK(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() ==
          Approx(10).margin(5));

    for (int i = 0; i < 100; i++) execTime.DoBoundaryCrossing(tmp, 0, 0);

    CHECK(execTime.mean() == Approx(10 * 1000).margin(2 * 1000));

    CHECK(execTime.sumTime() == Approx(10 * 100 * 1000).margin((10 * 100) * 1000));

    CHECK(fabs(execTime.var()) < 100000);
  }

  SECTION("Continuous") {
    ExecTime<DummyContinuousProcess<50>> execTime;
    auto start = std::chrono::steady_clock::now();
    CHECK(execTime.DoContinuous(tmp, tmp) == EProcessReturn::eOk);
    auto end = std::chrono::steady_clock::now();
    CHECK(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() ==
          Approx(50).margin(5));

    for (int i = 0; i < 100; i++) execTime.DoContinuous(tmp, tmp);

    CHECK(execTime.mean() == Approx(50 * 1000).margin(2 * 1000));

    CHECK(execTime.sumTime() == Approx(50 * 100 * 1000).margin((10 * 100) * 1000));

    CHECK(fabs(execTime.var()) < 100000);
  }

  SECTION("Decay") {
    ExecTime<DummyDecayProcess<10>> execTime;
    auto start = std::chrono::steady_clock::now();
    CHECK(execTime.DoDecay(tmp) == EProcessReturn::eOk);
    auto end = std::chrono::steady_clock::now();
    CHECK(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() ==
          Approx(10).margin(5));

    for (int i = 0; i < 100; i++) execTime.DoDecay(tmp);

    CHECK(execTime.mean() == Approx(10 * 1000).margin(2 * 1000));

    CHECK(execTime.sumTime() == Approx(10 * 100 * 100).margin((10 * 100) * 1000));

    CHECK(fabs(execTime.var()) < 100000);
  }

  SECTION("Interaction") {
    ExecTime<DummyInteractionProcess<10>> execTime;
    auto start = std::chrono::steady_clock::now();
    CHECK(execTime.DoInteraction(tmp) == EProcessReturn::eOk);
    auto end = std::chrono::steady_clock::now();
    CHECK(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() ==
          Approx(10).margin(5));

    for (int i = 0; i < 100; i++) execTime.DoInteraction(tmp);

    CHECK(execTime.mean() == Approx(10 * 1000).margin(2 * 1000));

    CHECK(execTime.sumTime() == Approx(10 * 100 * 1000).margin((10 * 100) * 1000));

    CHECK(fabs(execTime.var()) < 100000);
  }

  SECTION("Secondaries") {
    ExecTime<DummySecondariesProcess<10>> execTime;
    auto start = std::chrono::steady_clock::now();
    CHECK(execTime.DoSecondaries(tmp) == EProcessReturn::eOk);
    auto end = std::chrono::steady_clock::now();
    CHECK(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() ==
          Approx(10).margin(5));

    for (int i = 0; i < 100; i++) execTime.DoSecondaries(tmp);

    CHECK(execTime.mean() == Approx(10 * 1000).margin(2 * 1000));

    CHECK(execTime.sumTime() == Approx(10 * 100 * 1000).margin((10 * 100) * 1000));

    CHECK(fabs(execTime.var()) < 100000);
  }

  SECTION("TestMeanAlgo") {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(10000.0, 200.0);

    double fElapsedSum = 0;
    double fMean = 0;
    double fMean2 = 0;
    long long fMin = std::numeric_limits<long long>::max();
    long long fMax = std::numeric_limits<long long>::min();
    int fN = 0;

    std::vector<double> elems;

    for (int i = 0; i < 1000; i++) {
      double timeDiv = distribution(generator);

      elems.push_back(timeDiv);

      fElapsedSum += timeDiv;
      fN = fN + 1;

      if (fMax < timeDiv) fMax = timeDiv;

      if (timeDiv < fMin) fMin = timeDiv;

      double delta = timeDiv - fMean;
      fMean += delta / static_cast<double>(fN);

      double delta2 = timeDiv - fMean;

      fMean2 += delta * delta2;
    }

    CHECK(fN == 1000);

    double mean = 0;
    std::for_each(elems.begin(), elems.end(), [&](double i) { mean += i; });
    mean = mean / fN;

    double var = 0;
    std::for_each(elems.begin(), elems.end(),
                  [&](double i) { var += (mean - i) * (mean - i); });
    var = var / fN;

    CHECK(mean == Approx(10000.0).margin(10));
    CHECK(var == Approx(200.0 * 200).margin(2000));

    CHECK(fMean2 / fN == Approx(200 * 200).margin(2000)); // Varianz
    CHECK(fMean == Approx(10000).margin(10));
  }
}
