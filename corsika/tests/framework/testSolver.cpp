/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/utility/QuadraticSolver.hpp>
#include <corsika/framework/utility/LinearSolver.hpp>
#include <corsika/framework/utility/CubicSolver.hpp>
#include <corsika/framework/utility/QuarticSolver.hpp>

using namespace corsika;
using namespace std;
using Catch::Approx;

double pol4(long double x, long double a, long double b, long double c, long double d,
            long double e) {
  return std::pow(x, 4) * a + std::pow(x, 3) * b + std::pow(x, 2) * c + x * d + e;
}

void remove_duplicates(vector<double>& v, double const eps) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end(),
                      [eps](auto v1, auto v2) {
                        return (std::abs(v2) > eps ? std::abs(2 * (v1 - v2) / (v1 + v2)) <
                                                         eps                   // relative
                                                   : std::abs(v1 - v2) < eps); // absolute
                      }),
          v.end());
}

TEST_CASE("Solver") {

  logging::set_level(logging::level::info);

  double epsilon_check = 1e-3; // for catch2 asserts

  SECTION("linear") {

    std::vector<std::pair<double, double>> vals = {{13.33, 8338.3},
                                                   {-333.99, -8.15},
                                                   {-58633.3, 2343.54},
                                                   {87632.87, -982.37},
                                                   {7e8, -1e-7}};

    for (auto v : vals) {
      {
        double a = v.first;
        double b = v.second;

        vector<double> s1 = solve_linear_real(a, b);

        CORSIKA_LOG_INFO("linear: a={}, b={}, N={}, s1[0]={}", a, b, s1.size(), s1[0]);

        double expected = -b / a;

        CHECK(s1.size() == 1);
        CHECK((s1[0] == Approx(expected).epsilon(epsilon_check)));

        vector<complex<double>> s2 = solve_linear(a, b);
        CHECK(s2.size() == 1);
        CHECK((s2[0].real() == Approx(expected).epsilon(epsilon_check)));
      }
    }

    CHECK(solve_linear_real(0, 55.).size() == 0);
    CHECK(solve_linear(0, 55.).size() == 0);

  } // linear

  SECTION("quadratic") {

    // tests of type:  (x-z1)*(x-z2) = 0  --> x1=z1, x2=z2   and    a=1, b=-z2-z1,
    // c=z1*z2

    std::vector<std::vector<double>> zeros = {
        {13.33, 8338.3}, {-333.99, -8.15}, {-58633.3, 2343.54}, {87632.87, -982.37},
        {1.3e-6, 5e7},   {-4.2e-7, 65e6},  {0.1, -13e-5},       {7e8, -1e-7}};

    for (unsigned int idegree = 0; idegree < 2; ++idegree) {
      CORSIKA_LOG_INFO("-------------     quadratic idegree={}", idegree);
      for (auto z : zeros) {

        {
          long double const z1 = z[0];
          long double const z2 = idegree <= 0 ? z1 : z[1];

          long double a = 1;
          long double b = -z1 - z2;
          long double c = z1 * z2;

          {
            std::vector<double> s1 = solve_quadratic_real(a, b, c);
            remove_duplicates(s1, epsilon_check);

            CORSIKA_LOG_INFO("quadratic: z1={}, z2={}, N={}, s1=[{}]", z1, z2, s1.size(),
                             fmt::join(s1, ", "));

            CHECK(s1.size() == idegree + 1);
            for (auto value : s1) {
              if (std::abs(value) < epsilon_check) {
                CHECK(((value == Approx(z1).margin(epsilon_check)) ||
                       (value == Approx(z2).margin(epsilon_check))));
              } else {
                CHECK(((value == Approx(z1).epsilon(epsilon_check)) ||
                       (value == Approx(z2).epsilon(epsilon_check))));
              }
            }
          }

          // cubic with x^3 * 0 should result in the same
          {
            std::vector<double> s1 = solve_cubic_real(0, a, b, c);
            remove_duplicates(s1, epsilon_check);

            CORSIKA_LOG_INFO("quadratic/cubic: z1={}, z2={}, N={}, s1=[{}]", z1, z2,
                             s1.size(), fmt::join(s1, ", "));

            CHECK(s1.size() == idegree + 1);
            for (auto value : s1) {
              if (std::abs(value) < epsilon_check) {
                CHECK(((value == Approx(z1).margin(epsilon_check)) ||
                       (value == Approx(z2).margin(epsilon_check))));
              } else {
                CHECK(((value == Approx(z1).epsilon(epsilon_check)) ||
                       (value == Approx(z2).epsilon(epsilon_check))));
              }
            }
          }

          // quartic with x^4 * 0 + x^3 * 0 must be the same
          {
            std::vector<double> s1 = solve_quartic_real(0, 0, a, b, c);
            remove_duplicates(s1, epsilon_check);

            CORSIKA_LOG_INFO("quadratic/quartic: z1={}, z2={}, N={}, s1=[{}]", z1, z2,
                             s1.size(), fmt::join(s1, ", "));

            CHECK(s1.size() == idegree + 1);
            for (auto value : s1) {
              if (std::abs(value) < epsilon_check) {
                CHECK(((value == Approx(z1).margin(epsilon_check)) ||
                       (value == Approx(z2).margin(epsilon_check))));
              } else {
                CHECK(((value == Approx(z1).epsilon(epsilon_check)) ||
                       (value == Approx(z2).epsilon(epsilon_check))));
              }
            }
          }
        }
      }
    }
  } // quadratic

  SECTION("cubic") {

    // tests of type:
    // (x-z1) * (x-z2) * (x-z3) = 0  --> x1=z1, x2=z2, x3=z3   and

    // (x^2 - x*(z1+z2) + z1*z2) * (x-z3)   =
    //  x^3 - x^2*(z1+z2) + x*z1*z2 - x^2 z3 + x*(z1+z2)*z3 - z1*z2*z3   =
    //  x^3 + x^2 (-z3-z1-z2) + x (z1*z2 + (z1+z2)*z3) - z1*z2*z3

    epsilon_check = 1e-2; // for catch2 asserts

    std::vector<std::vector<double>> zeros = {
        {13.33, 838.3, 44.},     {-333.99, -8.15, -33e8}, {-58633.3, 2343.54, -1e-5},
        {87632.87, -982.37, 0.}, {1.3e-4, 5e7, 6.6e9},    {-4.2e-7, 65e6, 433.3334},
        {23e-1, -13e-2, 10.},    {7e8, -1e-7, 1e8}};

    for (unsigned int idegree = 0; idegree < 3; ++idegree) {
      CORSIKA_LOG_INFO("-------------     cubic idegree={}", idegree);
      for (auto z : zeros) {

        {
          long double const z1 = z[0];
          long double const z2 = idegree <= 0 ? z1 : z[1];
          long double const z3 = idegree <= 1 ? z1 : z[2];

          long double const a = 1;
          long double const b = -z1 - z2 - z3;
          long double const c = z1 * z2 + (z1 + z2) * z3;
          long double const d = -z1 * z2 * z3;
          //
          CORSIKA_LOG_INFO(
              "cubic: z1={}, z2={}, z3={}, "
              "a={}, b={}, c={}, d={}",
              z1, z2, z3, a, b, c, d);

          {
            vector<double> s1 = solve_cubic_real(a, b, c, d);
            remove_duplicates(s1, epsilon_check * 10);

            CORSIKA_LOG_INFO("N={}, s1=[{}]", s1.size(), fmt::join(s1, ", "));

            CHECK(s1.size() == idegree + 1);
            for (double value : s1) {
              CORSIKA_LOG_INFO("value={}, z1={} z2={} z3={} eps_check={}", value, z1, z2,
                               z3, epsilon_check);
              if (std::abs(value) < epsilon_check) {
                CHECK(((value == Approx(z1).margin(epsilon_check)) ||
                       (value == Approx(z2).margin(epsilon_check)) ||
                       (value == Approx(z3).margin(epsilon_check))));
              } else {
                CHECK(((value == Approx(z1).epsilon(epsilon_check)) ||
                       (value == Approx(z2).epsilon(epsilon_check)) ||
                       (value == Approx(z3).epsilon(epsilon_check))));
              }
            }
          }

          // quartic with x^4 *0  must be the same
          {
            vector<double> s1 = solve_quartic_real(0, a, b, c, d);
            remove_duplicates(s1, epsilon_check * 10);

            CORSIKA_LOG_INFO("(quartic) N={}, s1=[{}]", s1.size(), fmt::join(s1, ", "));

            CHECK(s1.size() == idegree + 1);
            for (double value : s1) {
              CORSIKA_LOG_INFO("value={}, z1={} z2={} z3={} eps_check={}", value, z1, z2,
                               z3, epsilon_check);
              if (std::abs(value) < epsilon_check) {
                CHECK(((value == Approx(z1).margin(epsilon_check)) ||
                       (value == Approx(z2).margin(epsilon_check)) ||
                       (value == Approx(z3).margin(epsilon_check))));
              } else {
                CHECK(((value == Approx(z1).epsilon(epsilon_check)) ||
                       (value == Approx(z2).epsilon(epsilon_check)) ||
                       (value == Approx(z3).epsilon(epsilon_check))));
              }
            }
          }
        }
      }
    }
  } // cubic

  SECTION("quartic") {

    epsilon_check = 1e-2; // for catch2 asserts

    // **clang-format-off**
    // tests of type:
    // (x-z1) * (x-z2) * (x-z3) * (x-z4) = 0  --> x1=z1, x2=z2, x3=z3, x4=z4   and

    // (x^2 - x   (z1+z2)     +      z1*z2) * (x-z3)                            *
    // (x-z4) = (x^3 - x^2 (z1+z2)     + x    z1*z2  - x^2 z3 + x*(z1+z2)*z3 -
    // z1*z2*z3) * (x-z4) = (x^3 + x^2 (-z1-z2-z3) + x   (z1*z2 + (z1+z2)*z3) -
    // z1*z2*z3) * (x-z4)
    // =
    //
    //  x^4 + x^3 (-z1-z2-z3) + x^2 (z1*z2 + (z1+z2)*z3) - x z1*z2*z3
    //      - x^3 z4          - x^2 (-z1-z2-z3)*z4       - x (z1*z2 + (z1+z2)*z3)*z4
    //      + z1*z2*z3*z4
    //
    // = x^4
    // + x^3 (-z1-z2-z3-z4)
    // + x^2 (z1*z2 + (z1+z2)*z3 + (z1+z2+z3)*z4)
    // - x (z1*z2*z3 + (z1*z2 + (z1+z2)*z3)*z4))
    // + z1*z2*z3*z4
    // **clang-format-on**

    std::vector<std::vector<double>> zeros = {
        {13.33, 838.3, 44., 2.3},       {-3333.99, -8.15, -33e4, 8.8e3},
        {-58633.3, 2343.54, -1e-1, 0.}, {87632.87, -982.37, 10., 1e-2},
        {1.3e2, 5e5, 6.6e5, 1e3},       {-4.9, 65e2, 433.3334, 6e5},
        {23e-1, -13e-2, 10., 3.4e6},    {7e6, -1e-1, 1e5, 2.55e4}};

    for (unsigned int idegree = 2; idegree < 4;
         ++idegree) { // idegree==1 is not very stable !!!!!!!
      CORSIKA_LOG_INFO("-------------     quartic idegree={}", idegree);
      for (auto z : zeros) {

        {
          long double const z1 = z[0];
          long double const z2 = idegree <= 0 ? z1 : z[1];
          long double const z3 = idegree <= 1 ? z1 : z[2];
          long double const z4 = idegree <= 2 ? z1 : z[3];

          long double const a = 1;
          long double const b = -z1 - z2 - z3 - z4;
          long double const c = z1 * z2 + (z1 + z2) * z3 + (z1 + z2 + z3) * z4;
          long double const d = -z1 * z2 * z3 - (z1 * z2 + (z1 + z2) * z3) * z4;
          long double const e = z1 * z2 * z3 * z4;

          //
          CORSIKA_LOG_INFO(
              "quartic: z1={}, z2={}, z3={}, z4={},  "
              "a={}, b={}, c={}, d={}, e={}",
              z1, z2, z3, z4, a, b, c, d, e);

          // make sure the tested cases are all ZERO (printout)
          CORSIKA_LOG_INFO("quartic trace: {} {} {} {}", pol4(z1, a, b, c, d, e),
                           pol4(z2, a, b, c, d, e), pol4(z3, a, b, c, d, e),
                           pol4(z4, a, b, c, d, e));

          vector<double> s1 = andre::solve_quartic_real(a, b, c, d, e);
          remove_duplicates(s1, epsilon_check * 10);
          vector<double> s2 = solve_quartic_real(a, b, c, d, e);
          remove_duplicates(s2, epsilon_check * 5);

          CORSIKA_LOG_INFO("N={}, s1=[{}]", s1.size(), fmt::join(s1, ", "));
          CORSIKA_LOG_INFO("N={}, s2=[{}]", s2.size(), fmt::join(s2, ", "));

          CHECK(s1.size() == idegree + 1);
          for (double value : s1) {
            CORSIKA_LOG_INFO("value={}, z1={} z2={} z3={} z4={} eps_check={}", value, z1,
                             z2, z3, z4, epsilon_check);
            if (std::abs(value) < epsilon_check) {
              CHECK(((value == Approx(z1).margin(epsilon_check)) ||
                     (value == Approx(z2).margin(epsilon_check)) ||
                     (value == Approx(z3).margin(epsilon_check)) ||
                     (value == Approx(z4).margin(epsilon_check))));
            } else {
              CHECK(((value == Approx(z1).epsilon(epsilon_check)) ||
                     (value == Approx(z2).epsilon(epsilon_check)) ||
                     (value == Approx(z3).epsilon(epsilon_check)) ||
                     (value == Approx(z4).epsilon(epsilon_check))));
            }
          }

          // this is a bit less precise
          CHECK(s2.size() == idegree + 1);
          for (double value : s2) {
            CORSIKA_LOG_INFO("value={}, z1={} z2={} z3={} z4={} eps_check={}", value, z1,
                             z2, z3, z4, epsilon_check);
            if (std::abs(value) < epsilon_check) {
              CHECK(((value == Approx(z1).margin(epsilon_check * 5)) ||
                     (value == Approx(z2).margin(epsilon_check * 5)) ||
                     (value == Approx(z3).margin(epsilon_check * 5)) ||
                     (value == Approx(z4).margin(epsilon_check * 5))));
            } else {
              CHECK(((value == Approx(z1).epsilon(epsilon_check * 5)) ||
                     (value == Approx(z2).epsilon(epsilon_check * 5)) ||
                     (value == Approx(z3).epsilon(epsilon_check * 5)) ||
                     (value == Approx(z4).epsilon(epsilon_check * 5))));
            }
          }
        }
      }
    }

  } // quartic
}
