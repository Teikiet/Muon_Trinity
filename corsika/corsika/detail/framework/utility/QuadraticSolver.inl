/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  inline std::vector<std::complex<double>> solve_quadratic(long double a, long double b,
                                                           long double c,
                                                           double const epsilon) {
    if (std::abs(a) < epsilon) { return solve_linear(b, c); }

    if (std::abs(c) < epsilon) {
      std::vector<std::complex<double>> lin_result = solve_linear(a, b);
      lin_result.push_back({0.});
      return lin_result;
    }

    long double const radicant = static_pow<2>(b) - a * c * 4;

    if (radicant < -epsilon) { // two complex solutions
      double const rpart = -b / 2 * a;
      double const ipart = std::sqrt(-radicant);
      return {{rpart, ipart}, {rpart, -ipart}};
    }

    if (radicant < epsilon) { // just one real solution
      return {{double(-b / 2 * a), 0}};
    }

    // two real solutions, use Citardauq formula and actively avoid cancellation in the
    // denominator

    const long double x1 =
        2 * c / (b > 0 ? -b - std::sqrt(radicant) : -b + std::sqrt(radicant));

    return {{double(x1), 0}, {double(c / (a * x1)), 0}};
  }

  inline std::vector<double> solve_quadratic_real(long double a, long double b,
                                                  long double c, double const epsilon) {

    CORSIKA_LOG_TRACE("quadratic: a={} b={} c={}", a, b, c);

    if (std::abs(a) < epsilon) { return solve_linear_real(b, c); }
    if (std::abs(c) < epsilon) {
      std::vector<double> lin_result = solve_linear_real(a, b);
      lin_result.push_back(0.);
      return lin_result;
    }

    // long double const radicant = std::pow(b, 2) - a * c * 4;
    // Thanks!
    // https://math.stackexchange.com/questions/2434354/numerically-stable-scheme-for-the-3-real-roots-of-a-cubic
    long double w = a * 4 * c;
    long double e = std::fma(a * 4, c, -w);
    long double f = std::fma(b, b, -w);
    long double radicant = f + e;

    CORSIKA_LOG_TRACE("radicant={} {} ", radicant, b * b - a * c * 4);

    if (std::abs(radicant) < epsilon) { // just one real solution
      return {double(-b / (2 * a))};
    }

    if (radicant < 0) { return {}; }

    // two real solutions, use Citardauq formula and actively avoid cancellation in the
    // denominator

    long double const x1 =
        c * 2 / (b > 0 ? -b - std::sqrt(radicant) : -b + std::sqrt(radicant));

    CORSIKA_LOG_TRACE("quadratic x1={} x2={}", double(x1), double(c / (a * x1)));

    return {double(x1), double(c / (a * x1))};
  }
} // namespace corsika
