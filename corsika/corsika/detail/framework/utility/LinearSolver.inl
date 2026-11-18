/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <vector>
#include <cmath>
#include <complex>

namespace corsika {

  inline std::vector<double> solve_linear_real(double a, double b) {

    if (a == 0) {
      return {}; // no (b!=0), or infinite number (b==0) of solutions....
    }

    return {-b / a};
  }

  inline std::vector<std::complex<double>> solve_linear(double a, double b) {

    if (std::abs(a) == 0) {
      return {}; // no (b!=0), or infinite number (b==0) of solutions....
    }

    return {std::complex<double>(-b / a, 0)};
  }

} // namespace corsika
