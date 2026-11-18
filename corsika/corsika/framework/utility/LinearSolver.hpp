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

  std::vector<double> solve_linear_real(double a, double b);

  std::vector<std::complex<double>> solve_linear(double a, double b);

} // namespace corsika

#include <corsika/detail/framework/utility/LinearSolver.inl>
