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

#include <corsika/framework/utility/LinearSolver.hpp>

namespace corsika {

  std::vector<std::complex<double>> solve_quadratic(long double a, long double b,
                                                    long double c,
                                                    double const epsilon = 1e-12);

  std::vector<double> solve_quadratic_real(long double a, long double b, long double c,
                                           double const epsilon = 1e-12);
} // namespace corsika

#include <corsika/detail/framework/utility/QuadraticSolver.inl>
