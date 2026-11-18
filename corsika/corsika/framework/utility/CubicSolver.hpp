/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>

#include <corsika/framework/utility/QuadraticSolver.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>

/**
 *   @file CubicSolver.hpp
 */

namespace corsika {

  /**
   *  @ingroup Utility
   * @{
   */

  namespace andre {
    /**
       Solve a x^3 + b x^2 + c x + d = 0
    */
    std::vector<double> solveP3(long double a, long double b, long double c,
                                long double d, double const epsilon = 1e-12);
  } // namespace andre

  /**
       Solve depressed cubic: x^3 + p x + q = 0
   */
  inline std::vector<double> solve_cubic_depressed_disciminant_real(
      long double p, long double q, long double const disc, double const epsilon = 1e-12);

  std::vector<double> solve_cubic_depressed_real(long double p, long double q,
                                                 double const epsilon = 1e-12);

  /**
     Solve a x^3 + b x^2 + c x + d = 0

     Analytical approach. Not very stable in all conditions.
   */
  std::vector<double> solve_cubic_real_analytic(long double a, long double b,
                                                long double c, long double d,
                                                double const epsilon = 1e-12);

  /**
   * Cubic function.
   *
   * T must be a floating point type.
   *
   *  @return a x^3 + b x^2 + c x + d
   */
  template <typename T>
  T cubic_function(T x, T a, T b, T c, T d);

  /**
     @return 3 a x^2 + 2 b x + c
  */
  template <typename T>
  T cubic_function_dfdx(T x, T a, T b, T c);

  /**
     @return 6 a x + 2 b
  */
  template <typename T>
  T cubic_function_d2fd2x(T x, T a, T b);

  /**
     Iterative approach to solve: a x^3 + b x^2 + c x + d = 0

     This is often fastest and most precise.
   */
  std::vector<double> solve_cubic_real(long double a, long double b, long double c,
                                       long double d, double const epsilon = 1e-12);

  //! @}

} // namespace corsika

#include <corsika/detail/framework/utility/CubicSolver.inl>
