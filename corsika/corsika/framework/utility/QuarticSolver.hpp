/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/utility/CubicSolver.hpp>
#include <corsika/framework/utility/QuadraticSolver.hpp>

#include <complex>
#include <algorithm>
#include <cmath>
#include <numeric>

/**
 * @file QuarticSolver.hpp
 */

namespace corsika {

  /**
   * @ingroup Utilities
   * @{
   */

  namespace andre {

    /**
     solve quartic equation a*x^4 + b*x^3 + c*x^2 + d*x + e
     Attention - this function returns dynamically allocated array. It has to be
     released afterwards.
    */
    std::vector<double> solve_quartic_real(long double a, long double b, long double c,
                                           long double d, long double e,
                                           double const epsilon = 1e-12);
  } // namespace andre

  /**
     solve quartic equation a*x^4 + b*x^3 + c*x^2 + d*x + e
  */
  std::vector<double> solve_quartic_real(long double a, long double b, long double c,
                                         long double d, long double e,
                                         double const epsilon = 1e-12);

  //! @}

} // namespace corsika

#include <corsika/detail/framework/utility/QuarticSolver.inl>
