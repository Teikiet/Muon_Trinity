/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/CubicSolver.hpp>
#include <cmath>

namespace corsika {

  namespace andre {

    inline std::vector<double> solve_quartic_real(long double a, long double b,
                                                  long double c, long double d,
                                                  long double e, double const epsilon) {

      if (std::abs(a) < epsilon) { return solve_cubic_real(b, c, d, e, epsilon); }

      b /= a;
      c /= a;
      d /= a;
      e /= a;

      long double a3 = -c;
      long double b3 = b * d - 4. * e;
      long double c3 = -b * b * e - d * d + 4. * c * e;

      // cubic resolvent
      // y^3 − c*y^2 + (bd−4e)*y − b^2*e−d^2+4*c*e = 0

      std::vector<double> x3 = solve_cubic_real(1, a3, b3, c3, epsilon);
      if (!x3.size()) {
        return {}; // no solution, numeric problem (LCOV_EXCL_LINE)
      }
      long double y = x3[0]; // there is always at least one solution
      // The essence - choosing Y with maximal absolute value.
      if (x3.size() == 3) {
        if (std::abs(x3[1]) > std::abs(y)) y = x3[1];
        if (std::abs(x3[2]) > std::abs(y)) y = x3[2];
      }

      long double q1, q2, p1, p2;
      // h1+h2 = y && h1*h2 = e  <=>  h^2 -y*h + e = 0    (h === q)

      long double Det = y * y - 4 * e;
      CORSIKA_LOG_TRACE("Det={}", Det);
      if (std::abs(Det) < epsilon) // in other words - D==0
      {
        q1 = q2 = y * 0.5;
        // g1+g2 = b && g1+g2 = c-y   <=>   g^2 - b*g + c-y = 0    (p === g)
        Det = b * b - 4 * (c - y);
        CORSIKA_LOG_TRACE("Det={}", Det);
        if (std::abs(Det) < epsilon) { // in other words - D==0
          p1 = p2 = b * 0.5;
        } else {
          if (Det < 0) return {};
          long double sqDet = sqrt(Det);
          p1 = (b + sqDet) * 0.5;
          p2 = (b - sqDet) * 0.5;
        }
      } else {
        if (Det < 0) return {};
        long double sqDet1 = sqrt(Det);
        q1 = (y + sqDet1) * 0.5;
        q2 = (y - sqDet1) * 0.5;
        // g1+g2 = b && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
        p1 = (b * q1 - d) / (q1 - q2);
        p2 = (d - b * q2) / (q1 - q2);
      }

      // solving quadratic eqs.  x^2 + p1*x + q1 = 0
      //                         x^2 + p2*x + q2 = 0

      std::vector<double> quad1 = solve_quadratic_real(1, p1, q1, 1e-5);
      std::vector<double> quad2 = solve_quadratic_real(1, p2, q2, 1e-5);
      if (quad2.size() > 0) {
        for (auto const val : quad2) quad1.push_back(val);
      }
      return quad1;
    }
  } // namespace andre

  inline std::vector<double> solve_quartic_depressed_real(long double p, long double q,
                                                          long double r,
                                                          double const epsilon) {

    CORSIKA_LOG_TRACE("quartic-depressed: p={:f}, q={:f}, r={:f},  epsilon={}", p, q, r,
                      epsilon);

    long double const p2 = static_pow<2>(p);
    long double const q2 = static_pow<2>(q);

    std::vector<double> const resolve_cubic =
        solve_cubic_real(1, p, p2 / 4 - r, -q2 / 8, epsilon);

    CORSIKA_LOG_TRACE("resolve_cubic: N={}, m=[{}]", resolve_cubic.size(),
                      fmt::join(resolve_cubic, ", "));

    if (!resolve_cubic.size()) return {};

    long double m = 0;
    for (auto const& v : resolve_cubic) {
      CORSIKA_LOG_TRACE("check pol3(v)={}", (static_pow<3>(v) + static_pow<2>(v) * p +
                                             v * (p2 / 4 - r) - q2 / 8));
      if (std::abs(v) > epsilon && std::abs(v) > m) { m = v; }
    }
    CORSIKA_LOG_TRACE("check m={}", m);
    if (m == 0) { return {0}; }
    if (m < 0) {
      // this is a rare numerical instability
      // first: try analytical solution, second: discard (curved->straight tracking)
      std::vector<double> const resolve_cubic_analytic =
          andre::solve_cubic_real_analytic(1, p, p2 / 4 - r, -q2 / 8, epsilon);

      CORSIKA_LOG_TRACE("andre::resolve_cubic_analytic: N={}, m=[{}]",
                        resolve_cubic_analytic.size(),
                        fmt::join(resolve_cubic_analytic, ", "));

      if (!resolve_cubic_analytic.size()) return {};

      for (auto const& v : resolve_cubic_analytic) {
        CORSIKA_LOG_TRACE("check pol3(v)={}", (static_pow<3>(v) + static_pow<2>(v) * p +
                                               v * (p2 / 4 - r) - q2 / 8));
        if (std::abs(v) > epsilon && std::abs(v) > m) { m = v; }
      }
      CORSIKA_LOG_TRACE("check m={}", m);
      if (m == 0) { return {0}; }
      if (m < 0) {
        return {}; // now we are out of options, cannot solve: curved->straight tracking
      }
    }
    long double const quad_term1 = p / 2 + m;
    long double const quad_term2 = std::sqrt(2 * m);
    long double const quad_term3 = q / (2 * quad_term2);

    std::vector<double> z_quad1 =
        solve_quadratic_real(1, quad_term2, quad_term1 - quad_term3, 1e-5);
    std::vector<double> z_quad2 =
        solve_quadratic_real(1, -quad_term2, quad_term1 + quad_term3, 1e-5);
    for (auto const& z : z_quad2) z_quad1.push_back(z);
    return z_quad1;
  }

  inline std::vector<double> solve_quartic_real(long double a, long double b,
                                                long double c, long double d,
                                                long double e, double const epsilon) {

    CORSIKA_LOG_TRACE("quartic: a={:f}, b={:f}, c={:f}, d={:f}, e={:f}, epsilon={}", a, b,
                      c, d, e, epsilon);

    if (std::abs(a) < epsilon) { // this is just a quadratic
      return solve_cubic_real(b, c, d, e, epsilon);
    }

    if ((std::abs(a - 1) < epsilon) &&
        (std::abs(b) < epsilon)) { // this is a depressed quartic
      return solve_quartic_depressed_real(c, d, e, epsilon);
    }

    long double const b2 = static_pow<2>(b);
    long double const b3 = static_pow<3>(b);
    long double const b4 = static_pow<4>(b);
    long double const a2 = static_pow<2>(a);
    long double const a3 = static_pow<3>(a);
    long double const a4 = static_pow<4>(a);

    long double const p = (c * a * 8 - b2 * 3) / (a4 * 8);
    long double const q = (b3 - b * c * a * 4 + d * a2 * 8) / (a4 * 8);
    long double const r =
        (-b4 * 3 + e * a3 * 256 - b * d * a2 * 64 + b2 * c * a * 16) / (a4 * 256);

    std::vector<double> zs = solve_quartic_depressed_real(p, q, r, epsilon);
    CORSIKA_LOG_TRACE("quartic: solve_depressed={}, b/4a={}", fmt::join(zs, ", "),
                      b / (4 * a));
    for (auto& z : zs) { z -= b / (4 * a); }
    CORSIKA_LOG_TRACE("quartic: solve_quartic_real returns={}", fmt::join(zs, ", "));
    return zs;
  }
} // namespace corsika
