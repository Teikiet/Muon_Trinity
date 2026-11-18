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

    //---------------------------------------------------------------------------
    // solve cubic equation A x^3 + B*x^2 + C*x + D = 0
    //                        x^3 + a*x^2 + b*x + c = 0
    // mainly along WolframAlpha formulas
    inline std::vector<double> solve_cubic_real_analytic(long double A, long double B,
                                                         long double C, long double D,
                                                         double const epsilon) {

      if (std::abs(A) < epsilon) { return solve_quadratic_real(B, C, epsilon); }

      long double a = B / A;
      long double b = C / A;
      long double c = D / A;

      long double a2 = a * a;
      long double q = (3 * b - a2) / 9;
      long double r = (a * (9 * b - 2 * a2) - 27 * c) / 54;
      long double q3 = q * q * q;

      // disc = q**3 + r**2
      // w:e = r*r exactly
      long double w = r * r;
      long double e = std::fma(r, r, -w);
      // s:t =  q*q exactly
      long double s = q * q;
      long double t = std::fma(q, q, -s);
      // s:t * q + w:e = s*q + w + t*q +e = s*q+w + u:v + e = f + u:v + e
      long double f = std::fma(s, q, w);
      long double u = t * q;
      long double v = std::fma(t, q, -u);
      // error-free sum f+u. See Knuth, TAOCP vol. 2
      long double uf = u + f;
      long double au = uf - u;
      long double ab = uf - au;
      au = f - au;
      ab = u - ab;
      // sum all terms into final result
      long double const disc = (((e + uf) + au) + ab) + v;

      CORSIKA_LOG_TRACE("disc={} {}", disc, q3 + r * r);

      if (std::abs(disc) < epsilon) {

        a /= 3;
        long double const cbrtR = std::cbrt(r);
        return {double(2 * cbrtR - a), double(-cbrtR - a)}; // 2nd solution is doublet

      } else if (disc > 0) {

        long double const S = std::cbrt(r + std::sqrt(disc));
        long double const T = std::cbrt(r - std::sqrt(disc));
        a /= 3;
        return {double((S + T) - a)}; // plus two imaginary solution

      } else { // disc < 0

        long double t = r / std::sqrt(-q3);
        if (t < -1) t = -1;
        if (t > 1) t = 1;
        t = std::acos(t);
        a /= 3;
        q = 2 * std::sqrt(-q);
        return {double(q * std::cos(t / 3) - a),
                double(q * std::cos((t + 2 * M_PI) / 3) - a),
                double(q * std::cos((t + 4 * M_PI) / 3) - a)};
      }
    }
  } // namespace andre

  inline std::vector<double> solve_cubic_depressed_disciminant_real(
      long double p, long double q, long double const disc, double const epsilon) {

    CORSIKA_LOG_TRACE("p={}, q={}, disc={}", p, q, disc);

    if (std::abs(disc) < epsilon) { // disc==0   multiple roots !
      if (std::abs(p) < epsilon) {  // tripple root
        return {0};
      }
      // double root, single root
      CORSIKA_LOG_TRACE("cubic double root");
      return {double(-3 * q / (2 * p)), double(3 * q / p)};
    }

    if (std::abs(p) < epsilon) { // p==0  --> x^3 + q = 0
      return {double(std::cbrt(-q))};
    }

    if (disc > 0) { // three real roots
      CORSIKA_LOG_TRACE("cubic three roots");
      long double const p_third = p / 3;
      long double const sqrt_minus_p_third = std::sqrt(-p_third);

      long double const arg = std::acos(q / (2 * p_third) / sqrt_minus_p_third) / 3;

      long double constexpr pi = M_PI;
      return {double(2 * sqrt_minus_p_third * std::cos(arg)),
              double(2 * sqrt_minus_p_third * std::cos(arg - 2 * pi / 3)),
              double(2 * sqrt_minus_p_third * std::cos(arg - 4 * pi / 3))};
    }

    // thus disc < 0  ->  one real root
    if (p < 0) {
      CORSIKA_LOG_TRACE("cubic cosh");
      long double const abs_q = std::abs(q);
      long double const p_third = p / 3;
      long double const sqrt_minus_p_third = std::sqrt(-p_third);
      CORSIKA_LOG_TRACE("sqrt_minus_p_third={}, arcosh({})={}", sqrt_minus_p_third,
                        -abs_q / (2 * p_third) / sqrt_minus_p_third,
                        std::acosh(-abs_q / (2 * p_third) / sqrt_minus_p_third));
      CORSIKA_LOG_TRACE(
          "complex: {}",
          -2 * abs_q / q * sqrt_minus_p_third *
              std::cosh(std::acosh(-abs_q / (2 * p_third * sqrt_minus_p_third)) / 3));
      double const z =
          -2 * abs_q / q * sqrt_minus_p_third *
          std::cosh(std::acosh(-abs_q / (2 * p_third * sqrt_minus_p_third)) / 3);
      return {z};
    } else { // p > 0
      CORSIKA_LOG_TRACE("cubic sinh");
      long double const p_third = p / 3;
      long double const sqrt_p_third = std::sqrt(p_third);
      return {double(-2 * sqrt_p_third *
                     std::sinh(std::asinh(q / (2 * p_third * sqrt_p_third)) / 3))};
    }
  }

  inline std::vector<double> solve_cubic_depressed_real(long double p, long double q,
                                                        double const epsilon) {

    // thanks!
    // https://math.stackexchange.com/questions/2434354/numerically-stable-scheme-for-the-3-real-roots-of-a-cubic
    // long double const disc = -(4 * p * p * p + 27 * q * q);
    // disc = (p/3)**3 + (q/2)**2
    long double p_third = p / 3;
    long double q_half = q / 2;
    // w:e = (q/2)*(q/2) exactly
    long double w = q_half * q_half;
    long double e = std::fma(q_half, q_half, -w);
    // s:t =  (p/3)*(p/3) exactly
    long double s = p_third * p_third;
    long double t = std::fma(p_third, p_third, -s);
    // s:t * p + w:e = s*p + w + t*p +e = s*p+w + u:v + e = f + u:v + e
    long double f = std::fma(s, p_third, w);
    long double u = t * p_third;
    long double v = std::fma(t, p_third, -u);
    // error-free sum f+u. See Knuth, TAOCP vol. 2
    long double a = u + f;
    long double b = a - u;
    long double c = a - b;
    b = f - b;
    c = u - c;
    // sum all terms into final result
    long double const disc = -(((e + a) + b) + c) + v;
    return solve_cubic_depressed_disciminant_real(p, q, disc, epsilon);
  }

  /**
   * Analytical approach. Not very stable in all conditions.
   */
  inline std::vector<double> solve_cubic_real_analytic(long double a, long double b,
                                                       long double c, long double d,
                                                       double const epsilon) {

    CORSIKA_LOG_TRACE("cubic: a={:f}, b={:f}, c={:f}, d={:f}, epsilon={} {} {}", a, b, c,
                      d, epsilon, (std::abs(a - 1) < epsilon), (std::abs(b) < epsilon));

    if (std::abs(a) < epsilon) { // this is just a quadratic
      return solve_quadratic_real(b, c, d, epsilon);
    }

    if ((std::abs(a - 1) < epsilon) &&
        (std::abs(b) < epsilon)) { // this is a depressed cubic
      return solve_cubic_depressed_real(c, d, epsilon);
    }

    // p = (3ac - b^2) / (3a^2) = 3 * ( c/(3*a) - b^2/(9*a^2) )
    long double b_over_a = b / a;
    long double const p_third = std::fma(-b_over_a, b_over_a / 9, c / (a * 3));

    // p = std::fma(a * 3, c, -b2) / (3 * a2);
    // q = (2*b^3 - 9*abc + 27*a^2*d) / (27a^3) = 2 * ( b^3/(27*a^3) - bc/(6*a^2) +
    // d/(2*a) )
    long double const q_half_term1 = std::fma(b_over_a, b_over_a / 27, -c / (a * 6));
    long double const q_half = std::fma(b_over_a, q_half_term1, d / (a * 2));

    std::vector<double> zs = solve_cubic_depressed_real(3 * p_third, 2 * q_half, epsilon);
    CORSIKA_LOG_TRACE("cubic: solve_depressed={}, b/3a={}", fmt::join(zs, ", "),
                      b / (a * 3));
    for (auto& z : zs) {
      z -= b / (a * 3);
      double const b1 = z + b / a;
      double const b0 = b1 * z + c / a;
      std::vector<double> quad_check = solve_quadratic_real(1, b1, b0, epsilon);
      CORSIKA_LOG_TRACE("quad_check=[{}], f(z)={}", fmt::join(quad_check, ", "),
                        static_pow<3>(z) * a + static_pow<2>(z) * b + z * c + d);
    }
    CORSIKA_LOG_TRACE("cubic: solve_cubic_real returns={}", fmt::join(zs, ", "));
    return zs;
  }

  template <typename T> // T must be floating point type
  inline T cubic_function(T x, T a, T b, T c, T d) {
    T x2 = x * x;
    return x2 * x * a + x2 * b + x * c + d;
  }

  template <typename T> // T must be floating point type
  inline T cubic_function_dfdx(T x, T a, T b, T c) {
    T x2 = x * x;
    return x2 * a * 3 + x * b * 2 + c;
  }

  template <typename T> // T must be floating point type
  inline T cubic_function_d2fd2x(T x, T a, T b) {
    return x * a * 6 + b * 2;
  }

  /**
   * Iterative approach. https://en.wikipedia.org/wiki/Halley%27s_method
   *  Halley's method
   */

  inline std::vector<double> solve_cubic_real(long double a, long double b, long double c,
                                              long double d, double const epsilon) {

    CORSIKA_LOG_TRACE("cubic_iterative: a={:f}, b={:f}, c={:f}, d={:f}, epsilon={} {} {}",
                      a, b, c, d, epsilon, (std::abs(a - 1) < epsilon),
                      (std::abs(b) < epsilon));

    if (std::abs(a) < epsilon) { // this is just a quadratic
      return solve_quadratic_real(b, c, d, epsilon);
    }

    auto pre_opt = andre::solve_cubic_real_analytic(a, b, c, d, epsilon);
    long double x1 = 0; // start value

    if (pre_opt.size()) {
      x1 = pre_opt[0]; //*std::max_element(pre_opt.begin(), pre_opt.end());
#ifdef _C8_DEBUG_
      for (long double test_v : pre_opt) {
        CORSIKA_LOG_TRACE("test,andre x={} f(x)={}", test_v,
                          cubic_function(test_v, a, b, c, d));
      }
#endif
    } else {
      // this is only if the former solve_cubic_real_analytic would not result
      // in any solution. We have no test case for this. This is excluded from tests:
      // LCOV_EXCL_START
      long double const dist = std::fma(b / a, b / a, -3 * c / a);
      long double const xinfl = -b / (a * 3);

      x1 = xinfl;
      long double f_test = cubic_function(xinfl, a, b, c, d);

      if (std::abs(f_test) > epsilon) {
        if (std::abs(dist) < epsilon) {
          x1 = xinfl - std::cbrt(f_test);
        } else if (dist > 0) {
          if (f_test > 0)
            x1 = xinfl - 2 / 3 * std::sqrt(dist);
          else
            x1 = xinfl + 2 / 3 * std::sqrt(dist);
        }
      }
      // LCOV_EXCL_STOP
    }

    long double f_x1 = cubic_function(x1, a, b, c, d);
    long double dx1 = 0;

    int niter = 0;
    const int maxiter = 100;
    do {
      long double const f_prime_x1 = cubic_function_dfdx(x1, a, b, c);
      long double const f_prime2_x1 = cubic_function_d2fd2x(x1, a, b);
      // if (potential) saddle point... avoid
      if (std::abs(f_prime_x1) < epsilon) {
        dx1 = std::cbrt(f_x1);
      } else {
        dx1 = f_x1 * f_prime_x1 * 2 / (f_prime_x1 * f_prime_x1 * 2 - f_x1 * f_prime2_x1);
      }
      x1 -= dx1;
      f_x1 = cubic_function(x1, a, b, c, d);
      CORSIKA_LOG_TRACE(
          "niter={} x1={:.20f} f_x1={:.20f} f_prime={:.20f} f_prime2={:.20f} dx1={}",
          niter, x1, f_x1, f_prime_x1, f_prime2_x1,
          f_x1 * f_prime_x1 / (f_prime_x1 * f_prime_x1 - f_x1 * f_prime2_x1 * 0.5));
    } while ((++niter < maxiter) && (std::abs(f_x1) > epsilon * 1000) &&
             (std::abs(dx1) > epsilon));

    CORSIKA_LOG_TRACE("niter={}", niter);
    if (niter >= maxiter) {
      CORSIKA_LOG_DEBUG("niter reached max iterations {}", niter);
      return andre::solve_cubic_real_analytic(a, b, c, d, epsilon);
    }

    CORSIKA_LOG_TRACE("x1={} f_x1={}", x1, f_x1);

    double const b1 = x1 + b / a;
    double const b0 = b1 * x1 + c / a;
    std::vector<double> quad_check = solve_quadratic_real(1, b1, b0, 1e-3);
    CORSIKA_LOG_TRACE("quad_check=[{}], f(z)={}", fmt::join(quad_check, ", "),
                      cubic_function(x1, a, b, c, d));

    quad_check.push_back(x1);
    CORSIKA_LOG_TRACE("cubic: solve_cubic_real returns={}", fmt::join(quad_check, ", "));
    return quad_check;
  } // namespace corsika

} // namespace corsika
