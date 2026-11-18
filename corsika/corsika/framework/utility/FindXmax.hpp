/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>

#include <Eigen/Dense>

#include <cassert>
#include <cmath>
#include <optional>
#include <tuple>
#include <vector>

namespace corsika {

  /**
   * Interpolates profiles around maximum.
   *
   * The maximum of profiles can be robustly estimated from the three maximal points by
   * quadratic interpolation. This code is copied from CONEX and is awaiting full adaption
   * into CORSIKA8. This is a temporary solution (just copied).
   */
  namespace find_xmax {

    namespace detail {

      using Matrix3RowMajor = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

      // Solves A * x = y for 3x3 A. Returns false if A is singular or numerically
      // non-invertible.
      [[nodiscard]] inline bool solve_3by3(double const y[3], double A[3][3],
                                           double x[3]) {
        Eigen::Map<Matrix3RowMajor> mat(A[0]);
        Eigen::Map<Eigen::Vector3d const> rhs(y);
        Eigen::Map<Eigen::Vector3d> sol(x);

        Eigen::FullPivLU<Matrix3RowMajor> lu(mat);
        if (!lu.isInvertible()) {
          CORSIKA_LOG_WARN(
              "solve_3by3: matrix is singular or numerically non-invertible");
          return false;
        }

        sol = lu.solve(rhs);
        return true;
      }

    } // namespace detail

    [[nodiscard]] inline std::optional<std::tuple<double, double>> interpolate_profile(
        double x1, double x2, double x3, double y1, double y2, double y3) {

      // quadratic "fit" around maximum to get dEdXmax and Xmax
      double x[3] = {x1, x2, x3};
      double y[3] = {y1, y2, y3};

      double A[3][3];
      A[0][0] = x[0] * x[0];
      A[0][1] = x[0];
      A[0][2] = 1.;
      A[1][0] = x[1] * x[1];
      A[1][1] = x[1];
      A[1][2] = 1.;
      A[2][0] = x[2] * x[2];
      A[2][1] = x[2];
      A[2][2] = 1.;

      double a[3];

      if (!detail::solve_3by3(y, A, a)) { return std::nullopt; }

      if (a[0] < 0.) {
        double const Xmax = -a[1] / (2. * a[0]);
        return std::make_tuple(Xmax, a[0] * Xmax * Xmax + a[1] * Xmax + a[2]);
      }
      return std::nullopt;
    }

    [[nodiscard]] inline std::optional<std::tuple<double, double>>
    estimate_parabola_parameters(std::vector<double> const& x,
                                 std::vector<double> const& y) {
      assert(x.size() == y.size());
      assert(x.size() >= 3);
      size_t const N = x.size();

      double x0 = 0;
      for (size_t i = 0; i < N; i++) { x0 += x[i]; }
      x0 /= static_cast<double>(N);

      // accumulate sums in mean-centered coordinates (dx = x - x0)
      double Sx2 = 0;
      double Sx3 = 0;
      double Sx4 = 0;
      double Sy = 0;
      double Sxy = 0;
      double Sx2y = 0;
      for (size_t i = 0; i < N; i++) {
        double const dx = x[i] - x0;
        double const dx2 = dx * dx;
        double const yi = y[i];
        Sx2 += dx2;
        Sx3 += dx2 * dx;
        Sx4 += dx2 * dx2;
        Sy += yi;
        Sxy += dx * yi;
        Sx2y += dx2 * yi;
      }

      // solve normal equations for y = a*dx^2 + b'*dx + c  (dx = x - x0)
      // | Sx4 Sx3 Sx2 |   |a |   |Sx2y|
      // | Sx3 Sx2  0  | * |b'| = |Sxy |   (Sx' = 0 at mean)
      // | Sx2  0   N  |   |c |   |Sy  |

      double A[3][3] = {
          {Sx4, Sx3, Sx2}, {Sx3, Sx2, 0.}, {Sx2, 0., static_cast<double>(N)}};
      double const rhs[3] = {Sx2y, Sxy, Sy};
      double coef[3];
      if (!detail::solve_3by3(rhs, A, coef)) { return std::nullopt; }

      double const a = coef[0];
      double const b = coef[1];
      double const c = coef[2];

      if (a < 0) {
        double const Xmax = x0 - b / (2. * a);
        double const dxMax = Xmax - x0;
        double const Ymax = a * dxMax * dxMax + b * dxMax + c;
        return std::make_tuple(Xmax, Ymax);
      }
      return std::nullopt;
    }

  } // namespace find_xmax
} // namespace corsika
