/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/media/density_and_composition/FlatInterpolatedMedium.hpp>

#include <boost/math/tools/minima.hpp>

#include <string>

namespace corsika {
  namespace media {

    template <typename TEnvModel>
    inline ShowerAxis::ShowerAxis(Point const& pStart, Point const& pEnd,
                                  media::Environment<TEnvModel> const& env,
                                  bool const doThrow, int const steps)
        : ShowerAxis(pStart, pEnd - pStart, env, doThrow, steps) {}

    template <typename TEnvModel>
    inline ShowerAxis::ShowerAxis(Point const& pStart, Vector<length_d> const& length,
                                  media::Environment<TEnvModel> const& env,
                                  bool const doThrow, int const steps)
        : pointStart_(pStart)
        , length_(length)
        , throw_(doThrow)
        , max_length_(length_.getNorm())
        , steplength_(max_length_ / steps)
        , axis_normalized_(length / max_length_)
        , X_(steps + 1) {

      auto const* const universe = env.getUniverse().get();

      // lambda to get the density at a fractional distance, x, along axis
      auto rho = [&](double x) {
        auto const p = pointStart_ + length_ * x;
        auto const* node = universe->getContainingNode(p);
        if (!node->hasModelProperties()) {
          CORSIKA_LOG_CRITICAL(
              "Unable to construct ShowerAxis. media::ShowerAxis includes volume "
              "with no model properties at point {}.",
              p);
          throw std::runtime_error("Unable to construct ShowerAxis. Missing properties");
        }
        return node->getModelProperties().getMassDensity(p).magnitude();
      };

      double error;
      X_[0] = GrammageType::zero();
      auto sum = GrammageType::zero();

      // note: the default calculation of shower axis bins using gauss_kronrod integration
      // is very slow for a FlatInterpolatedMedium, therefore an alternative calculation
      // is implemented

      // get top-level node in universe and try to cast its medium
      // to an FlatInterpolatedMedium
      auto const* node = universe->getContainingNode(pointStart_);
      auto const* interpMediumPtr =
          dynamic_cast<FlatInterpolatedMedium<TEnvModel> const*>(
              &node->getModelProperties());

      // check universe and its containing nodes
      if (universe->getChildNodes().size() == 1 &&
          universe->getChildNodes().front()->getChildNodes().size() == 0 &&
          interpMediumPtr != nullptr) {
        // universe contains a single interpolated medium, use the alternative calculation

        // initialize density
        MassDensityType densityStart = rho(0) * 1_kg / static_pow<3>(1_m);
        MassDensityType densityEnd;

        for (int i = 1; i <= steps; ++i) {
          auto const frac = static_cast<double>(i) / steps;

          // get step endpoint density
          densityEnd = rho(frac) * 1_kg / static_pow<3>(1_m);

          // calculate integrated density
          auto const result = steplength_ * (densityStart + densityEnd) / 2;

          // grammage in a step cannot be negative
          assert(result > GrammageType::zero());

          sum += result;
          X_[i] = sum;

          // current endpoint is next step's startpoint
          densityStart = densityEnd;
        }
      } else {
        // standard atmosphere with multiple volumes
        for (int i = 1; i <= steps; ++i) {
          auto const frac_prev = (i - 1.0) / steps;
          auto const frac = double(i) / steps;
          auto const r = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
              rho, frac_prev, frac, 15, 1e-9, &error);
          auto const result =
              MassDensityType(phys::units::detail::magnitude_tag, r) * max_length_;
          sum += result;
          X_[i] = sum;
        }
      }

      if (!std::is_sorted(X_.cbegin(), X_.cend())) {
        CORSIKA_LOG_ERROR("Shower axis bins are not sorted");
        throw std::runtime_error("Bad media::ShowerAxis construction");
      }
    }

    inline GrammageType ShowerAxis::getX(LengthType l) const {
      double const fractionalBin = l / steplength_;
      int const lower = std::floor(fractionalBin); // indices of nearest X support points
      double const fraction = fractionalBin - lower;
      unsigned int const upper = lower + 1;

      if (fractionalBin < 0) {
        CORSIKA_LOG_TRACE("cannot extrapolate to points behind point of injection l={} m",
                          l / 1_m);
        if (throw_) {
          throw std::runtime_error(
              "cannot extrapolate to points behind point of injection");
        }
        return getMinimumX();
      }

      if (upper >= X_.size()) {
        CORSIKA_LOG_TRACE(
            "shower axis too short, cannot extrapolate (l / max_length_ = {} )",
            l / max_length_);
        if (throw_) {
          const std::string err = fmt::format(
              "shower axis too short, cannot extrapolate (l / max_length_ = {} )",
              l / max_length_);
          throw std::runtime_error(err.c_str());
        }
        return getMaximumX();
      }
      CORSIKA_LOG_TRACE("showerAxis::X frac={}, fractionalBin={}, lower={}, upper={}",
                        fraction, fractionalBin, lower, upper);

      if (!(0 <= fraction && fraction <= 1.)) {
        CORSIKA_LOG_ERROR("Fraction must be 0 <= {} <= 1", fraction);
        throw std::runtime_error("Fraction out of range");
      }

      CORSIKA_LOG_TRACE("ShowerAxis::getX l={} m, lower={}, fraction={}, upper={}",
                        l / 1_m, lower, fraction, upper);

      // linear interpolation between getX[lower] and X[upper]
      return X_[upper] * fraction + X_[lower] * (1 - fraction);
    }

    inline LengthType ShowerAxis::getSteplength() const { return steplength_; }

    inline GrammageType ShowerAxis::getMaximumX() const { return *X_.rbegin(); }

    inline GrammageType ShowerAxis::getMinimumX() const { return GrammageType::zero(); }

    inline GrammageType ShowerAxis::getProjectedX(Point const& p) const {
      auto const projectedLength = (p - pointStart_).dot(axis_normalized_);
      return getX(projectedLength);
    }

    inline DirectionVector const& ShowerAxis::getDirection() const {
      return axis_normalized_;
    }

    inline Point const& ShowerAxis::getStart() const { return pointStart_; }

  } // namespace media
} // namespace corsika
