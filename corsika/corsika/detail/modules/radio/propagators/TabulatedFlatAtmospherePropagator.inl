/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/propagators/TabulatedFlatAtmospherePropagator.hpp>

namespace corsika {

  template <typename TEnvironment>
  inline TabulatedFlatAtmospherePropagator<
      TEnvironment>::TabulatedFlatAtmospherePropagator(TEnvironment const& env,
                                                       Point const& upperLimit,
                                                       Point const& lowerLimit,
                                                       LengthType const step)
      : RadioPropagator<TabulatedFlatAtmospherePropagator, TEnvironment>(env)
      , upperLimit_(upperLimit)
      , lowerLimit_(lowerLimit)
      , step_(step)
      , inverseStep_(1 / step)
      , maxHeight_(upperLimit_.getCoordinates().getZ())
      , minHeight_(lowerLimit_.getCoordinates().getZ() - 1_km) {
    auto const minX_ = lowerLimit_.getCoordinates().getX();
    auto const minY_ = lowerLimit_.getCoordinates().getY();
    std::size_t const nBins_ = (maxHeight_ - minHeight_) * inverseStep_ + 1;
    refractivityTable_.reserve(nBins_);
    heightTable_.reserve(nBins_);
    integratedRefractivityTable_.reserve(nBins_);

    // get the root coordinate system of this environment
    auto const& rootCS = env.getCoordinateSystem();
    // get the universe for this environment
    auto const* const universe{Base::env_.getUniverse().get()};

    for (std::size_t i = 0; i < nBins_; i++) {
      Point point_{rootCS, minX_, minY_, minHeight_ + i * step_};
      auto const* const node{universe->getContainingNode(point_)};
      auto const ri_ = node->getModelProperties().getRefractiveIndex(point_);
      refractivityTable_.push_back(ri_ - 1);
      auto const height_ = minHeight_ + i * step_;
      heightTable_.push_back(height_);
    }

    double const stepOverMeter_{inverseStep_ * 1_m};
    auto const intRerfZero_ = refractivityTable_.at(0) * stepOverMeter_;
    integratedRefractivityTable_.push_back(intRerfZero_);
    for (std::size_t i = 1; i < nBins_; i++) {
      auto const intRefrI_ =
          integratedRefractivityTable_[i - 1] + refractivityTable_[i] * stepOverMeter_;
      integratedRefractivityTable_.push_back(intRefrI_);
    }

    // this is used to interpolate refractivity and integrated refractivity for particles
    // below the lower limit
    slopeRefrLower_ = (refractivityTable_.at(10) - refractivityTable_.front()) / 10.;
    slopeIntRefrLower_ =
        (integratedRefractivityTable_.at(10) - integratedRefractivityTable_.front()) /
        10.;

    // this is used to interpolate refractivity and integrated refractivity for particles
    // above the upper limit
    lastElement_ = integratedRefractivityTable_.size() - 1;
    slopeRefrUpper_ =
        (refractivityTable_.at(lastElement_) - refractivityTable_.at(lastElement_ - 10)) /
        10.;
    slopeIntRefrUpper_ = (integratedRefractivityTable_.at(lastElement_) -
                          integratedRefractivityTable_.at(lastElement_ - 10)) /
                         10.;
  }

  template <typename TEnvironment>
  template <typename Particle>
  inline typename TabulatedFlatAtmospherePropagator<TEnvironment>::SignalPathCollection
  TabulatedFlatAtmospherePropagator<TEnvironment>::propagate(
      [[maybe_unused]] Particle const& particle, Point const& source,
      Point const& destination) {

    /**
     * This is a simple case of straight propagator where
     * tabulated values of refractive index are called assuming
     * a flat atmosphere.
     *
     */

    // these are used for the direction of emission and reception of signal at the
    // observer
    auto const emit_{(destination - source).normalized()};
    auto const receive_{-emit_};

    // the geometrical distance from the point of emission to an observer
    auto const distance_{(destination - source).getNorm()};

    // clear the refractive index vector and points deque for this signal propagation.
    rindex.clear();
    points.clear();
    auto const sourceZ_{source.getCoordinates().getZ()};
    auto const checkMaxZ_{(maxHeight_ - sourceZ_) / 1_m};
    auto ri_source{1.};
    auto intRerfr_source{1.};
    auto ri_destination{1.};
    auto intRerfr_destination{1.};
    auto height_{1.};

    double const sourceHeight_{(sourceZ_ - heightTable_.front()) * inverseStep_};
    double const destinationHeight_{
        (destination.getCoordinates().getZ() - heightTable_.front()) * inverseStep_};

    if ((sourceHeight_ + 0.5) >=
        lastElement_) { // source particle is above maximum height
      ri_source =
          refractivityTable_.back() +
          slopeRefrUpper_ * std::abs((sourceZ_ - heightTable_.back()) * inverseStep_) + 1;
      intRerfr_source =
          integratedRefractivityTable_.back() + slopeIntRefrUpper_ * std::abs(checkMaxZ_);
      rindex.push_back(ri_source);
      points.push_back(source);

      // add the refractive index of last point 'destination' and store it.
      std::size_t const indexDestination_{
          static_cast<std::size_t>(destinationHeight_ + 0.5)};
      ri_destination = refractivityTable_.at(indexDestination_) + 1;
      intRerfr_destination = integratedRefractivityTable_.at(indexDestination_);
      rindex.push_back(ri_destination);
      points.push_back(destination);

      height_ = sourceHeight_ - destinationHeight_;

    } else if ((sourceHeight_ + 0.5 < lastElement_) &&
               sourceHeight_ > 0) { // source particle in the table
      // get and store the refractive index of the first point 'source'.
      std::size_t const indexSource_{static_cast<std::size_t>(sourceHeight_ + 0.5)};
      ri_source = refractivityTable_.at(indexSource_) + 1;
      intRerfr_source = integratedRefractivityTable_.at(indexSource_);
      rindex.push_back(ri_source);
      points.push_back(source);

      // add the refractive index of last point 'destination' and store it.
      std::size_t const indexDestination_{
          static_cast<std::size_t>(destinationHeight_ + 0.5)};
      ri_destination = refractivityTable_.at(indexDestination_) + 1;
      intRerfr_destination = integratedRefractivityTable_.at(indexDestination_);
      rindex.push_back(ri_destination);
      points.push_back(destination);

      height_ =
          (heightTable_.at(indexSource_) - heightTable_.at(indexDestination_)) / 1_m;
      if (height_ == 0) { height_ = 1.; }
    } else if (sourceHeight_ ==
               0) { // source particle is exactly at the lower edge of the table
      // get and store the refractive index of the first point 'source'.
      std::size_t const indexSource_{static_cast<std::size_t>(sourceHeight_ + 0.5)};
      auto const ri_source{refractivityTable_.at(indexSource_) + 1};
      intRerfr_source = integratedRefractivityTable_.at(indexSource_);
      rindex.push_back(ri_source);
      points.push_back(source);

      // add the refractive index of last point 'destination' and store it.
      std::size_t const indexDestination_{
          static_cast<std::size_t>(destinationHeight_ + 0.5)};
      auto const ri_destination{refractivityTable_.at(indexDestination_) + 1};
      intRerfr_destination = integratedRefractivityTable_.at(indexDestination_);
      rindex.push_back(ri_destination);
      points.push_back(destination);

      height_ = destinationHeight_ - sourceHeight_;
    } else if (sourceHeight_ < 0) { // source particle is in the ground.
      // get and store the refractive index of the first point 'source'.
      ri_source =
          refractivityTable_.front() + slopeRefrLower_ * std::abs(sourceHeight_) + 1;
      intRerfr_source = integratedRefractivityTable_.front() +
                        slopeIntRefrLower_ * std::abs(sourceHeight_);
      rindex.push_back(ri_source);
      points.push_back(source);

      // add the refractive index of last point 'destination' and store it.
      std::size_t const indexDestination_{
          static_cast<std::size_t>(destinationHeight_ + 0.5)};
      auto const ri_destination{refractivityTable_.at(indexDestination_) + 1};
      intRerfr_destination = integratedRefractivityTable_.at(indexDestination_);
      rindex.push_back(ri_destination);
      points.push_back(destination);

      height_ = destinationHeight_ - std::abs(sourceHeight_);
    }

    // compute the average refractive index.
    auto const averageRefractiveIndex_ = (ri_source + ri_destination) * 0.5;

    // compute the total time delay.
    TimeType const time = (1 + (intRerfr_source - intRerfr_destination) / height_) *
                          (distance_ / constants::c);

    return {SignalPath(time, averageRefractiveIndex_, ri_source, ri_destination, emit_,
                       receive_, distance_, points)};

  } // END: propagate()

} // namespace corsika
