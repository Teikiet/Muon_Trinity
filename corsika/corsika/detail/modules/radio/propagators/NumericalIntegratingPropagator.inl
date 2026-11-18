/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/propagators/NumericalIntegratingPropagator.hpp>

namespace corsika {

  template <typename TEnvironment>
  // TODO: maybe the constructor doesn't take any arguments for the environment (?)
  inline NumericalIntegratingPropagator<TEnvironment>::NumericalIntegratingPropagator(
      TEnvironment const& env, LengthType const stepsize)
      : RadioPropagator<NumericalIntegratingPropagator, TEnvironment>(env)
      , stepsize_(stepsize) {}

  template <typename TEnvironment>
  template <typename Particle>
  inline typename NumericalIntegratingPropagator<TEnvironment>::SignalPathCollection
  NumericalIntegratingPropagator<TEnvironment>::propagate(
      Particle const& particle, Point const& source, Point const& destination) const {

    /*
     * get the normalized (unit) vector from `source` to `destination'.
     * this is also the `emit` and `receive` vectors in the SignalPath class.
     * in this case emit and receive unit vectors should be the same
     * so they are both called direction
     */

    // these are used for the direction of emission and reception of signal at the
    // observer
    auto const emit{(destination - source).normalized()};
    auto const receive{-emit};

    // the distance from the point of emission to an observer
    auto const distance{(destination - source).getNorm()};

    try {
      if (stepsize_ <= 0.5 * distance) {

        // "step" is the direction vector with length `stepsize`
        auto const step{emit * stepsize_};

        // calculate the number of points (roughly) for the numerical integration
        auto const n_points{(destination - source).getNorm() / stepsize_};

        // get the universe for this environment
        auto const* const universe{Base::env_.getUniverse().get()};

        // the points that we build along the way for the numerical integration
        std::deque<Point> points;

        // store value of the refractive index at points
        std::vector<double> rindex;
        rindex.reserve(n_points);

        // get and store the refractive index of the first point 'source'
        auto const* const nodeSource{particle.getNode()};
        auto const ri_source{nodeSource->getModelProperties().getRefractiveIndex(source)};
        rindex.push_back(ri_source);
        points.push_back(source);

        // loop from `source` to `destination` to store values before Simpson's rule.
        // this loop skips the last point 'destination' and "misses" the extra point
        for (auto point = source + step;
             (point - destination).getNorm() > 0.6 * stepsize_; point = point + step) {

          // get the environment node at this specific 'point'
          auto const* const node{universe->getContainingNode(point)};

          // get the associated refractivity at 'point'
          auto const refractive_index{
              node->getModelProperties().getRefractiveIndex(point)};
          rindex.push_back(refractive_index);

          // add this 'point' to our deque collection
          points.push_back(point);
        }

        // Get the extra points that the for loop misses until the destination
        auto const extrapoint_{points.back() + step};

        // add the refractive index of last point 'destination' and store it
        auto const* const node{universe->getContainingNode(destination)};
        auto const ri_destination{
            node->getModelProperties().getRefractiveIndex(destination)};
        rindex.push_back(ri_destination);
        points.push_back(destination);

        auto N = rindex.size();
        std::size_t index = 0;
        double sum = rindex.at(index);
        auto refra = rindex.at(index);
        TimeType time{0_s};

        if ((N - 1) % 2 == 0) {
          // Apply the standard Simpson's rule
          auto const h = ((destination - source).getNorm()) / (N - 1);

          for (std::size_t index = 1; index < (N - 1); index += 2) {
            sum += 4 * rindex.at(index);
            refra += rindex.at(index);
          }
          for (std::size_t index = 2; index < (N - 1); index += 2) {
            sum += 2 * rindex.at(index);
            refra += rindex.at(index);
          }
          index = N - 1;
          sum = sum + rindex.at(index);
          refra += rindex.at(index);

          // compute the total time delay.
          time = sum * (h / (3 * constants::c));
        } else {
          // Apply Simpson's rule for one "extra" point and then subtract the difference
          points.pop_back();
          rindex.pop_back();
          auto const* const node{universe->getContainingNode(extrapoint_)};
          auto const ri_extrapoint{
              node->getModelProperties().getRefractiveIndex(extrapoint_)};
          rindex.push_back(ri_extrapoint);
          points.push_back(extrapoint_);
          auto const extrapoint2_{extrapoint_ + step};
          auto const* const node2{universe->getContainingNode(extrapoint2_)};
          auto const ri_extrapoint2{
              node2->getModelProperties().getRefractiveIndex(extrapoint2_)};
          rindex.push_back(ri_extrapoint2);
          points.push_back(extrapoint2_);
          N = rindex.size();
          auto const h = ((extrapoint2_ - source).getNorm()) / (N - 1);
          for (std::size_t index = 1; index < (N - 1); index += 2) {
            sum += 4 * rindex.at(index);
            refra += rindex.at(index);
          }
          for (std::size_t index = 2; index < (N - 1); index += 2) {
            sum += 2 * rindex.at(index);
            refra += rindex.at(index);
          }
          index = N - 1;
          sum = sum + rindex.at(index);
          refra += rindex.at(index);

          // compute the total time delay including the correction
          time =
              sum * (h / (3 * constants::c)) -
              (ri_extrapoint2 * ((extrapoint2_ - destination).getNorm()) / constants::c);
        }

        // compute the average refractive index.
        auto const averageRefractiveIndex = refra / N;

        return std::vector<SignalPath>(
            1, SignalPath(time, averageRefractiveIndex, ri_source, ri_destination, emit,
                          receive, distance, points));
      } else {
        throw stepsize_;
      }
    } catch (const LengthType& s) {
      CORSIKA_LOG_ERROR("Please choose a smaller stepsize for the numerical integration");
    }

    std::deque<Point> const defaultPoints;
    return std::vector<SignalPath>(
        1, SignalPath(0_s, 0, 0, 0, emit, receive, distance,
                      defaultPoints)); // Dummy return that is never called (combat
                                       // compile warnings)
  }                                    // END: propagate()

} // namespace corsika
