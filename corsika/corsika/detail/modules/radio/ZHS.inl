/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/radio/ZHS.hpp>

namespace corsika {

  template <typename TRadioDetector, typename TPropagator>
  template <typename Particle>
  inline ProcessReturn ZHS<TRadioDetector, TPropagator>::simulate(
      Step<Particle> const& step) {
    auto const startTime{step.getTimePre()};
    auto const endTime{step.getTimePost()};

    // LCOV_EXCL_START
    if (startTime == endTime) {
      return ProcessReturn::Ok;
      // LCOV_EXCL_STOP
    } else {

      auto const startPoint{step.getPositionPre()};
      auto const endPoint{step.getPositionPost()};
      LengthType const trackLength{(startPoint - endPoint).getNorm()};

      auto const betaModule{(endPoint - startPoint).getNorm() /
                            (constants::c * (endTime - startTime))};
      auto const beta{(endPoint - startPoint).normalized() * betaModule};

      auto const charge{get_charge(step.getParticlePre().getPID())};

      // // get "mid" position of the track geometrically

      auto const halfVector{(startPoint - endPoint) / 2};
      auto const midPoint{endPoint + halfVector};

      // get thinning weight
      auto const thinningWeight{step.getParticlePre().getWeight()};

      auto const constants{charge * emConstant_ * thinningWeight};

      // we loop over each observer in the collection
      for (auto& observer : observers_.getObservers()) {
        auto midPaths{this->propagator_.propagate(step.getParticlePre(), midPoint,
                                                  observer.getLocation())};
        // Loop over midPaths, first check Fraunhoffer limit
        for (size_t i{0}; i < midPaths.size(); i++) {
          double const uTimesK{beta.dot(midPaths[i].emit_) / betaModule};
          double const sinTheta2{1. - uTimesK * uTimesK};
          LengthType const lambda{constants::c / observer.getSampleRate()};
          double const fraunhLimit{sinTheta2 * trackLength * trackLength /
                                   midPaths[i].R_distance_ / lambda * 2 * M_PI};
          // Checks if we are in fraunhoffer domain
          if (fraunhLimit > 1.0) {
            /// code for dividing track and calculating field.
            int const nSubTracks = sqrt(fraunhLimit) + 1.0;
            auto const step_{(endPoint - startPoint) / double(nSubTracks)};
            TimeType const timeStep{(endTime - startTime) / double(nSubTracks)};
            // energy should be divided up when it is possible to get the energy at end of
            // track!!!!
            auto point1{startPoint};
            TimeType time1{startTime};
            for (int j{0}; j < nSubTracks; j++) {
              auto const point2{point1 + step_};
              TimeType const time2{time1 + timeStep};
              auto const newHalfVector{(point1 - point2) / 2.};
              auto const newMidPoint{point2 + newHalfVector};
              auto const newMidPaths{this->propagator_.propagate(
                  step.getParticlePre(), newMidPoint, observer.getLocation())};
              // A function for calculating the field should be made since it is repeated
              // later
              for (size_t k{0}; k < newMidPaths.size(); k++) {

                double const n_source{newMidPaths[k].refractive_index_source_};

                double const betaTimesK{beta.dot(newMidPaths[k].emit_)};
                TimeType const midTime{(time1 + time2) / 2.};
                TimeType detectionTime1{time1 + newMidPaths[k].propagation_time_ -
                                        n_source * betaTimesK * (time1 - midTime)};
                TimeType detectionTime2{time2 + newMidPaths[k].propagation_time_ -
                                        n_source * betaTimesK * (time2 - midTime)};

                // make detectionTime1_ be the smallest time =>
                // changes step function order so sign is changed to account for it
                double sign = 1.;
                if (detectionTime1 > detectionTime2) {
                  detectionTime1 = time2 + newMidPaths[k].propagation_time_ -
                                   n_source * betaTimesK * (time2 - midTime);
                  detectionTime2 = time1 + newMidPaths[k].propagation_time_ -
                                   n_source * betaTimesK * (time1 - midTime);
                  sign = -1.;
                } // end if statement for time structure

                double const startBin{
                    std::floor((detectionTime1 - observer.getStartTime()) *
                                   observer.getSampleRate() +
                               0.5)};
                double const endBin{
                    std::floor((detectionTime2 - observer.getStartTime()) *
                                   observer.getSampleRate() +
                               0.5)};

                auto const betaPerp{
                    newMidPaths[k].emit_.cross(beta.cross(newMidPaths[k].emit_))};
                double const denominator{1. - n_source * betaTimesK};

                if (startBin == endBin) {
                  // track contained in bin
                  // if not in Cerenkov angle then
                  if (std::fabs(denominator) > 1.e-15) {
                    double const f{
                        std::fabs((detectionTime2 * observer.getSampleRate() -
                                   detectionTime1 * observer.getSampleRate()))};
                    VectorPotential const Vp = betaPerp * sign * constants * f /
                                               denominator / newMidPaths[k].R_distance_;
                    observer.receive(detectionTime2, newMidPaths[k], Vp);
                  } else { // If emission in Cerenkov angle => approximation
                    double const f{time2 * observer.getSampleRate() -
                                   time1 * observer.getSampleRate()};
                    VectorPotential const Vp =
                        betaPerp * sign * constants * f / newMidPaths[k].R_distance_;
                    observer.receive(detectionTime2, newMidPaths[k], Vp);
                  } // end if Cerenkov angle approx
                } else {
                  /*Track is contained in more than one bin*/
                  int const numberOfBins{static_cast<int>(endBin - startBin)};
                  // first contribution/ plus 1 bin minus 0.5 from new observer ruonding
                  double f{std::fabs(startBin + 0.5 -
                                     (detectionTime1 - observer.getStartTime()) *
                                         observer.getSampleRate())};
                  VectorPotential Vp = betaPerp * sign * constants * f / denominator /
                                       newMidPaths[k].R_distance_;
                  observer.receive(detectionTime1, newMidPaths[k], Vp);
                  // intermidiate contributions
                  for (int it{1}; it < numberOfBins; ++it) {
                    Vp = betaPerp * constants / denominator / newMidPaths[k].R_distance_;
                    observer.receive(detectionTime1 + static_cast<double>(it) /
                                                          observer.getSampleRate(),
                                     newMidPaths[k], Vp);
                  } // end loop over bins in which potential vector is not zero
                  // final contribution// f +0.5 from new observer rounding
                  f = std::fabs((detectionTime2 - observer.getStartTime()) *
                                    observer.getSampleRate() +
                                0.5 - endBin);
                  Vp = betaPerp * sign * constants * f / denominator /
                       newMidPaths[k].R_distance_;
                  observer.receive(detectionTime2, newMidPaths[k], Vp);
                } // end if statement for track in multiple bins

              } // end of loop over newMidPaths
              // update points for next sub track
              point1 = point1 + step_;
              time1 = time1 + timeStep;
            }

          } else // Calculate vector potential of whole track
          {
            double const n_source{midPaths[i].refractive_index_source_};

            double const betaTimesK{beta.dot(midPaths[i].emit_)};
            TimeType const midTime{(startTime + endTime) / 2};
            TimeType detectionTime1{startTime + midPaths[i].propagation_time_ -
                                    n_source * betaTimesK * (startTime - midTime)};
            TimeType detectionTime2{endTime + midPaths[i].propagation_time_ -
                                    n_source * betaTimesK * (endTime - midTime)};

            // make detectionTime1_ be the smallest time =>
            // changes step function order so sign is changed to account for it
            double sign = 1.;
            if (detectionTime1 > detectionTime2) {
              detectionTime1 = endTime + midPaths[i].propagation_time_ -
                               n_source * betaTimesK * (endTime - midTime);
              detectionTime2 = startTime + midPaths[i].propagation_time_ -
                               n_source * betaTimesK * (startTime - midTime);
              sign = -1.;
            } // end if statement for time structure

            double const startBin{std::floor((detectionTime1 - observer.getStartTime()) *
                                                 observer.getSampleRate() +
                                             0.5)};
            double const endBin{std::floor((detectionTime2 - observer.getStartTime()) *
                                               observer.getSampleRate() +
                                           0.5)};

            auto const betaPerp{midPaths[i].emit_.cross(beta.cross(midPaths[i].emit_))};
            double const denominator{1. -
                                     midPaths[i].refractive_index_source_ * betaTimesK};

            if (startBin == endBin) {
              // track contained in bin
              // if not in Cerenkov angle then
              if (std::fabs(denominator) > 1.e-15) {
                double const f{std::fabs((detectionTime2 * observer.getSampleRate() -
                                          detectionTime1 * observer.getSampleRate()))};

                VectorPotential const Vp = betaPerp * sign * constants * f / denominator /
                                           midPaths[i].R_distance_;
                observer.receive(detectionTime2, midPaths[i], Vp);
              } else { // If emission in Cerenkov angle => approximation
                double const f{endTime * observer.getSampleRate() -
                               startTime * observer.getSampleRate()};
                VectorPotential const Vp =
                    betaPerp * sign * constants * f / midPaths[i].R_distance_;
                observer.receive(detectionTime2, midPaths[i], Vp);
              } // end if Cerenkov angle approx
            } else {
              /*Track is contained in more than one bin*/
              int const numberOfBins{static_cast<int>(endBin - startBin)};
              // TODO: should we check for Cerenkov angle?
              // first contribution
              double f{std::fabs(startBin + 0.5 -
                                 (detectionTime1 - observer.getStartTime()) *
                                     observer.getSampleRate())};
              VectorPotential Vp =
                  betaPerp * sign * constants * f / denominator / midPaths[i].R_distance_;
              observer.receive(detectionTime1, midPaths[i], Vp);
              // intermediate contributions
              for (int it{1}; it < numberOfBins; ++it) {
                Vp = betaPerp * sign * constants / denominator / midPaths[i].R_distance_;
                observer.receive(
                    detectionTime1 + static_cast<double>(it) / observer.getSampleRate(),
                    midPaths[i], Vp);
              } // end loop over bins in which potential vector is not zero
              // final contribution
              f = std::fabs((detectionTime2 - observer.getStartTime()) *
                                observer.getSampleRate() +
                            0.5 - endBin);
              Vp =
                  betaPerp * sign * constants * f / denominator / midPaths[i].R_distance_;
              observer.receive(detectionTime2, midPaths[i], Vp);
            } // end if statement for track in multiple bins

          } // finish if statement of track in fraunhoffer or not

        } // end loop over mid paths

      } // END: loop over observers
      return ProcessReturn::Ok;
    }
  } // end simulate

} // namespace corsika
