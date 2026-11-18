/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/radio/CoREAS.hpp>

namespace corsika {

  template <typename TRadioDetector, typename TPropagator>
  template <typename Particle>
  inline ProcessReturn CoREAS<TRadioDetector, TPropagator>::simulate(
      Step<Particle> const& step) {

    // get the global simulation time for that track.
    auto const startTime_{
        step.getTimePre()}; // time at the start point of the track. I
                            // should use something similar to fCoreHitTime (?)
    auto const endTime_{step.getTimePost()}; // time at end point of track.

    if (startTime_ == endTime_) {
      return ProcessReturn::Ok;
    } else {

      // get start and end position of the track
      Point const startPoint_{step.getPositionPre()};
      Point const endPoint_{step.getPositionPost()};

      auto const currDirection{(endPoint_ - startPoint_).normalized()};
      // calculate the track length
      auto const tracklength_{(endPoint_ - startPoint_).getNorm()};

      // beta is velocity / speed of light. Start & end should be the same in endpoints!
      auto const corrBetaValue{(endPoint_ - startPoint_).getNorm() /
                               (constants::c * (endTime_ - startTime_))};
      auto const beta_{currDirection * corrBetaValue};

      // get particle charge
      auto const charge_{get_charge(step.getParticlePre().getPID())};

      // get thinning weight
      auto const thinningWeight{step.getParticlePre().getWeight()};

      // constants for electric field vector calculation
      auto const constants_{charge_ * emConstant_ * thinningWeight};

      // set threshold for application of ZHS-like approximation.
      const double approxThreshold_{1.0e-3};

      // loop over each observer in the observer collection (detector)
      for (auto& observer : observers_.getObservers()) {

        // get the SignalPathCollection (path1) from the start "endpoint" to the observer.
        auto paths1{this->propagator_.propagate(step.getParticlePre(), startPoint_,
                                                observer.getLocation())};

        // get the SignalPathCollection (path2) from the end "endpoint" to the observer.
        auto paths2{this->propagator_.propagate(step.getParticlePre(), endPoint_,
                                                observer.getLocation())};

        // LCOV_EXCL_START
        // This should never happen unless someone implements a bad propagator
        // good to catch this, hard to test without writing a custom and broken propagator
        if (paths1.size() != paths2.size()) {
          CORSIKA_LOG_CRITICAL(
              "Signal Paths do not have the same size in CoREAS. Path starts: {} and "
              "path ends: {}",
              paths1.size(), paths2.size());
        }
        // LCOV_EXCL_STOP

        // loop over both paths at once and directly compare 'start' and 'end'
        // attributes
        for (size_t i = 0; (i < paths1.size() && i < paths2.size()); i++) {

          // calculate preDoppler factor
          double preDoppler_{1.0 - paths1[i].refractive_index_source_ *
                                       beta_.dot(paths1[i].emit_)};

          // check if preDoppler has become zero in case of refractive index of unity
          // because of numerical limitations here you might need std::fabs(preDoppler)
          // in the if statement - same with post & mid
          // LCOV_EXCL_START
          if (preDoppler_ == 0) {
            CORSIKA_LOG_WARN("preDoppler factor numerically zero in COREAS");
            // redo calculation with higher precision
            auto const& beta_vec = beta_.getComponents().getEigenVector();
            auto const& emit_vec = paths1[i].emit_.getComponents().getEigenVector();
            long double const indexL_{paths1[i].refractive_index_source_};
            long double const betaX_{static_cast<double>(beta_vec(0))};
            long double const betaY_{static_cast<double>(beta_vec(1))};
            long double const betaZ_{static_cast<double>(beta_vec(2))};
            long double const startX_{static_cast<double>(emit_vec(0))};
            long double const startY_{static_cast<double>(emit_vec(1))};
            long double const startZ_{static_cast<double>(emit_vec(2))};
            long double const doppler =
                1.0l - indexL_ * (betaX_ * startX_ + betaY_ * startY_ + betaZ_ * startZ_);
            preDoppler_ = doppler;
          }
          // LCOV_EXCL_STOP

          // calculate postDoppler factor
          double postDoppler_{1.0 - paths2[i].refractive_index_source_ *
                                        beta_.dot(paths2[i].emit_)};

          // check if postDoppler has become zero in case of refractive index of unity
          // because of numerical limitations
          // LCOV_EXCL_START
          if (postDoppler_ == 0) {
            CORSIKA_LOG_WARN("postDoppler factor numerically zero in CoREAS");
            // redo calculation with higher precision
            auto const& beta_vec = beta_.getComponents().getEigenVector();
            auto const& emit_vec = paths2[i].emit_.getComponents().getEigenVector();
            long double const indexL_{paths2[i].refractive_index_source_};
            long double const betaX_{static_cast<double>(beta_vec(0))};
            long double const betaY_{static_cast<double>(beta_vec(1))};
            long double const betaZ_{static_cast<double>(beta_vec(2))};
            long double const endX_{static_cast<double>(emit_vec(0))};
            long double const endY_{static_cast<double>(emit_vec(1))};
            long double const endZ_{static_cast<double>(emit_vec(2))};
            long double const doppler =
                1.0l - indexL_ * (betaX_ * endX_ + betaY_ * endY_ + betaZ_ * endZ_);
            postDoppler_ = doppler;
          }
          // LCOV_EXCL_STOP

          // calculate receive time for startpoint (aka time delay)
          auto startPointReceiveTime_{
              startTime_ +
              paths1[i].propagation_time_}; // TODO: time 0 is when the imaginary
                                            // primary hits the ground

          // calculate receive time for endpoint
          auto endPointReceiveTime_{endTime_ + paths2[i].propagation_time_};

          // get unit vector for startpoint at observer location
          auto ReceiveVectorStart_{paths1[i].receive_};

          // get unit vector for endpoint at observer location
          auto ReceiveVectorEnd_{paths2[i].receive_};

          // perform ZHS-like calculation close to Cherenkov angle and for refractive
          // index at observer location greater than 1
          if ((paths1[i].refractive_index_destination_ > 1) &&
              ((std::fabs(preDoppler_) < approxThreshold_) ||
               (std::fabs(postDoppler_) < approxThreshold_))) {
            CORSIKA_LOG_DEBUG("Used ZHS-like approximation in CoREAS - radio");

            auto const halfVector_{(startPoint_ - endPoint_) * 0.5};
            auto const midPoint_{endPoint_ + halfVector_};

            // get global simulation time for the middle point of that track.
            TimeType const midTime_{(startTime_ + endTime_) * 0.5};

            // get the SignalPathCollection (path3) from the middle "endpoint" to the
            // observer.
            auto paths3{this->propagator_.propagate(step.getParticlePre(), midPoint_,
                                                    observer.getLocation())};

            // now loop over the paths for endpoint that we got above
            for (auto const& path : paths3) {

              auto const midPointReceiveTime_{midTime_ + path.propagation_time_};
              double midDoppler_{1.0 -
                                 path.refractive_index_source_ * beta_.dot(path.emit_)};

              // check if midDoppler has become zero because of numerical limitations
              // LCOV_EXCL_START
              if (midDoppler_ == 0) {
                CORSIKA_LOG_WARN("midDoppler factor numerically zero in COREAS");
                // redo calculation with higher precision
                auto const& beta_vec = beta_.getComponents().getEigenVector();
                auto const& emit_vec = path.emit_.getComponents().getEigenVector();
                long double const indexL_{path.refractive_index_source_};
                long double const betaX_{static_cast<double>(beta_vec(0))};
                long double const betaY_{static_cast<double>(beta_vec(1))};
                long double const betaZ_{static_cast<double>(beta_vec(2))};
                long double const midX_{static_cast<double>(emit_vec(0))};
                long double const midY_{static_cast<double>(emit_vec(1))};
                long double const midZ_{static_cast<double>(emit_vec(2))};
                long double const doppler =
                    1.0l - indexL_ * (betaX_ * midX_ + betaY_ * midY_ + betaZ_ * midZ_);
                midDoppler_ = doppler;
              }
              // LCOV_EXCL_STOP

              // CoREAS calculation -> get ElectricFieldVector for "midPoint"
              ElectricFieldVector EVmid_ = (path.emit_.cross(path.emit_.cross(beta_))) /
                                           midDoppler_ / path.R_distance_ * constants_ *
                                           observer.getSampleRate();

              ElectricFieldVector EV1_{EVmid_};
              ElectricFieldVector EV2_{EVmid_ * (-1.0)};

              TimeType deltaT_{tracklength_ / (constants::c * corrBetaValue) *
                               std::fabs(midDoppler_)}; // TODO: Caution with this!

              if (startPointReceiveTime_ <
                  endPointReceiveTime_) // EVstart_ arrives earlier
              {
                startPointReceiveTime_ = midPointReceiveTime_ - 0.5 * deltaT_;
                endPointReceiveTime_ = midPointReceiveTime_ + 0.5 * deltaT_;
              } else // EVend_ arrives earlier
              {
                startPointReceiveTime_ = midPointReceiveTime_ + 0.5 * deltaT_;
                endPointReceiveTime_ = midPointReceiveTime_ - 0.5 * deltaT_;
              }

              TimeType const gridResolution_{1 / observer.getSampleRate()};
              deltaT_ = endPointReceiveTime_ - startPointReceiveTime_;

              // redistribute contributions over time scale defined by the observation
              // time resolution
              if (abs(deltaT_) < (gridResolution_)) {

                EV1_ *= std::fabs((deltaT_ / gridResolution_));
                EV2_ *= std::fabs((deltaT_ / gridResolution_));

                // ToDO: be careful with times in C8!!! where is the zero (time). Is it
                // close-by?
                long const startBin = static_cast<long>(
                    std::floor(startPointReceiveTime_ / gridResolution_ + 0.5l));
                long const endBin = static_cast<long>(
                    std::floor(endPointReceiveTime_ / gridResolution_ + 0.5l));
                double const startBinFraction =
                    (startPointReceiveTime_ / gridResolution_) -
                    std::floor(startPointReceiveTime_ / gridResolution_);
                double const endBinFraction =
                    (endPointReceiveTime_ / gridResolution_) -
                    std::floor(endPointReceiveTime_ / gridResolution_);

                // only do timing modification if contributions would land in same bin
                if (startBin == endBin) {

                  // if startE arrives before endE
                  if ((deltaT_) >= 0_s) {
                    if ((startBinFraction >= 0.5) &&
                        (endBinFraction >= 0.5)) // both points left of bin center
                    {
                      startPointReceiveTime_ -=
                          gridResolution_; // shift EV1_ to previous gridpoint
                    } else if ((startBinFraction < 0.5) &&
                               (endBinFraction < 0.5)) // both points right of bin center
                    {
                      endPointReceiveTime_ +=
                          gridResolution_; // shift EV2_ to next gridpoint
                    } else                 // points on both sides of bin center
                    {
                      double const leftDist = 1.0 - startBinFraction;
                      double const rightDist = endBinFraction;
                      // check if asymmetry to right or left
                      if (rightDist >= leftDist) {
                        endPointReceiveTime_ +=
                            gridResolution_; // shift EV2_ to next gridpoint
                      } else {
                        startPointReceiveTime_ -=
                            gridResolution_; // shift EV1_ to previous gridpoint
                      }
                    }
                  } else // if endE arrives before startE
                  {
                    if ((startBinFraction >= 0.5) &&
                        (endBinFraction >= 0.5)) // both points left of bin center
                    {
                      endPointReceiveTime_ -=
                          gridResolution_; // shift EV2_ to previous gridpoint
                    } else if ((startBinFraction < 0.5) &&
                               (endBinFraction < 0.5)) // both points right of bin center
                    {
                      startPointReceiveTime_ +=
                          gridResolution_; // shift EV1_ to next gridpoint
                    } else                 // points on both sides of bin center
                    {
                      double const leftDist = 1.0 - endBinFraction;
                      double const rightDist = startBinFraction;
                      // check if asymmetry to right or left
                      if (rightDist >= leftDist) {
                        startPointReceiveTime_ +=
                            gridResolution_; // shift EV1_ to next gridpoint
                      } else {
                        endPointReceiveTime_ -=
                            gridResolution_; // shift EV2_ to previous gridpoint
                      }
                    }
                  } // End of else statement
                }   // End of if for startbin == endbin
              }     // End of if deltaT < gridresolution

              // TODO: Be very careful with this. Maybe the EVs should be fed after the
              // for loop of paths3
              observer.receive(startPointReceiveTime_, paths1[i], EV1_);
              observer.receive(endPointReceiveTime_, paths2[i], EV2_);
            } // End of looping over paths3

          } // end of ZHS-like approximation
          else {

            // calculate electric field vector for startpoint
            ElectricFieldVector EV1_ =
                (paths1[i].emit_.cross(paths1[i].emit_.cross(beta_))) / preDoppler_ /
                paths1[i].R_distance_ * constants_ * observer.getSampleRate();

            // calculate electric field vector for endpoint
            ElectricFieldVector EV2_ =
                (paths2[i].emit_.cross(paths2[i].emit_.cross(beta_))) / postDoppler_ /
                paths2[i].R_distance_ * constants_ * (-1.0) * observer.getSampleRate();

            if ((preDoppler_ < 1.e-9) || (postDoppler_ < 1.e-9)) {

              TimeType const gridResolution_{1 / observer.getSampleRate()};
              TimeType deltaT_{endPointReceiveTime_ - startPointReceiveTime_};

              if (abs(deltaT_) < (gridResolution_)) {

                EV1_ *= std::fabs(deltaT_ / gridResolution_); // Todo: rename EV1 and 2
                EV2_ *= std::fabs(deltaT_ / gridResolution_);

                long const startBin = static_cast<long>(
                    std::floor(startPointReceiveTime_ / gridResolution_ + 0.5l));
                long const endBin = static_cast<long>(
                    std::floor(endPointReceiveTime_ / gridResolution_ + 0.5l));
                double const startBinFraction =
                    (startPointReceiveTime_ / gridResolution_) -
                    std::floor(startPointReceiveTime_ / gridResolution_);
                double const endBinFraction =
                    (endPointReceiveTime_ / gridResolution_) -
                    std::floor(endPointReceiveTime_ / gridResolution_);

                // only do timing modification if contributions would land in same bin
                if (startBin == endBin) {

                  if ((startBinFraction >= 0.5) &&
                      (endBinFraction >= 0.5)) // both points left of bin center
                  {
                    startPointReceiveTime_ -=
                        gridResolution_; // shift EV1_ to previous gridpoint
                  } else if ((startBinFraction < 0.5) &&
                             (endBinFraction < 0.5)) // both points right of bin center
                  {
                    endPointReceiveTime_ +=
                        gridResolution_; // shift EV2_ to next gridpoint
                  } else                 // points on both sides of bin center
                  {
                    double const leftDist = 1.0 - startBinFraction;
                    double const rightDist = endBinFraction;
                    // check if asymmetry to right or left
                    if (rightDist >= leftDist) {
                      endPointReceiveTime_ +=
                          gridResolution_; // shift EV2_ to next gridpoint
                    } else {
                      startPointReceiveTime_ -=
                          gridResolution_; // shift EV1_ to previous gridpoint
                    }
                  }

                } // End of if for startbin == endbin
              }   // End of if deltaT < gridresolution
            }     // End of if that checks small doppler factors
            observer.receive(startPointReceiveTime_, paths1[i], EV1_);
            observer.receive(endPointReceiveTime_, paths2[i], EV2_);
          } // End of else that does not perform ZHS-like approximation

        } // End of loop over both paths to get signal info
      }   // End of looping over observer

      return ProcessReturn::Ok;
    }
  } // End of simulate method

} // namespace corsika
