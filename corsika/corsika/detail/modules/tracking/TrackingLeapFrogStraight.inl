/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/StraightTrajectory.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/modules/tracking/TrackingStraight.hpp>

#include <type_traits>
#include <utility>

namespace corsika {

  namespace tracking_leapfrog_straight {

    template <typename Particle>
    inline auto Tracking::getTrack(Particle& particle) {
      VelocityVector initialVelocity =
          particle.getMomentum() / particle.getEnergy() * constants::c;

      Point const initialPosition = particle.getPosition();

      CORSIKA_LOG_DEBUG(
          "TrackingLeapFrogStraight pid: {}"
          " , E = {} GeV \n"
          "\tTracking pos: {} \n"
          "\tTracking   p: {} GeV \n"
          "\tTracking   v: {} ",
          particle.getPID(), particle.getEnergy() / 1_GeV,
          initialPosition.getCoordinates(),
          particle.getMomentum().getComponents() / 1_GeV,
          initialVelocity.getComponents());

      typedef decltype(particle.getNode()) node_type;
      node_type const volumeNode = particle.getNode();

      // check if particle is moving at all
      auto const absVelocity = initialVelocity.getNorm();
      if (absVelocity * 1_s == 0_m) {
        return std::make_tuple(
            StraightTrajectory(Line(initialPosition, initialVelocity), 0_s), volumeNode);
      }

      // charge of the particle, and magnetic field
      auto const charge = particle.getCharge();
      auto magneticfield =
          volumeNode->getModelProperties().getMagneticField(initialPosition);
      auto const magnitudeB = magneticfield.getNorm();
      CORSIKA_LOG_DEBUG("field={} uT, charge={}, magnitudeB={} uT",
                        magneticfield.getComponents() / 1_uT, charge, magnitudeB / 1_T);
      bool const no_deflection = (charge == 0 * constants::e) || magnitudeB == 0_T;

      // check, where the first halve-step direction has geometric intersections
      auto const [initialTrack, initialTrackNextVolume] =
          tracking_line::Tracking::getTrack(particle);
      //{ [[maybe_unused]] auto& initialTrackNextVolume_dum = initialTrackNextVolume; }
      auto const initialTrackLength = initialTrack.getLength(1);

      CORSIKA_LOG_DEBUG("initialTrack(0)={}, initialTrack(1)={}, initialTrackLength={}",
                        initialTrack.getPosition(0).getCoordinates(),
                        initialTrack.getPosition(1).getCoordinates(), initialTrackLength);

      // if particle is non-deflectable, we are done:
      if (no_deflection) {
        CORSIKA_LOG_DEBUG("no deflection. tracking finished");
        return std::make_tuple(initialTrack, initialTrackNextVolume);
      }

      HEPMomentumType const p_perp =
          (particle.getMomentum() -
           particle.getMomentum().getParallelProjectionOnto(magneticfield))
              .getNorm();

      if (p_perp == 0_GeV) {
        // particle travels along, parallel to magnetic field. Rg is
        // "0", return straight track here.
        CORSIKA_LOG_TRACE("p_perp is 0_GeV --> parallel");
        return std::make_tuple(initialTrack, initialTrackNextVolume);
      }

      LengthType const gyroradius = (convert_HEP_to_SI<MassType::dimension_type>(p_perp) *
                                     constants::c / (abs(charge) * magnitudeB));

      // we need to limit maximum step-length since we absolutely
      // need to follow strongly curved trajectories segment-wise,
      // at least if we don't employ concepts as "Helix
      // Trajectories" or similar
      double const maxRadians = 0.001;
      LengthType const steplimit = 2 * cos(maxRadians) * sin(maxRadians) * gyroradius;
      CORSIKA_LOG_DEBUG("gyroradius {}, Steplimit: {}", gyroradius, steplimit);

      // calculate first halve step for "steplimit"
      DirectionVector const initialDirection = particle.getDirection();

      // avoid any intersections within first halve steplength
      // aim for intersection in second halve step
      LengthType const firstHalveSteplength =
          std::min(steplimit / 2, initialTrackLength * firstFraction_);

      CORSIKA_LOG_DEBUG("first halve step length {}, steplimit={}, initialTrackLength={}",
                        firstHalveSteplength, steplimit, initialTrackLength);
      // perform the first halve-step
      Point const position_mid =
          initialPosition + initialDirection * firstHalveSteplength;
      auto const k = charge / (constants::c * convert_HEP_to_SI<MassType::dimension_type>(
                                                  particle.getMomentum().getNorm()));
      DirectionVector const new_direction =
          initialDirection +
          initialDirection.cross(magneticfield) * firstHalveSteplength * 2 * k;
      auto const new_direction_norm = new_direction.getNorm(); // by design this is >1
      CORSIKA_LOG_DEBUG(
          "position_mid={}, new_direction={}, (new_direction_norm)={}, deflection={}",
          position_mid.getCoordinates(), new_direction.getComponents(),
          new_direction_norm,
          acos(std::min(1.0, initialDirection.dot(new_direction) / new_direction_norm)) *
              180 / M_PI);

      // check, where the second halve-step direction has geometric intersections
      particle.setPosition(position_mid);
      particle.setDirection(new_direction);
      auto const [finalTrack, finalTrackNextVolume] =
          tracking_line::Tracking::getTrack(particle);
      particle.setPosition(initialPosition);   // this is not nice...
      particle.setDirection(initialDirection); // this is not nice...

      LengthType const finalTrackLength = finalTrack.getLength(1);
      LengthType const secondLeapFrogLength = firstHalveSteplength * new_direction_norm;

      // check if volume transition is obvious, OR
      // for numerical reasons, particles slighly bend "away" from a
      // volume boundary have a very hard time to cross the border,
      // thus, if secondLeapFrogLength is just slighly shorter (1e-4m) than
      // finalTrackLength we better just [extend the
      // secondLeapFrogLength slightly and] force the volume
      // crossing:
      bool const switch_volume = finalTrackLength - 0.0001_m <= secondLeapFrogLength;
      LengthType const secondHalveStepLength =
          std::min(secondLeapFrogLength, finalTrackLength);

      CORSIKA_LOG_DEBUG(
          "finalTrack(0)={}, finalTrack(1)={}, finalTrackLength={}, "
          "secondLeapFrogLength={}, secondHalveStepLength={}, "
          "secondLeapFrogLength-finalTrackLength={}, "
          "secondHalveStepLength-finalTrackLength={}, "
          "nextVol={}, transition={}",
          finalTrack.getPosition(0).getCoordinates(),
          finalTrack.getPosition(1).getCoordinates(), finalTrackLength,
          secondLeapFrogLength, secondHalveStepLength,
          secondLeapFrogLength - finalTrackLength,
          secondHalveStepLength - finalTrackLength, fmt::ptr(finalTrackNextVolume),
          switch_volume);

      // perform the second halve-step
      auto const new_direction_normalized = new_direction.normalized();
      Point const finalPosition =
          position_mid + new_direction_normalized * secondHalveStepLength;

      LengthType const totalStep = firstHalveSteplength + secondHalveStepLength;
      auto const delta_pos = finalPosition - initialPosition;
      auto const distance = delta_pos.getNorm();

      return std::make_tuple(
          StraightTrajectory(
              Line(initialPosition,
                   (distance == 0_m ? initialVelocity
                                    : delta_pos.normalized() * absVelocity)),
              distance / absVelocity,  // straight distance
              totalStep / absVelocity, // bend distance
              initialVelocity,
              new_direction_normalized * absVelocity), // trajectory
          (switch_volume ? finalTrackNextVolume : volumeNode));
    }

  } // namespace tracking_leapfrog_straight

} // namespace corsika
