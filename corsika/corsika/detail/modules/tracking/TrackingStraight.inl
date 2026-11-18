/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Intersections.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/StraightTrajectory.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/modules/tracking/Intersect.hpp>

#include <cmath>
#include <type_traits>
#include <utility>

namespace corsika::tracking_line {

  template <typename TParticle>
  inline auto Tracking::getTrack(TParticle const& particle) {
    VelocityVector const initialVelocity =
        particle.getMomentum() / particle.getEnergy() * constants::c;

    auto const& initialPosition = particle.getPosition();
    CORSIKA_LOG_DEBUG(
        "TrackingStraight pid: {}"
        " , E = {} GeV \n"
        "\tTracking pos: {} \n"
        "\tTracking   p: {} GeV \n"
        "\tTracking   v: {}",
        particle.getPID(), particle.getEnergy() / 1_GeV, initialPosition.getCoordinates(),
        particle.getMomentum().getComponents() / 1_GeV, initialVelocity.getComponents());

    // traverse the environment volume tree and find next
    // intersection
    auto [minTime, minNode] = nextIntersect(particle);

    return std::make_tuple(StraightTrajectory(Line(initialPosition, initialVelocity),
                                              minTime), // trajectory
                           minNode);                    // next volume node
  }

  template <typename TParticle>
  inline Intersections Tracking::intersect(TParticle const& particle,
                                           Sphere const& sphere) {
    auto const& position = particle.getPosition();
    auto const delta = position - sphere.getCenter();
    auto const velocity = particle.getMomentum() / particle.getEnergy() * constants::c;
    auto const vSqNorm = velocity.getSquaredNorm();
    auto const R = sphere.getRadius();

    auto const vDotDelta = velocity.dot(delta);
    auto const discriminant =
        vDotDelta * vDotDelta - vSqNorm * (delta.getSquaredNorm() - R * R);

    if (discriminant.magnitude() > 0) {
      auto const sqDisc = sqrt(discriminant);
      auto const invDenom = 1 / vSqNorm;

      CORSIKA_LOG_TRACE("numericallyInside={}", sphere.contains(position));
      return Intersections((-vDotDelta - sqDisc) * invDenom,
                           (-vDotDelta + sqDisc) * invDenom);
    }
    return Intersections();
  }

  template <typename TParticle>
  inline Intersections Tracking::intersect(TParticle const& particle, Box const& box) {
    Point const& position = particle.getPosition();
    VelocityVector const velocity =
        particle.getMomentum() / particle.getEnergy() * constants::c;

    // Transform once, access components multiple times
    // This reduces 6 transformations (3 for position + 3 for velocity) to just 2
    auto const& cs = box.getCoordinateSystem();
    auto const posInBox = position.getCoordinates(cs);
    auto const velInBox = velocity.getComponents(cs);

    LengthType x0 = posInBox.getX();
    LengthType y0 = posInBox.getY();
    LengthType z0 = posInBox.getZ();
    SpeedType vx = velInBox.getX();
    SpeedType vy = velInBox.getY();
    SpeedType vz = velInBox.getZ();

    CORSIKA_LOG_TRACE(
        "particle in box coordinate: position: ({:.3f}, {:.3f}, "
        "{:.3f}) m, veolocity: ({:.3f}, {:.3f}, {:.3f}) m/ns",
        x0 / 1_m, y0 / 1_m, z0 / 1_m, vx / (1_m / 1_ns), vy / (1_m / 1_ns),
        vz / (1_m / 1_ns));

    auto get_intersect_min_max = [](LengthType x0, SpeedType v0, LengthType dx) {
      auto t1 = (dx - x0) / v0;
      auto t2 = (-dx - x0) / v0;
      if (t1 > t2)
        return std::make_pair(t1, t2);
      else
        return std::make_pair(t2, t1);
    };

    auto [tx_max, tx_min] = get_intersect_min_max(x0, vx, box.getX());
    auto [ty_max, ty_min] = get_intersect_min_max(y0, vy, box.getY());
    auto [tz_max, tz_min] = get_intersect_min_max(z0, vz, box.getZ());

    TimeType t_exit = std::min(std::min(tx_max, ty_max), tz_max);
    TimeType t_enter = std::max(std::max(tx_min, ty_min), tz_min);

    CORSIKA_LOG_DEBUG("t_enter: {} ns, t_exit: {} ns", t_enter / 1_ns, t_exit / 1_ns);
    if ((t_exit > t_enter)) {
      if (t_enter < 0_s && t_exit > 0_s)
        CORSIKA_LOG_DEBUG("numericallyInside={}", box.contains(position));
      else if (t_enter < 0_s && t_exit < 0_s)
        CORSIKA_LOG_DEBUG("oppisite direction");
      return Intersections(std::move(t_enter), std::move(t_exit));
    } else
      return Intersections();
  }

  template <typename TParticle, typename TBaseNodeType>
  inline Intersections Tracking::intersect(TParticle const& particle,
                                           TBaseNodeType const& volumeNode) {
    if (Sphere const* sphere = dynamic_cast<Sphere const*>(&volumeNode.getVolume());
        sphere) {
      return Tracking::intersect<TParticle>(particle, *sphere);
    } else if (Box const* box = dynamic_cast<Box const*>(&volumeNode.getVolume()); box) {
      return Tracking::intersect<TParticle>(particle, *box);
    } else if (SeparationPlane const* sepPlane =
                   dynamic_cast<SeparationPlane const*>(&volumeNode.getVolume());
               sepPlane) {
      return Tracking::intersect<TParticle>(particle, *sepPlane);
    } else {
      throw std::runtime_error(
          "The Volume type provided is not supported in "
          "TrackingStraight::intersect(particle, node)");
    }
  }

  template <typename TParticle>
  inline Intersections Tracking::intersect(TParticle const& particle,
                                           Plane const& plane) {
    auto const delta = plane.getCenter() - particle.getPosition();
    auto const velocity = particle.getMomentum() / particle.getEnergy() * constants::c;
    auto const n = plane.getNormal();
    auto const n_dot_v = n.dot(velocity);

    CORSIKA_LOG_TRACE("n_dot_v={}, delta={}, momentum={}", n_dot_v, delta,
                      particle.getMomentum());

    if (n_dot_v.magnitude() == 0)
      return Intersections();
    else
      return Intersections(n.dot(delta) / n_dot_v);
  }

  template <typename TParticle>
  inline Intersections Tracking::intersect(TParticle const& particle,
                                           SeparationPlane const& sepPlane) {
    return intersect(particle, sepPlane.getPlane());
  }

} // namespace corsika::tracking_line