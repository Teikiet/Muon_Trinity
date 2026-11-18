/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/logging/Logging.h>
#include <corsika/stack/history/HistoryObservationPlane.hpp>

#include <boost/histogram/ostream.hpp>

namespace corsika::history {

  template <typename TStack>
  inline HistoryObservationPlane<TStack>::HistoryObservationPlane(TStack const& stack,
                                                                  Plane const& obsPlane,
                                                                  bool deleteOnHit)
      : stack_{stack}
      , plane_{obsPlane}
      , deleteOnHit_{deleteOnHit} {}

  template <typename TStack>
  template <typename TParticle, typename TTrajectory>
  inline ProcessReturn HistoryObservationPlane<TStack>::doContinuous(
      TParticle const& particle, TTrajectory const& trajectory) {
    TimeType const timeOfIntersection =
        (plane_.getCenter() - trajectory.getR0()).dot(plane_.getNormal()) /
        trajectory.getV0().dot(plane_.getNormal());

    if (timeOfIntersection < TimeType::zero()) { return ProcessReturn::Ok; }

    if (plane_.isAbove(trajectory.getR0()) == plane_.isAbove(trajectory.getPosition(1))) {
      return ProcessReturn::Ok;
    }

    CORSIKA_LOG_DEBUG(fmt::format("HistoryObservationPlane: Particle detected: pid={}",
                                  particle.getPID()));

    auto const pid = particle.getPID();
    if (particles::isMuon(pid)) { fillHistoryHistogram(particle); }

    if (deleteOnHit_) {
      return ProcessReturn::ParticleAbsorbed;
    } else {
      return ProcessReturn::Ok;
    }
  }

  template <typename TStack>
  template <typename TParticle, typename TTrajectory>
  inline LengthType HistoryObservationPlane<TStack>::getMaxStepLength(
      TParticle const&, TTrajectory const& trajectory) {
    TimeType const timeOfIntersection =
        (plane_.getCenter() - trajectory.getR0()).dot(plane_.getNormal()) /
        trajectory.getV0().dot(plane_.getNormal());

    if (timeOfIntersection < TimeType::zero()) {
      return std::numeric_limits<double>::infinity() * 1_m;
    }

    auto const pointOfIntersection = trajectory.getPosition(timeOfIntersection);
    return (trajectory.getR0() - pointOfIntersection).norm() * 1.0001;
  }

  template <typename TStack>
  template <typename TParticle>
  inline void HistoryObservationPlane<TStack>::fillHistoryHistogram(
      TParticle const& muon) {
    double const muon_energy = muon.getEnergy() / 1_GeV;

    int genctr{0};
    Event const* event = muon.getEvent().get();
    while (event) {
      auto const projectile = stack_.cfirst() + event->projectileIndex();
      if (event->eventType() == EventType::Interaction) {
        genctr++;
        double const projEnergy = projectile.getEnergy() / 1_GeV;
        int const pdg = static_cast<int>(particles::getPDG(projectile.getPID()));

        histogram_(muon_energy, projEnergy, pdg);
      }
      event = event->parentEvent().get(); // projectile.getEvent().get();
    }
  }
} // namespace corsika::history
