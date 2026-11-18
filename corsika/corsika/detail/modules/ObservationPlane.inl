/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

namespace corsika {

  template <typename TTracking, typename TOutput>
  template <typename... TArgs>
  ObservationPlane<TTracking, TOutput>::ObservationPlane(Plane const& obsPlane,
                                                         DirectionVector const& x_axis,
                                                         bool const deleteOnHit,
                                                         LengthType const padding,
                                                         TArgs&&... args)
      : TOutput(std::forward<TArgs>(args)...)
      , plane_(obsPlane)
      , xAxis_(x_axis.normalized())
      , yAxis_(obsPlane.getNormal().cross(xAxis_))
      , deleteOnHit_(deleteOnHit)
      , padding_(padding) {}

  template <typename TTracking, typename TOutput>
  inline void ObservationPlane<TTracking, TOutput>::setGlobalOutput(
      Point const& origin) {
    globalOutputOrigin_ = origin;
  }

  template <typename TTracking, typename TOutput>
  template <typename TParticle>
  inline ProcessReturn ObservationPlane<TTracking, TOutput>::doContinuous(
      Step<TParticle>& step, bool const stepLimit) {
    /*
       The current step did not yet reach the ObservationPlane, do nothing now and wait:
     */
    if (!stepLimit) {
      // @todo this is actually needed to fix small instabilities of the leap-frog
      // tracking: Note, this is NOT a general solution and should be clearly revised with
      // a more robust tracking. #ifdef DEBUG
      if (deleteOnHit_) {
        // since this is basically a bug, it cannot be tested LCOV_EXCL_START
        LengthType const check =
            (step.getPositionPost() - plane_.getCenter()).dot(plane_.getNormal());
        if (check < 0_m) {
          CORSIKA_LOG_WARN("PARTICLE AVOIDED OBSERVATIONPLANE {}", check);
          CORSIKA_LOG_WARN("Temporary fix: write and remove particle.");
        } else
          return ProcessReturn::Ok;
        // LCOV_EXCL_STOP
      } else
        // #endif
        return ProcessReturn::Ok;
    }

    HEPEnergyType const kineticEnergy = step.getEkinPost();
    Point const pointOfIntersection = step.getPositionPost();
    DirectionVector const direction = step.getDirectionPost();

    // add our particles to the output file stream
    double const weight = step.getParticlePre().getWeight();

    if (globalOutputOrigin_) {
      // Global coordinates: position and direction are expressed in the coordinate
      // system of the observation point (the root system), translated so that the
      // configured origin sits at (0, 0, 0).
      Vector const globalDisplacement = pointOfIntersection - *globalOutputOrigin_;
      QuantityVector<length_d> const& gc = globalDisplacement.getComponents();
      auto const dc = direction.getComponents();
      this->write(step.getParticlePre().getPID(), kineticEnergy, gc[0], gc[1], gc[2],
                  static_cast<double>(dc[0]), static_cast<double>(dc[1]),
                  static_cast<double>(dc[2]), step.getTimePost(), weight);
    } else {
      Vector const displacement = pointOfIntersection - plane_.getCenter();
      this->write(step.getParticlePre().getPID(), kineticEnergy,
                  displacement.dot(xAxis_), displacement.dot(yAxis_), 0_m,
                  direction.dot(xAxis_), direction.dot(yAxis_),
                  direction.dot(plane_.getNormal()), step.getTimePost(), weight);
    }

    CORSIKA_LOG_TRACE("Particle detected absorbed={}", deleteOnHit_);

    if (deleteOnHit_) {
      return ProcessReturn::ParticleAbsorbed;
    } else {
      // due to numerical precision of tracker, can fall into a Zeno's paradox
      // ensure that the particle is on the far side to avoid continuously hitting this
      // plane
      step.add_displacement(step.getDisplacement().normalized() * padding_);
      return ProcessReturn::Ok;
    }
  } // namespace corsika

  template <typename TTracking, typename TOutput>
  template <typename TParticle, typename TTrajectory>
  inline LengthType ObservationPlane<TTracking, TOutput>::getMaxStepLength(
      TParticle const& particle, TTrajectory const& trajectory) {

    CORSIKA_LOG_TRACE("getMaxStepLength, particle={}, pos={}, dir={}, plane={}",
                      particle.asString(), particle.getPosition(),
                      particle.getDirection(), plane_.asString());

    auto const intersection = TTracking::intersect(particle, plane_);

    TimeType const timeOfIntersection = intersection.getEntry();
    CORSIKA_LOG_TRACE("timeOfIntersection={}", timeOfIntersection);
    if (timeOfIntersection <= TimeType::zero()) {
      return std::numeric_limits<double>::infinity() * 1_m;
    }
    if (timeOfIntersection > trajectory.getDuration()) {
      return std::numeric_limits<double>::infinity() * 1_m;
    }
    double const fractionOfIntersection = timeOfIntersection / trajectory.getDuration();
    CORSIKA_LOG_TRACE("ObservationPlane: getMaxStepLength dist={} m, pos={}",
                      trajectory.getLength(fractionOfIntersection) / 1_m,
                      trajectory.getPosition(fractionOfIntersection));
    return trajectory.getLength(fractionOfIntersection);
  }

  template <typename TTracking, typename TOutput>
  inline YAML::Node ObservationPlane<TTracking, TOutput>::getConfig() const {
    using namespace units::si;

    // construct the top-level node
    YAML::Node node;

    // basic info
    node["type"] = "ObservationPlane";
    node["units"]["length"] = "m"; // add default units for values
    node["units"]["energy"] = "GeV";
    node["units"]["time"] = "s";

    // the center of the plane
    auto const center{plane_.getCenter()};

    // coordinates are stored in root CS
    auto const center_coords{center.getCoordinates()};
    node["plane"]["center"].push_back(center_coords.getX() / 1_m);
    node["plane"]["center"].push_back(center_coords.getY() / 1_m);
    node["plane"]["center"].push_back(center_coords.getZ() / 1_m);

    // the normal vector of the plane
    auto const normal{plane_.getNormal().getComponents()};
    node["plane"]["normal"].push_back(normal.getX().magnitude());
    node["plane"]["normal"].push_back(normal.getY().magnitude());
    node["plane"]["normal"].push_back(normal.getZ().magnitude());

    // the x-axis vector (stored in root CS)
    auto const xAxis_coords{xAxis_.getComponents()};
    node["x-axis"].push_back(xAxis_coords.getX().magnitude());
    node["x-axis"].push_back(xAxis_coords.getY().magnitude());
    node["x-axis"].push_back(xAxis_coords.getZ().magnitude());

    // the y-axis vector (stored in root CS)
    auto const yAxis_coords{yAxis_.getComponents()};
    node["y-axis"].push_back(yAxis_coords.getX().magnitude());
    node["y-axis"].push_back(yAxis_coords.getY().magnitude());
    node["y-axis"].push_back(yAxis_coords.getZ().magnitude());

    node["delete_on_hit"] = deleteOnHit_;
    node["padding"] = padding_ / 1_m;

    // coordinate frame of the recorded values
    if (globalOutputOrigin_) {
      auto const origin_coords{globalOutputOrigin_->getCoordinates()};
      node["global_output"]["enabled"] = true;
      node["global_output"]["origin"].push_back(origin_coords.getX() / 1_m);
      node["global_output"]["origin"].push_back(origin_coords.getY() / 1_m);
      node["global_output"]["origin"].push_back(origin_coords.getZ() / 1_m);
    } else {
      node["global_output"]["enabled"] = false;
    }

    return node;
  }

} // namespace corsika