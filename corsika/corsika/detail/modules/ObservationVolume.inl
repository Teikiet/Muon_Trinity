/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

namespace corsika {

  template <typename TTracking, typename TVolume, typename TOutput>
  ObservationVolume<TTracking, TVolume, TOutput>::ObservationVolume(TVolume vol)
      : vol_(vol)
      , energy_(0_GeV)
      , count_(0) {}

  template <typename TTracking, typename TVolume, typename TOutput>
  template <typename TParticle>
  inline ProcessReturn ObservationVolume<TTracking, TVolume, TOutput>::doContinuous(
      Step<TParticle>& step, bool const stepLimit) {
    /*
       The current step did not yet reach the ObservationVolume, do nothing now and
       wait:
     */
    if (!stepLimit) { return ProcessReturn::Ok; }

    HEPEnergyType const kineticEnergy = step.getEkinPost();
    Point const pointOfIntersection = step.getPositionPost();
    DirectionVector const dirction = step.getDirectionPost();
    double const weight = step.getParticlePre().getWeight();

    // add particles to the output file stream
    // Optimization: transform once instead of 6 times
    auto const& cs = vol_.getCoordinateSystem();
    auto const pos = pointOfIntersection.getCoordinates(cs);
    auto const dir = dirction.getComponents(cs);
    this->write(step.getParticlePre().getPID(), kineticEnergy, pos.getX(), pos.getY(),
                pos.getZ(), dir.getX(), dir.getY(), dir.getZ(), step.getTimePost(),
                weight);

    // always absorb particles
    count_++;
    energy_ += kineticEnergy + get_mass(step.getParticlePre().getPID());
    return ProcessReturn::ParticleAbsorbed;
  }

  template <typename TTracking, typename TVolume, typename TOutput>
  template <typename TParticle, typename TTrajectory>
  inline LengthType ObservationVolume<TTracking, TVolume, TOutput>::getMaxStepLength(
      TParticle const& particle, TTrajectory const& trajectory) {

    CORSIKA_LOG_TRACE("getMaxStepLength, particle={}, pos={}, dir={}, Box={}",
                      particle.asString(), particle.getPosition(),
                      particle.getDirection(), vol_.asString());

    auto const intersection = TTracking::intersect(particle, vol_);

    TimeType const timeOfEntry = intersection.getEntry();
    TimeType const timeOfExit = intersection.getExit();
    CORSIKA_LOG_TRACE("timeOfEntry={}, timeOfExit={}", timeOfEntry, timeOfExit);
    if (timeOfEntry < TimeType::zero()) {
      if (timeOfExit < TimeType::zero()) {
        // opposite direction
        return std::numeric_limits<double>::infinity() * 1_m;
      } else {
        // inside box: not zero but 1 pm, to allow short lived particles to decay
        return 1e-12 * 1_m;
      }
    }
    if (timeOfEntry > trajectory.getDuration()) {
      // can not reach
      return std::numeric_limits<double>::infinity() * 1_m;
    }
    double const fractionOfIntersection = timeOfEntry / trajectory.getDuration();
    CORSIKA_LOG_TRACE("ObservationVolume: getMaxStepLength dist={} m, pos={}",
                      trajectory.getLength(fractionOfIntersection) / 1_m,
                      trajectory.getPosition(fractionOfIntersection));
    return trajectory.getLength(fractionOfIntersection);
  }

  template <typename TTracking, typename TVolume, typename TOutput>
  inline void ObservationVolume<TTracking, TVolume, TOutput>::showResults() const {
    CORSIKA_LOG_INFO(
        " ******************************\n"
        " ObservationVolume: \n"
        " energy at Box (GeV)     :  {}\n"
        " no. of particles at Box :  {}\n"
        " ******************************",
        energy_ / 1_GeV, count_);
  }

  template <typename TTracking, typename TVolume, typename TOutput>
  inline YAML::Node ObservationVolume<TTracking, TVolume, TOutput>::getConfig() const {
    using namespace units::si;
    auto const& cs = vol_.getCoordinateSystem();
    auto center = vol_.getCenter();

    // construct the top-level node
    YAML::Node node;

    // basic info
    node["type"] = "ObservationVolume";
    node["units"]["length"] = "m";
    node["units"]["energy"] = "GeV";
    node["units"]["time"] = "s";

    // center coordinates in volume CS
    auto const center_coords = center.getCoordinates(cs);
    node["center"].push_back(center_coords.getX() / 1_m);
    node["center"].push_back(center_coords.getY() / 1_m);
    node["center"].push_back(center_coords.getZ() / 1_m);

    // volume axes in volume CS
    DirectionVector const x_axis = DirectionVector{cs, {1, 0, 0}};
    auto const x_coords = x_axis.getComponents(cs);
    node["x-axis"].push_back(x_coords.getX().magnitude());
    node["x-axis"].push_back(x_coords.getY().magnitude());
    node["x-axis"].push_back(x_coords.getZ().magnitude());

    DirectionVector const y_axis = DirectionVector{cs, {0, 1, 0}};
    auto const y_coords = y_axis.getComponents(cs);
    node["y-axis"].push_back(y_coords.getX().magnitude());
    node["y-axis"].push_back(y_coords.getY().magnitude());
    node["y-axis"].push_back(y_coords.getZ().magnitude());

    DirectionVector const z_axis = DirectionVector{cs, {0, 0, 1}};
    auto const z_coords = z_axis.getComponents(cs);
    node["z-axis"].push_back(z_coords.getX().magnitude());
    node["z-axis"].push_back(z_coords.getY().magnitude());
    node["z-axis"].push_back(z_coords.getZ().magnitude());

    return node;
  }

  template <typename TTracking, typename TVolume, typename TOutput>
  inline void ObservationVolume<TTracking, TVolume, TOutput>::reset() {
    energy_ = 0_GeV;
    count_ = 0;
  }

} // namespace corsika