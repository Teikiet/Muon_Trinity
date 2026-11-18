
/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <PROPOSAL/PROPOSAL.h>

#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/modules/proposal/ContinuousProcess.hpp>
#include <corsika/modules/proposal/InteractionModel.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/COMBoost.hpp>
#include <corsika/framework/core/Logging.hpp>

namespace corsika::proposal {

  template <typename TOutput>
  inline void ContinuousProcess<TOutput>::buildCalculator(Code code,
                                                          size_t const& component_hash) {
    // search crosssection builder for given particle
    auto p_cross = cross.find(code);
    if (p_cross == cross.end())
      throw std::runtime_error("PROPOSAL could not find corresponding builder");
    if (code == Code::Photon) return; // no continuous builders needed for photons

    // interpolate the crosssection for given media and energy cut. These may
    // take some minutes if you have to build the tables and cannot read the tables
    // from disk
    auto c = p_cross->second(media.at(component_hash), proposal_energycutsettings[code]);

    // choose multiple scattering model
    static constexpr auto ms_type = PROPOSAL::MultipleScatteringType::MoliereInterpol;

    // Build displacement integral and scattering object and interpolate them too and
    // saved in the calc map by a key build out of a hash of composed of the component and
    // particle code.
    auto calculator = Calculator{PROPOSAL::make_displacement(c, true),
                                 PROPOSAL::make_multiple_scattering(
                                     ms_type, particle[code], media.at(component_hash))};
    calc[std::make_pair(component_hash, code)] = std::move(calculator);
  }

  template <typename TOutput>
  template <typename TEnvironment, typename... TOutputArgs>
  inline ContinuousProcess<TOutput>::ContinuousProcess(TEnvironment const& _env,
                                                       TOutputArgs&&... args)
      : ProposalProcessBase(_env)
      , TOutput(args...) {
    //! Initialize PROPOSAL tables for all media and all particles
    for (auto const& medium : media) {
      for (auto const particle_code : tracked) {
        buildCalculator(particle_code, medium.first);
      }
    }
  }

  template <typename TOutput>
  template <typename TParticle>
  inline void ContinuousProcess<TOutput>::scatter(Step<TParticle>& step,
                                                  HEPEnergyType const& loss,
                                                  GrammageType const& grammage) {

    // get or build corresponding calculators
    auto c = getCalculator(step.getParticlePre(), calc);

    auto initial_particle_dir = step.getDirectionPre();

    auto E_i_total = step.getEkinPre() + step.getParticlePre().getMass();
    auto E_f_total = E_i_total - loss;

    // sample scattering_angle, which is the combination sqrt(theta_x^2 + theta_y^2),
    // where theta_x and theta_y itself are independent angles drawn from the multiple
    // scattering distribution
    std::uniform_real_distribution<double> distr(0., 1.);
    auto scattering_angle = (c->second).scatter->CalculateScatteringAngle2D(
        grammage / 1_g * square(1_cm), E_i_total / 1_MeV, E_f_total / 1_MeV, distr(RNG_),
        distr(RNG_));

    // construct vector that is normal to initial direction (in root CS).
    DirectionVector normal_vec{{0, 0, 0}};
    if (std::abs(initial_particle_dir.getX()) > 0.1) {
      normal_vec = {{-initial_particle_dir.getY(), initial_particle_dir.getX(), 0}};
    } else {
      // if x is small, use y and z to construct normal vector
      normal_vec = {{0, -initial_particle_dir.getZ(), initial_particle_dir.getY()}};
    }

    auto const axis1 = normal_vec.getComponents();
    if (axis1.getEigenVector().isZero()) {
      throw std::runtime_error("null-vector given as axis parameter");
    }

    // rotation of zenith by moliere_angle
    Eigen::Matrix3d const rotation1{
        Eigen::AngleAxisd(-scattering_angle, axis1.getEigenVector().normalized())
            .toRotationMatrix()};
    std::uniform_real_distribution<double> distr_azimuth(0., 2 * M_PI);

    // rotation of azimuth by random angle between 0 and 2*PI
    double const random_angle = distr_azimuth(RNG_);
    auto const axis2 = initial_particle_dir.getComponents();

    Eigen::Matrix3d const rotation2{
        Eigen::AngleAxisd(-random_angle, axis2.getEigenVector().normalized())
            .toRotationMatrix()};

    // Rotation to get the scattered_particle_dir
    QuantityVector<dimensionless_d> scattered_vector = QuantityVector<dimensionless_d>(
        rotation2 * rotation1 * initial_particle_dir.getComponents().getEigenVector());

    // update particle direction after continuous loss caused by multiple
    // scattering
    DirectionVector diff_dir_{scattered_vector - axis2.getEigenVector()};
    step.add_dU(diff_dir_);
  }

  template <typename TOutput>
  template <typename TParticle>
  inline ProcessReturn ContinuousProcess<TOutput>::doContinuous(Step<TParticle>& step,
                                                                bool const) {
    if (!canInteract(step.getParticlePre().getPID())) return ProcessReturn::Ok;
    if (step.getDisplacement().getSquaredNorm() == static_pow<2>(0_m))
      return ProcessReturn::Ok;
    if (step.getParticlePre().getPID() == Code::Photon)
      return ProcessReturn::Ok; // no continuous energy losses, no scattering for photons
    // calculate passed grammage
    auto dX = step.getParticlePre().getNode()->getModelProperties().getIntegratedGrammage(
        step.getStraightTrack());

    // get or build corresponding track integral calculator and solve the
    // integral
    auto c = getCalculator(step.getParticlePre(), calc);
    auto E_i_total = (step.getEkinPre() + step.getParticlePre().getMass());
    auto E_f_total = (c->second).disp->UpperLimitTrackIntegral(
                         E_i_total * (1 / 1_MeV), dX * ((1 / 1_g) * 1_cm * 1_cm)) *
                     1_MeV;
    auto dE = E_i_total - E_f_total;

    // if the particle has a charge take multiple scattering into account
    if (step.getParticlePre().getChargeNumber() != 0) scatter(step, dE, dX);
    step.add_dEkin(-dE); // on the stack, this is just kinetic energy, E-m

    // also send to output
    TOutput::write(step.getPositionPre(), step.getPositionPost(),
                   step.getParticlePre().getPID(),
                   step.getParticlePre().getWeight() * dE);

    return ProcessReturn::Ok;
  }

  template <typename TOutput>
  template <typename TParticle, typename TTrajectory>
  inline LengthType ContinuousProcess<TOutput>::getMaxStepLength(
      TParticle const& vP, TTrajectory const& track) {
    auto const code = vP.getPID();
    if (!canInteract(code)) return meter * std::numeric_limits<double>::infinity();
    if (code == Code::Photon)
      return meter *
             std::numeric_limits<double>::infinity(); // no step limitation for photons

    // Limit the step size of a conitnuous loss. The maximal continuous loss seems to be
    // a hyper parameter which must be adjusted.
    //
    auto const energy = vP.getEnergy();
    auto const energy_lim = std::max(
        energy * 0.9, // either 10% relative loss max., or
        (get_kinetic_energy_propagation_threshold(code) +
         get_mass(code)) // energy thresholds globally defined for individual particles
            * 0.9999     // need to go slightly below global e-cut to assure removal
                         // in ParticleCut. This does not matter since at cut-time
                         // the entire energy is removed.
    );

    // solving the track integral for giving energy lim
    auto c = getCalculator(vP, calc);
    auto grammage =
        (c->second).disp->SolveTrackIntegral(energy / 1_MeV, energy_lim / 1_MeV) * 1_g /
        square(1_cm);

    // return it in distance aequivalent
    auto dist =
        vP.getNode()->getModelProperties().getArclengthFromGrammage(track, grammage);
    CORSIKA_LOG_TRACE("PROPOSAL::getMaxStepLength X={} g/cm2, l={} m ",
                      grammage / 1_g * square(1_cm), dist / 1_m);
    if (dist < 0_m) {
      CORSIKA_LOG_WARN(
          "PROPOSAL::getMaxStepLength calculated a negative step length of l={} m. "
          "Return 0_m instead.",
          dist / 1_m);
      return 0_m;
    }
    return dist;
  }

  template <typename TOutput>
  inline YAML::Node ContinuousProcess<TOutput>::getConfig() const {
    return YAML::Node();
  }

} // namespace corsika::proposal