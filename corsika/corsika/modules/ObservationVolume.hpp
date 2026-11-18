/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/writers/ParticleWriterParquet.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/core/Step.hpp>

namespace corsika {

  /**
     @ingroup Modules
     @{

     The ObservationVolume writes PDG codes, kinetic energy, position, and direction
     of particles in the observation frame into its output file. The particles
     are considered "absorbed" afterwards.

     @note Output coordinates are in the volume's coordinate system (Box CS or root CS for
     Sphere).
   */
  template <typename TTracking, typename TVolume,
            typename TOutputWriter = ParticleWriterParquet>
  class ObservationVolume
      : public ContinuousProcess<ObservationVolume<TTracking, TVolume, TOutputWriter>>,
        public TOutputWriter {

  public:
    ObservationVolume(TVolume vol);

    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle>&, bool const stepLimit);

    template <typename TParticle, typename TTrajectory>
    LengthType getMaxStepLength(TParticle const&, TTrajectory const& vTrajectory);

    void showResults() const;
    void reset();
    HEPEnergyType getEnergy() const { return energy_; }
    YAML::Node getConfig() const;

  private:
    TVolume vol_;
    HEPEnergyType energy_;
    unsigned int count_;
  };
  //! @}
} // namespace corsika

#include <corsika/detail/modules/ObservationVolume.inl>