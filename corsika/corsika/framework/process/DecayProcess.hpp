/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/detail/framework/process/DecayProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
     @ingroup Processes
     @{

     Process decribing the decay of particles

     Create a new DecayProcess, e.g. for XYModel, via
     @code
     class XYModel : public DecayProcess<XYModel> {};
     @endcode

     and provide the two necessary interface methods
     @code
     template <typename TSecondaryView>
     void XYModel::doDecay(TSecondaryView&);

     template <typename TParticle>
     TimeType getLifetime(TParticle const&)
     @endcode

     Where, of course, SecondaryView and Particle are the valid
     classes to access particles on the Stack. In user code those two methods do
     not need to be templated, they could use the types
     e.g. corsika::setup::Stack::particle_type -- but by the cost of
     loosing all flexibility otherwise provided.

     SecondaryView allows to retrieve the properties of the projectile
     particles, AND to store new particles (secondaries) which then
     subsequently can be processes by SecondariesProcess. This is how
     the output of decays can be studied right away.
   */

  template <typename TDerived>
  struct DecayProcess : BaseProcess<TDerived> {
  public:
    using BaseProcess<TDerived>::getRef;

    template <typename TParticle>
    InverseTimeType getInverseLifetime(TParticle const& particle) {

      // interface checking on TProcess1
      static_assert(has_method_getLifetime_v<TDerived, TimeType, TParticle const&>,
                    "TDerived has no method with correct signature \"GrammageType "
                    "getInteractionLength(TParticle const&)\" required for "
                    "InteractionProcess<TDerived>. ");

      return 1. / getRef().getLifetime(particle);
    }
  };

  /**
   * ProcessTraits specialization to flag DecayProcess objects
   **/
  template <typename TProcess>
  struct is_decay_process<
      TProcess,
      std::enable_if_t<std::is_base_of_v<DecayProcess<typename std::decay_t<TProcess>>,
                                         typename std::decay_t<TProcess>>>>
      : std::true_type {};

  /** @} */

} // namespace corsika
