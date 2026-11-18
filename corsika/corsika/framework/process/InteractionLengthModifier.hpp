/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <functional>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>

namespace corsika {
  /*!
    @ingroup Processes
    @{

   * Wrapper around an InteractionProcess that allows modifying the interaction
   * length returned by the underlying process in a user-defined way.   *
   */
  template <class TUnderlyingProcess>
  class InteractionLengthModifier
      : public InteractionProcess<InteractionLengthModifier<TUnderlyingProcess>> {

    //! identity function as default modifier
    static auto constexpr non_modifying_functor = [](GrammageType original, corsika::Code,
                                                     HEPEnergyType) { return original; };

  public:
    //! The signature of the modifying functor. Arguments are original int. length, PID,
    //! energy
    using functor_signature = GrammageType(GrammageType, corsika::Code, HEPEnergyType);

    /**
     * Create wrapper around InteractionProcess. Note that the passed process object
     * itself may no longer be used, only through this class.
     */
    InteractionLengthModifier(TUnderlyingProcess&& process,
                              std::function<functor_signature> modifier);

    //! wrapper around internal process doInteraction
    template <typename TSecondaryView>
    void doInteraction(TSecondaryView& view);

    ///! returns underlying process getInteractionLength modified
    template <typename TParticle>
    GrammageType getInteractionLength(TParticle const& particle);

    ///! obtain reference to wrapped process
    TUnderlyingProcess const& getProcess() const;
    TUnderlyingProcess& getProcess();

  private:
    TUnderlyingProcess process_;
    std::function<functor_signature> const modifier_{non_modifying_functor};
  };

  //! @}

} // namespace corsika

#include <corsika/detail/framework/process/InteractionLengthModifier.inl>
