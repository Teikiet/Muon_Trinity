/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

#include <corsika/detail/framework/process/InteractionProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
   * @ingroup Processes
   * @{
   *
   * Process describing the interaction of particles.
   *
   * Create a new InteractionProcess, e.g. for XYModel, via:
   * @code
   * class XYModel : public InteractionProcess<XYModel> {};
   * @endcode
   *
   * and provide the two necessary interface methods:
   * @code
   * template <typename TSecondaryView>
   * void XYModel::doInteraction(TSecondaryView&);
   *
   * template <typename TParticle>
   * GrammageType XYModel::getInteractionLength(TParticle const&)
   * @endcode
   *
   * Where, of course, SecondaryView and Particle are the valid
   * classes to access particles on the Stack. In user code, those two methods do
   * not need to be templated, they could use the types
   * e.g. corsika::setup::Stack::particle_type -- but by the cost of
   * loosing all flexibility otherwise provided.
   *
   * SecondaryView allows to retrieve the properties of the projectile
   * particles, AND to store new particles (secondaries) which then
   * subsequently can be processes by SecondariesProcess. This is how
   * the output of interactions can be studied right away.
   */

  template <typename TModel>
  class InteractionProcess : public BaseProcess<TModel> {

  public:
    using BaseProcess<TModel>::getRef;
  };

  /**
   * ProcessTraits specialization to flag InteractionProcess objects.
   */
  template <typename TProcess>
  struct is_interaction_process<
      TProcess, std::enable_if_t<
                    std::is_base_of_v<InteractionProcess<typename std::decay_t<TProcess>>,
                                      typename std::decay_t<TProcess>>>>
      : std::true_type {};

  /** @} */

} // namespace corsika
