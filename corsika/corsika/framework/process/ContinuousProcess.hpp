/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/process/ProcessReturn.hpp>
#include <corsika/framework/process/ProcessTraits.hpp>

#include <corsika/detail/framework/process/ContinuousProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
    * @ingroup Processes
    * @{
    * Processes with continuous effects along a particle Trajectory.
    *
    * Create a new ContinuousProcess, e.g. for XYModel, via:
    * @code{.cpp}
    * class XYModel : public ContinuousProcess<XYModel> {};
    * @endcode
    *
    * and provide two necessary interface methods:
    * @code{.cpp}
    * template <typename TParticle, typename TTrack>
    * LengthType getMaxStepLength(TParticle const& p, TTrack const& track) const;
    * @endcode

    * which allows any ContinuousProcess to tell to CORSIKA a maximum
    * allowed step length. Such step-length limitation, if it turns out
    * to be smaller/sooner than any other limit (decay length,
    * interaction length, other continuous processes, geometry, etc.)
    * will lead to a limited step length.

    * @code{.cpp}
    * template <typename TParticle, typename TTrack>
    * ProcessReturn doContinuous(TParticle& p, TTrack const& t, bool const stepLimit)
    * const;
    * @endcode

    * which applied any continuous effects on Particle p along
    * Trajectory t. The particle in all typical scenarios will be
    * altered by a doContinuous. The flag stepLimit will be true if the
    * preious evaluation of getMaxStepLength resulted in this
    * particular ContinuousProcess to be responsible for the step
    * length limit on the current track t. This information can be
    * expoited and avoid e.g. any uncessary calculations.

    * Particle and Track are the valid classes to
    * access particles and track (Trajectory) data on the Stack. Those two methods
    * do not need to be templated, they could use the types
    * e.g. corsika::setup::Stack::particle_type -- but by the cost of
    * loosing all flexibility otherwise provided.

   */

  template <typename TDerived>
  class ContinuousProcess : public BaseProcess<TDerived> {
  public:
  };

  /**
   * ProcessTraits specialization to flag ContinuousProcess objects.
   */
  template <typename TProcess>
  struct is_continuous_process<
      TProcess, std::enable_if_t<
                    std::is_base_of_v<ContinuousProcess<typename std::decay_t<TProcess>>,
                                      typename std::decay_t<TProcess>>>>
      : std::true_type {};

  /** @} */

} // namespace corsika
