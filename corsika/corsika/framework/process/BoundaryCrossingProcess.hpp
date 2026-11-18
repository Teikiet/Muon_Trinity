/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/ProcessReturn.hpp>
#include <corsika/media/Environment.hpp>

#include <type_traits>

#include <corsika/detail/framework/process/BoundaryCrossingProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
     @ingroup Processes
     @{

     Processes acting on the particles traversion from one volume into
     another volume.

     Create a new BoundaryCrossingProcess, e.g. for XYModel, via
     @code{.cpp}
     class XYModel : public BoundaryCrossingProcess<XYModel> {};
     @endcode

     and provide the necessary interface method:
     @code{.cpp}
     template <typename TParticle>
     ProcessReturn XYModel::doBoundaryCrossing(TParticle& Particle,
                                     typename TParticle::node_type const& from,
                                     typename TParticle::node_type const& to);
     @endcode

     where Particle is the object to read particle data from a
     Stack. The volume the particle is originating from is `from`, the
     volume where it goes to is `to`.
   */

  template <typename TDerived>
  class BoundaryCrossingProcess : public BaseProcess<TDerived> {
  public:
  };

  /**
   * ProcessTraits specialization to flag BoundaryProcess objects
   **/
  template <typename TProcess>
  struct is_boundary_process<TProcess,
                             std::enable_if_t<std::is_base_of_v<
                                 BoundaryCrossingProcess<typename std::decay_t<TProcess>>,
                                 typename std::decay_t<TProcess>>>> : std::true_type {};

  /** @} */

} // namespace corsika
