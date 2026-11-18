/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/detail/framework/process/CascadeEquationsProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
   * @ingroup Processes
   * @{
   *
   * Processes executing cascade-equations calculations.
   *
   * Create a new CascadeEquationsProcess, e.g. for XYModel, via:
   * @code{.cpp}
   * class XYModel : public CascadeEquationsProcess<XYModel> {};
   * @endcode
   *
   * and provide the necessary interface method:
   * @code{.cpp}
   * template <typename TStack>
   * void doCascadeEquations(TStack& stack);
   * @endcode
   *
   * Cascade equation processes may generate new particles on the stack. They also
   * typically generate output.
   */

  template <typename TDerived>
  class CascadeEquationsProcess : public BaseProcess<TDerived> {
  public:
  };

  /**
   * ProcessTraits specialization to flag CascadeEquationsProcess objects.
   */
  template <typename TProcess>
  struct is_cascade_equations_process<
      TProcess, std::enable_if_t<std::is_base_of_v<
                    CascadeEquationsProcess<typename std::decay_t<TProcess>>,
                    typename std::decay_t<TProcess>>>> : std::true_type {};

  //! @}

} // namespace corsika
