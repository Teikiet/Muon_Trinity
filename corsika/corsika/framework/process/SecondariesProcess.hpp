/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/detail/framework/process/SecondariesProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
   * @ingroup Processes
   * @{
   *
   * Processes acting on the secondaries produced by other processes.
   *
   * Create a new SecondariesProcess, e.g. for XYModel, via
   * @code{.cpp}
   * class XYModel : public SecondariesProcess<XYModel> {};
   * @endcode
   *
   * and provide the necessary interface method:
   * @code{.cpp}
   * template <typename TStackView>
   * void doSecondaries(TStackView& StackView);
   * @endcode
   *
   * where StackView is an object that can store secondaries on a
   * stack and also iterate over these secondaries.
   */

  template <typename TDerived>
  class SecondariesProcess : public BaseProcess<TDerived> {
  public:
  };

  /**
   * ProcessTraits specialization to flag SecondariesProcess objects.
   */
  template <typename TProcess>
  struct is_secondaries_process<
      TProcess, std::enable_if_t<
                    std::is_base_of_v<SecondariesProcess<typename std::decay_t<TProcess>>,
                                      typename std::decay_t<TProcess>>>>
      : std::true_type {};

  //! @}

} // namespace corsika
