/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/ProcessTraits.hpp>

#include <type_traits>
#include <cstddef>

//! @file BaseProcess.hpp

namespace corsika {

  /**
   * @ingroup Processes
   * @{
   *
   * Each process in C8 must derive from BaseProcess
   *
   * The structural base type of a process object in a
   * ProcessSequence. Both, the ProcessSequence and all its elements
   * are of type BaseProcess.
   *
   * @todo rename BaseProcess into just Process
   */

  template <typename TDerived>
  struct BaseProcess {
  protected:
    friend TDerived;

    BaseProcess() = default; // protected constructor will allow only
                             // derived classes to be created, not
                             // BaseProcess itself

    /** @name getRef Return reference to underlying type
     *  @{
     */
    TDerived& getRef() { return static_cast<TDerived&>(*this); }
    const TDerived& getRef() const { return static_cast<const TDerived&>(*this); }
    //! @}

  public:
    static bool const is_process_sequence = false;
    static bool const is_switch_process_sequence = false;

    //! Default number of processes is just one, obviously
    static size_t constexpr getNumberOfProcesses() { return 1; }

    //! Base processor type for use in other template classes
    using process_type = TDerived;
  };

  /**
   * is_process traits specialization to indicate inheritance from BaseProcess.
   */
  template <typename TProcess>
  struct is_process<
      TProcess,
      std::enable_if_t<std::is_base_of_v<BaseProcess<typename std::decay_t<TProcess>>,
                                         typename std::decay_t<TProcess>>>>
      : std::true_type {};

  /**
   * count_processes traits specialization to increase process count by one.
   */
  template <typename TProcess, int N>
  struct count_processes<
      TProcess, N,
      typename std::enable_if_t<is_process_v<std::decay_t<TProcess>> &&
                                !std::decay_t<TProcess>::is_process_sequence>> {
    static size_t constexpr count = N + 1;
  };

  //! @}

} // namespace corsika
