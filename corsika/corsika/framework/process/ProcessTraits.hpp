/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

/**
 * @file ProcessTraits.hpp
 */

#include <type_traits>
#include <cstddef>

namespace corsika {

  /**
   * A traits marker to identify BaseProcess, thus any type of process.
   */
  template <typename TProcess, typename TEnable = void>
  struct is_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_process_v = is_process<TProcess>::value;

  /**
   * A traits marker to identify ContinuousProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_continuous_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_continuous_process_v = is_continuous_process<TProcess>::value;

  /**
   * A traits marker to identify DecayProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_decay_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_decay_process_v = is_decay_process<TProcess>::value;

  /**
   * A traits marker to identify StackProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_stack_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_stack_process_v = is_stack_process<TProcess>::value;

  /**
   * A traits marker to identify CascadeEquationsProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_cascade_equations_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_cascade_equations_process_v =
      is_cascade_equations_process<TProcess>::value;

  /**
   * A traits marker to identify SecondariesProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_secondaries_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_secondaries_process_v = is_secondaries_process<TProcess>::value;

  /**
   * A traits marker to identify BoundaryProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_boundary_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_boundary_process_v = is_boundary_process<TProcess>::value;

  /**
   * A traits marker to identify InteractionProcess.
   */
  template <typename TProcess, typename Enable = void>
  struct is_interaction_process : std::false_type {};

  template <typename TProcess>
  bool constexpr is_interaction_process_v = is_interaction_process<TProcess>::value;

  /**
   * A traits marker to identify ProcessSequence that contain a StackProcess.
   */
  template <typename TClass>
  struct contains_stack_process : std::false_type {};

  template <typename TClass>
  bool constexpr contains_stack_process_v = contains_stack_process<TClass>::value;

  /**
   * traits class to count any type of Process, general version.
   */
  template <typename TProcess, int N = 0, typename Enable = void>
  struct count_processes {
    static size_t constexpr count = N;
  };

} // namespace corsika
