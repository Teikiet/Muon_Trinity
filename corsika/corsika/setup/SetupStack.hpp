/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/detail/setup/SetupStack.hpp>

#include <array>
#include <memory>

namespace corsika::setup {

  // ---------------------------------------
  // this is the stack we use in C8 executables:

#ifdef WITH_HISTORY

#include <corsika/stack/history/HistoryStackExtension.hpp>
#include <corsika/stack/history/HistorySecondaryProducer.hpp>

  /*
   * the version with history
   */
  template <typename TEnvironment>
  using Stack = typename detail::StackGenerator<TEnvironment>::StackWithHistory;

#else // WITH_HISTORY

  /*
   * the version without history (and geometry data and weights)
   */
  template <typename TEnvironment>
  using Stack = typename detail::StackGenerator<TEnvironment>::StackWithWeight;

#endif

  // the correct secondary stack view
  template <typename TEnvironment>
  using StackView = typename Stack<TEnvironment>::stack_view_type;

} // namespace corsika::setup
