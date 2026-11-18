/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <tests/common/TestStack.hpp>

#include <array>
#include <memory>

namespace corsika::test {

  // ---------------------------------------
  // this is the stack we use in C8 executables:

#ifdef WITH_HISTORY

#include <corsika/stack/history/HistoryStackExtension.hpp>
#include <corsika/stack/history/HistorySecondaryProducer.hpp>

  /*
   * the version with history
   */
  using Stack = detail::StackWithHistory;

#else // WITH_HISTORY

  /*
   * the version without history
   */
  using Stack = detail::StackWithWeight;

#endif

  // the correct secondary stack view
  using StackView = typename Stack::stack_view_type;

} // namespace corsika::test
