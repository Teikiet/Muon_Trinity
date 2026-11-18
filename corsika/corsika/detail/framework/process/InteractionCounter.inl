/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/InteractionHistogram.hpp>

namespace corsika {

  template <class TCountedProcess>
  inline InteractionCounter<TCountedProcess>::InteractionCounter(TCountedProcess& process)
      : process_(process) {}

  template <class TCountedProcess>
  template <typename TSecondaryView>
  inline void InteractionCounter<TCountedProcess>::doInteraction(
      TSecondaryView& view, Code const projectileId, Code const targetId,
      FourMomentum const& projectileP4, FourMomentum const& targetP4) {
    size_t const massNumber = is_nucleus(targetId) ? get_nucleus_A(targetId) : 1;
    auto const massTarget = massNumber * constants::nucleonMass;
    histogram_.fill(projectileId, projectileP4.getTimeLikeComponent(), massTarget);
    process_.doInteraction(view, projectileId, targetId, projectileP4, targetP4);
  }

  template <class TCountedProcess>
  inline CrossSectionType InteractionCounter<TCountedProcess>::getCrossSection(
      Code const projectileId, Code const targetId, FourMomentum const& projectileP4,
      FourMomentum const& targetP4) const {
    return process_.getCrossSection(projectileId, targetId, projectileP4, targetP4);
  }

  template <class TCountedProcess>
  inline InteractionHistogram const& InteractionCounter<TCountedProcess>::getHistogram()
      const {
    return histogram_;
  }

  template <class TCountedProcess>
  inline void InteractionCounter<TCountedProcess>::resetHistogram() {
    histogram_.reset();
  }

} // namespace corsika
