/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/geometry/Plane.h>
#include <corsika/process/ContinuousProcess.h>
#include <corsika/units/PhysicalUnits.h>

#include <boost/histogram.hpp>

#include <functional>

// the detail namespace: here the histrograms are defined
//! \todo add options/parameters to 'detail::hist_factory()'
#include <corsika/detail/stack/history/HistoryObservationPlane.hpp>

namespace corsika::history {

  template <typename TStack>
  class HistoryObservationPlane
      : public ContinuousProcess<HistoryObservationPlane<TStack>> {
  public:
    HistoryObservationPlane(TStack const&, Plane const&, bool = true);

    template <typename TParticle, typename TTrajectory>
    LengthType getMaxStepLength(TParticle const&, TTrajectory const& vTrajectory);

    template <typename TParticle, typename TTrajectory>
    ProcessReturn doContinuous(TParticle const& vParticle,
                               TTrajectory const& vTrajectory);

    auto const& histogram() const { return histogram_; }

  private:
    template <typename TParticle>
    void fillHistoryHistogram(TParticle const&);

    TStack const& stack_;
    Plane const plane_;
    bool const deleteOnHit_;

    decltype(detail::hist_factory()) histogram_ = detail::hist_factory();
  };
} // namespace corsika::history

#include <corsika/detail/stack/history/HistoryObservationPlane.inl>
