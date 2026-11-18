/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/InteractionHistogram.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/geometry/FourVector.hpp>

namespace corsika {

  /**
   * @ingroup Processes
   * @{
   *
   * Wrapper around an InteractionProcess that fills histograms of the number
   * of calls to `doInteraction()` binned in projectile energy (both in
   * lab and center-of-mass frame) and species.
   *
   * Use by wrapping a normal InteractionProcess:
   * @code{.cpp}
   * InteractionProcess collision1;
   * InteractionClounter<collision1> counted_collision1;
   * @endcode
   */

  template <class TCountedProcess>
  class InteractionCounter
      : public InteractionProcess<InteractionCounter<TCountedProcess>> {

  public:
    InteractionCounter(TCountedProcess& process);

    /**
     * Wrapper around internal process doInteraction.
     */
    template <typename TSecondaryView>
    void doInteraction(TSecondaryView& view, Code const, Code const, FourMomentum const&,
                       FourMomentum const&);

    /**
     * Wrapper around internal process getCrossSection.
     */
    CrossSectionType getCrossSection(Code const, Code const, FourMomentum const&,
                                     FourMomentum const&) const;

    /**
     * returns the filles histograms.
     *
     * @return InteractionHistogram, which contains the histogram data
     */
    InteractionHistogram const& getHistogram() const;

    /**
     * Clear accumulated histogram data (e.g. at the start of a new shower when using
     * per-shower histograms).
     */
    void resetHistogram();

  private:
    TCountedProcess& process_;
    InteractionHistogram histogram_;
  };

  //! @}

} // namespace corsika

#include <corsika/detail/framework/process/InteractionCounter.inl>
