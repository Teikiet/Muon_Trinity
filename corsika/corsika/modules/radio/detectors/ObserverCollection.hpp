/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

namespace corsika {

  /**
   * The base interface for radio detectors.
   * At the moment it is a collection of observers with the same implementation.
   */

  template <typename TObserverImpl>
  class ObserverCollection {

    /**
     * The collection of observers used in this simulation.
     */
    std::vector<TObserverImpl> observers_;

  public:
    /**
     * Add an observer to this radio process.
     *
     * @param observer    The observer to add
     */
    void addObserver(TObserverImpl const& observer);

    /**
     * Get the specific observer at that place in the collection
     *
     * @param index in the collection
     */
    TObserverImpl& at(std::size_t const i);

    TObserverImpl const& at(std::size_t const i) const;

    /**
     * Get the number of observerss in the collection
     */
    int size() const;

    /**
     * Get a *non*-const reference to the collection of observers.
     *
     * @returns    An iterable mutable reference to the observers.
     */
    std::vector<TObserverImpl>& getObservers();

    /**
     * Get a const reference to the collection of observers.
     *
     * @returns    An iterable mutable reference to the observers.
     */
    std::vector<TObserverImpl> const& getObservers() const;

    /**
     * Reset all the observer waveforms.
     */
    void reset();
  }; // END: class RadioDetector

} // namespace corsika

#include <corsika/detail/modules/radio/detectors/ObserverCollection.inl>
