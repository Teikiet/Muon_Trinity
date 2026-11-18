/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <boost/filesystem.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/core/PhysicalGeometry.hpp>
#include <corsika/output/ParquetStreamer.hpp>

namespace corsika {

  /**
   * A common abstract interface for radio oberservers.
   *
   * All concrete observer implementations should be of
   * type Observer<T> where T is a concrete observer implementation.
   *
   */
  template <typename TObserverImpl>
  class Observer {

  protected:
    std::string const name_;                   ///< The name/identifier of this observer.
    Point const location_;                     ///< The location of this observer.
    CoordinateSystem const& coordinateSystem_; ///< The coordinate system of the observer

  public:
    using axistype = std::vector<long double>;

    /**
     * \brief Construct a base observer instance.
     *
     * @param name    A name for this observer.
     * @param location    The location of this observer.
     *
     */
    Observer(std::string const& name, Point const& location,
             CoordinateSystem const& coordinateSystem);

    /**
     * Receive a signal at this observer.
     *
     * This is a general implementation call that must be specialized
     * for the particular observer implementation and usage.
     *
     */
    template <typename... TVArgs>
    void receive(TVArgs&&... args);

    /**
     * Get the location of this observer.
     */
    Point const& getLocation() const;

    /**
     * Get the name of this name observer.
     *
     * This is used in producing the output data file.
     */
    std::string const& getName() const;

    /**
     * Reset the observer before starting a new simulation.
     */
    void reset();

    /**
     * Return a reference to the x-axis labels (i.e. time or frequency).
     *
     * This should be an xtensor-convertible type with
     * a ->data() method that converts to a raw pointer.
     */
    axistype getAxis() const;

    /**
     * Return a reference to the underlying waveform data for X polarization.
     *
     * This is used when writing the observer information to disk
     * and will be converted to a 32-bit float before writing.
     */
    std::vector<double> const& getWaveformX() const;

    /**
     * Return a reference to the underlying waveform data for Y polarization.
     *
     * This is used when writing the observer information to disk
     * and will be converted to a 32-bit float before writing.
     */
    std::vector<double> const& getWaveformY() const;

    /**
     * Return a reference to the underlying waveform data for Z polarization.
     *
     * This is used when writing the observer information to disk
     * and will be converted to a 32-bit float before writing.
     */
    std::vector<double> const& getWaveformZ() const;

    /**
     * Get a reference to the underlying radio implementation.
     */
    TObserverImpl& implementation();

  }; // END: class Observer final

} // namespace corsika

#include <corsika/detail/modules/radio/observers/Observer.inl>
