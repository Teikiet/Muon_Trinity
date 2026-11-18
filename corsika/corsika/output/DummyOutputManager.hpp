/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

namespace corsika {

  /*!
   * An output manager that does nothing.
   */
  class DummyOutputManager {

  public:
    /**
     * Construct an OutputManager instance with a name in a given directory.
     *
     * @param name    The name of this output collection.
     * @param dir     The directory where the output directory will be stored.
     */
    DummyOutputManager();

    /**
     * Handle graceful closure of the outputs upon destruction.
     */
    ~DummyOutputManager();

    /**
     * Register an existing output to this manager.
     *
     * @param name    The unique name of this output.
     * @param args... These are perfect forwarded to the
     *                constructor of the output.
     */
    template <typename TOutput>
    void add(std::string const& name, TOutput& output);

    /**
     * Called at the start of each library.
     *
     * This iteratively calls startOfLibrary on each registered output.
     */
    void startOfLibrary();

    /**
     * Called at the start of each event/shower.
     * This iteratively calls startOfEvent on each registered output.
     */
    void startOfShower();

    /**
     * Called at the end of each event/shower.
     * This iteratively calls endOfEvent on each registered output.
     */
    void endOfShower();

    /**
     * Called at the end of each library.
     * This iteratively calls endOfLibrary on each registered output.
     */
    void endOfLibrary();

  }; // class DummyOutputManager

} // namespace corsika

#include <corsika/detail/output/DummyOutputManager.inl>
