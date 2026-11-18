/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  /**
   * This is the base class for all outputs so that they
   * can be stored in homogeneous containers.
   */
  class NoOutput : public BaseOutput {

  protected:
    NoOutput() {}

  public:
    /**
     * Called at the start of each run.
     */
    void startOfLibrary(boost::filesystem::path const&) final override {}

    /**
     * Called at the start of each event/shower.
     */
    void startOfShower(unsigned int const) final override {}

    /**
     * Called at the end of each event/shower.
     */
    void endOfShower(unsigned int const) final override {}

    /**
     * Called at the end of each run.
     */
    void endOfLibrary() final override {}

    /**
     * Get the configuration of this output.
     */
    YAML::Node getConfig() const override { return YAML::Node(); };

    /**
     * Get any summary information for the entire library.
     */
    YAML::Node getSummary() const final override { return YAML::Node(); };

  protected:
    void write(Code const&, units::si::HEPEnergyType const&, units::si::LengthType const&,
               units::si::LengthType const&, units::si::TimeType const&) {}
  };

} // namespace corsika
