/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/output/Configurable.hpp>
#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>

namespace corsika {

  /**
   * This is the base class for all classes that have
   * YAML representations of their configurations.
   */
  class Configurable {

  protected:
    Configurable() = default;
    virtual ~Configurable() = default;

  public:
    /**
     * Provide YAML configuration for this BaseOutput.
     */
    virtual YAML::Node getConfig() const = 0; // LCOV_EXCL_LINE
  };

} // namespace corsika
