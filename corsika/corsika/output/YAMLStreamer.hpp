
/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>
#include <string>

namespace corsika {

  /**
   * This class automates the construction of simple tabular
   * YAML files using the YAML::StreamWriter.
   */
  class YAMLStreamer {
  public:
    inline void writeYAML(YAML::Node const& node,
                          boost::filesystem::path const& path) const;

  }; // class YAMLStreamer
} // namespace corsika

#include <corsika/detail/output/YAMLStreamer.inl>
