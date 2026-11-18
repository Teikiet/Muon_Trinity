/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/filesystem/fstream.hpp>

namespace corsika {

  inline void YAMLStreamer::writeYAML(YAML::Node const& node,
                                      boost::filesystem::path const& path) const {

    // construct a YAML emitter for this config file
    YAML::Emitter out;

    // and write the node to the output
    out << node;

    // open the output file - this is <output name>.yaml
    boost::filesystem::ofstream file(path);

    // dump the YAML to the file
    file << out.c_str() << std::endl;
  }

} // namespace corsika
