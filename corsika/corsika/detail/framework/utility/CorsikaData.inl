/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/corsika.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <boost/filesystem/path.hpp>

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace corsika {

  inline boost::filesystem::path corsika_data(boost::filesystem::path const& key) {
    boost::filesystem::path fname =
        boost::filesystem::path(corsika::CORSIKA_DATA_DIR) / key;
    // LCOV_EXCL_START, this cannot be easily tested system-independently
    if (auto const* p = std::getenv("CORSIKA_DATA"); p != nullptr) {
      fname = boost::filesystem::path(p) / key;
    }
    // LCOV_EXCL_STOP
    CORSIKA_LOG_INFO("opening data file={}", fname);
    return fname;
  }
} // namespace corsika
