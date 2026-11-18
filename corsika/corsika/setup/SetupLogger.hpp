/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika::setup {

  /**
   * The default "corsika" logger.
   */
  static std::shared_ptr<spdlog::logger> corsika_logger = get_logger("corsika", true);
  // corsika_logger->set_pattern("[%n:%^%-8l%$] custom pattern: %v");

} // namespace corsika::setup
