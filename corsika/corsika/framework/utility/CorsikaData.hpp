/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/filesystem/path.hpp>

namespace corsika {
  /**
   * @file CorsikaData.hpp
   * @ingroup Utilities
   * @{
   * returns the full path of the file \p filename within the CORSIKA_DATA directory.
   */
  boost::filesystem::path corsika_data(boost::filesystem::path const& filename);

  //! @}
} // namespace corsika

#include <corsika/detail/framework/utility/CorsikaData.inl>
