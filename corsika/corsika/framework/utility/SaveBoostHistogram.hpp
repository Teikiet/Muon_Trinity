/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/histogram.hpp>

namespace corsika {

  /**
   * @file SaveBoostHistogram.hpp
   * @ingroup Utilities
   *
   * This functions saves a boost::histogram into a numpy file. Only rather basic axis
   * types are supported: regular, variable, integer, category<int>. Only "ordinary" bin
   * counts (i.e. a double or int) are supported, nothing fancy like profiles.
   *
   * Note that this function makes a temporary, dense copy of the histogram, which could
   * be an issue for huge sizes (e.g. for high dimensions)
   *
   * @param overwrite silently overwrite existing files if true, otherwise throw
   * runtime_error
   */
  template <class Axes, class Storage>
  void save_hist(boost::histogram::histogram<Axes, Storage> const& h,
                 std::string const& filename, bool overwrite = true);
} // namespace corsika

#include <corsika/detail/framework/utility/SaveBoostHistogram.inl>
