/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/histogram.hpp>

namespace corsika::history {

  namespace detail {
    inline auto hist_factory() {
      namespace bh = boost::histogram;
      namespace bha = bh::axis;
      auto h = bh::make_histogram(
          bha::regular<float, bha::transform::log>{11 * 5, 1e0, 1e11, "muon energy"},
          bha::regular<float, bha::transform::log>{11 * 5, 1e0, 1e11,
                                                   "projectile energy"},
          bha::category<int, bh::use_default, bha::option::growth_t>{
              {211, -211, 2212, -2212}, "projectile PDG"});
      return h;
    }
  } // namespace detail

} // namespace corsika::history
