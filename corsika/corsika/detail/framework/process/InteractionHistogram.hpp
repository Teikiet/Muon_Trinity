/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/histogram.hpp>

namespace corsika {

  namespace detail {
    inline auto hist_factory(unsigned int const bin_number, double const e_low,
                             double const e_high) {
      namespace bh = boost::histogram;
      namespace bha = bh::axis;

      auto h = bh::make_histogram(
          bha::category<int, bh::use_default, bha::option::growth_t>{{2212, 2112},
                                                                     "projectile PDG"},
          bha::regular<float, bha::transform::log>{bin_number, (float)e_low,
                                                   (float)e_high, "energy/eV"});
      return h;
    }
  } // namespace detail
} // namespace corsika
