/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <cnpy.hpp>

#include <boost/histogram.hpp>
#include <boost/filesystem.hpp> // can be changed to std::filesystem if compiler supports it

#include <functional>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>
#include <string>

namespace corsika {

  template <class Axes, class Storage>
  inline void save_hist(boost::histogram::histogram<Axes, Storage> const& h,
                        std::string const& filename, bool overwrite) {
    if (boost::filesystem::exists(filename)) {
      if (overwrite) {
        boost::filesystem::remove(filename);
      } else {
        using namespace std::literals;
        throw std::runtime_error(
            ("save_hist(): "s + filename + " already exists"s).c_str());
      }
    }

    unsigned const rank = h.rank();

    std::vector<size_t> axes_dims;
    axes_dims.reserve(rank);

    auto overflow = std::make_unique<bool[]>(rank);
    auto underflow = std::make_unique<bool[]>(rank);

    std::vector<char> axis_types;
    axis_types.reserve(rank);

    for (unsigned int i = 0; i < rank; ++i) {
      auto const& ax = h.axis(i);
      unsigned const has_underflow = (ax.options() & 0x01) ? 1 : 0;
      unsigned const has_overflow = (ax.options() & 0x02) ? 1 : 0;

      underflow[i] = has_underflow;
      overflow[i] = has_overflow;

      axes_dims.emplace_back(ax.size() + has_underflow + has_overflow);

      if (ax.continuous()) {
        axis_types.push_back('c');
        std::vector<double> ax_edges;
        ax_edges.reserve(ax.size() + 1);

        for (int j = 0; j < ax.size(); ++j) { ax_edges.push_back(ax.bin(j).lower()); }

        ax_edges.push_back(ax.bin(ax.size() - 1).upper());

        cnpy::npz_save(filename, std::string{"binedges_"} + std::to_string(i),
                       ax_edges.data(), {ax_edges.size()}, "a");
      } else {
        axis_types.push_back('d');
        std::vector<int64_t> bins; // we assume that discrete axes have integer bins
        bins.reserve(ax.size());

        for (int j = 0; j < ax.size(); ++j) { bins.push_back(ax.bin(j).lower()); }

        cnpy::npz_save(filename, std::string{"bins_"} + std::to_string(i), bins.data(),
                       {bins.size()}, "a");
      }
    }

    cnpy::npz_save(filename, std::string{"axistypes"}, axis_types.data(), {rank}, "a");
    cnpy::npz_save(filename, std::string{"overflow"}, overflow.get(), {rank}, "a");
    cnpy::npz_save(filename, std::string{"underflow"}, underflow.get(), {rank}, "a");

    auto const prod_axis_size = std::accumulate(axes_dims.cbegin(), axes_dims.cend(),
                                                unsigned{1}, std::multiplies<>());
    auto temp = std::make_unique<float[]>(prod_axis_size);

    assert(prod_axis_size == h.size());

    // reduce multi-dim. to 1-dim, row-major (i.e., last axis index is contiguous in
    // memory) take special care of underflow bins which have -1 as index and thus need
    // to be shifted by +1

    for (auto&& x : indexed(h, boost::histogram::coverage::all)) {
      int p = 0;

      for (unsigned axis_index = 0; axis_index < rank; ++axis_index) {
        int const offset_underflow = (h.axis(axis_index).options() & 0x01) ? 1 : 0;
        auto k = x.index(axis_index) + offset_underflow;
        p = k + p * axes_dims.at(axis_index);
      }

      temp[p] = *x;
    }

    cnpy::npz_save(filename, "data", temp.get(), axes_dims, "a");
    // In Python this array can directly be assigned to a histogram view if that
    // histogram has its axes correspondingly: hist.view(flow=True)[:] = file['data']
  } // namespace corsika

} // namespace corsika
