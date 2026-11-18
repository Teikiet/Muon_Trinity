/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/functional/identity.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {
  namespace InterpolationTransforms {
    using Identity = boost::identity;

    struct Log {
      template <typename T>
      auto operator()(T val) const {
        if constexpr (is_quantity_v<T>) {
          return std::log(val.magnitude());
        } else {
          return std::log(val);
        }
      }
    };
  } // namespace InterpolationTransforms

  template <typename Transform>
  class CrossSectionTable {
  public:
    template <typename InputIt1, typename InputIt2>
    CrossSectionTable(InputIt1 energiesFirst, InputIt1 energiesLast,
                      InputIt2 crosssectionsFirst) {
      static_assert(std::is_same_v<typename InputIt1::value_type, HEPEnergyType>);
      static_assert(std::is_same_v<typename InputIt2::value_type, CrossSectionType>);

      table_.reserve(std::distance(energiesFirst, energiesLast));

      auto itE = boost::make_transform_iterator(std::move(energiesFirst), transform_);
      decltype(itE) const itEEnd =
          boost::make_transform_iterator(std::move(energiesLast), transform_);

      for (auto itXS = std::move(crosssectionsFirst); itE != itEEnd; ++itE, ++itXS) {
        table_.emplace_back(*itE, *itXS);
      }

      std::sort(table_.begin(), table_.end(), less_datapoint);
    }

    CrossSectionType interpolate(HEPEnergyType energy) const {
      auto const transformed_val = transform_(energy);
      auto const lb_it =
          std::lower_bound(table_.cbegin(), table_.cend(), transformed_val, less_x);
      if (lb_it == table_.cbegin()) {
        throw std::runtime_error{"CrossSectionTable: value out of bounds (lower limit)"};
      }
      if (lb_it == table_.cend()) {
        throw std::runtime_error{"CrossSectionTable: value out of bounds (upper limit)"};
      }

      auto const prev_it = std::prev(lb_it);
      auto const lambda =
          (transformed_val - prev_it->first) / (lb_it->first - prev_it->first);
      return lambda * lb_it->second + (1 - lambda) * prev_it->second;
    }

  private:
    using key_type = std::invoke_result_t<Transform, HEPEnergyType>;
    using datapoint_type = std::pair<key_type, CrossSectionType>;
    std::vector<datapoint_type> table_;
    Transform transform_{};

    static bool less_datapoint(datapoint_type const& a, datapoint_type const& b) {
      return a.first < b.first;
    }
    static bool less_x(datapoint_type const& a, typename datapoint_type::first_type b) {
      return a.first < b;
    }
  };
} // namespace corsika
