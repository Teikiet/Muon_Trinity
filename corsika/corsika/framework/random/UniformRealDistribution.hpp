/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <random>

namespace corsika {

  template <typename Quantity>
  class UniformRealDistribution {
    static_assert(is_quantity_v<Quantity>, "usable only with Quantity types");

    typedef typename Quantity::value_type real_type;
    typedef std::uniform_real_distribution<real_type> distribution_type;

  public:
    typedef Quantity value_type;

    UniformRealDistribution() = delete;

    UniformRealDistribution(value_type b)
        : min_{value_type::zero()}
        , diff_{b} {}

    UniformRealDistribution(value_type pmin, value_type pmax)
        : min_(pmin)
        , diff_{pmax - pmin} {}

    UniformRealDistribution(UniformRealDistribution<value_type> const& other)
        : min_{other.min_}
        , diff_{other.diff_} {}

    UniformRealDistribution<value_type>& operator=(
        UniformRealDistribution<value_type> const& other) {
      if (this == &other) return *this;
      min_ = other.min_;
      diff_ = other.diff_;
      return *this;
    }

    /**
     * Get the upper limit.
     *
     * @return value_type
     */
    value_type getMax() const { return min_ + diff_; }

    /**
     * Set the upper bound.
     *
     * @param pmax: new upper bound
     */
    void setMax(value_type const pmax) { diff_ = pmax - min_; }

    /**
     * Get the lower bound.
     *
     * @return The minimum possible value.
     */
    value_type getMin() const { return min_; }

    /**
     * Set the lower bound.
     *
     * @param pmin is the new lower bound.
     */
    void setMin(value_type const pmin) { min_ = pmin; }

    /**
     * Generate a random number in the range [min, max).
     *
     * @tparam TGenerator
     * @param g
     * @return value_type
     */
    template <class TGenerator>
    value_type operator()(TGenerator& g) {
      return min_ + dist_(g) * diff_;
    }

  private:
    distribution_type dist_{real_type(0.), real_type(1.)};

    value_type min_;
    value_type diff_;
  };

} // namespace corsika
