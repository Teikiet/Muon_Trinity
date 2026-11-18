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

  /**
   * Describes a random distribution with \f[ \beta e^{-\beta X} \f] for a physical
   * quantity of type Quantity.
   *
   * @tparam Quantity is the type of the physical quantity.
   */

  template <typename Quantity>
  class ExponentialDistribution {
    static_assert(is_quantity_v<Quantity>, "usable only with Quantity types");

    typedef typename Quantity::value_type real_type;
    typedef std::exponential_distribution<real_type> distribution_type;

  public:
    typedef Quantity value_type;

    ExponentialDistribution() = delete;

    ExponentialDistribution(value_type beta)
        : beta_(beta) {}

    ExponentialDistribution(ExponentialDistribution<value_type> const& other)
        : beta_(other.getBeta()) {}

    ExponentialDistribution<value_type>& operator=(
        ExponentialDistribution<value_type> const& other) {
      if (this == &other) return *this;
      beta_ = other.getBeta();
      return *this;
    }

    /**
     * @fn value_type getBeta() const
     * @brief Get parameter of exponential distribution \f[ \beta e^{-\beta X}\f].
     *
     * @return  value_type
     */
    value_type getBeta() const { return beta_; }

    /**
     * @fn void setBeta(value_type)
     * @brief Set parameter of exponential distribution \f[ \beta e^{-\beta X}\f].
     *
     * @param vBeta
     */
    void setBeta(value_type const& beta) { beta_ = beta; }

    /**
     * @fn value_type operator ()(Generator&)
     * @brief Generate a random number distributed like \f[ \beta e^{-\beta X}\f].
     *
     * @tparam Generator
     * @param g
     * @return
     */
    template <class Generator>
    value_type operator()(Generator& g) {
      return beta_ * dist_(g);
    }

  private:
    distribution_type dist_{1.};
    value_type beta_;
  };

} // namespace corsika
