/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <random>

namespace corsika {

  /**
   * Describes a random distribution with \f[ p(x) \sim x^{\alpha} \f] with spectral index
   * \f[ \alpha \f] for a physical quantity \f[x\f] of type Quantity within an interval
   * \f[ [x_0, x_1) \f].
   *
   * @tparam Quantity is the type of the physical quantity.
   */

  template <typename Quantity>
  class PowerLawDistribution {
    static_assert(is_quantity_v<Quantity>, "usable only with Quantity types");

    using real_type = typename Quantity::value_type;
    using uniform_distribution_type = std::uniform_real_distribution<real_type>;

  public:
    using value_type = Quantity;

    PowerLawDistribution(real_type index /** spectral index */,
                         value_type low /** lower bound */,
                         value_type high /** higher bound */)
        : low_{low}
        , high_{high}
        , index_{index} {}

    PowerLawDistribution(PowerLawDistribution<value_type> const& other) = default;
    PowerLawDistribution<value_type>& operator=(
        PowerLawDistribution<value_type> const& other) = default;

    template <class Generator>
    value_type operator()(Generator& g) {
      real_type const u = dist_(g);
      if (index_ == static_cast<real_type>(-1)) {
        return low_ * pow(high_ / low_, u);
      } else {
        return value_type{phys::units::detail::magnitude_tag,
                          pow(pow(low_.magnitude(), index_ + 1) * (1 + u) -
                                  pow(high_.magnitude(), index_ + 1),
                              1 / (index_ + 1))};
      }
    }

  private:
    uniform_distribution_type dist_{};
    value_type low_, high_;
    real_type index_;
  };

} // namespace corsika
