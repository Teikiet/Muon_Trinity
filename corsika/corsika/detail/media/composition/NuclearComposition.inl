/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace corsika {
  namespace media {

    inline NuclearComposition::NuclearComposition(std::vector<Code> const& pComponents,
                                                  std::vector<double> const& pFractions)
        : numberFractions_(pFractions)
        , components_(pComponents)
        , avgMassNumber_(getWeightedSum([](Code const compID) -> double {
          if (is_nucleus(compID)) {
            return get_nucleus_A(compID);
          } else {
            return get_mass(compID) / convert_SI_to_HEP(constants::u);
          }
        })) {
      if (pComponents.size() != pFractions.size()) {
        throw std::runtime_error(
            "Cannot construct NuclearComposition from vectors of different sizes.");
      }
      auto const sumFractions =
          std::accumulate(pFractions.cbegin(), pFractions.cend(), 0.);

      if (!(0.999 < sumFractions && sumFractions < 1.001)) {
        throw std::runtime_error("element fractions do not add up to 1");
      }
      this->updateHash();
    }

    template <typename TFunction>
    inline auto NuclearComposition::getWeighted(TFunction func) const {
      using ResultQuantity = decltype(func(std::declval<Code>()));
      auto const product = [&](auto const compID, auto const fraction) {
        return func(compID) * fraction;
      };

      if constexpr (phys::units::is_quantity_v<ResultQuantity>) {
        std::vector<ResultQuantity> result(components_.size(), ResultQuantity::zero());
        std::transform(components_.cbegin(), components_.cend(),
                       numberFractions_.cbegin(), result.begin(), product);
        return result;
      } else {
        std::vector<ResultQuantity> result(components_.size(), ResultQuantity(0));
        std::transform(components_.cbegin(), components_.cend(),
                       numberFractions_.cbegin(), result.begin(), product);
        return result;
      }
    } // namespace corsika

    template <typename TFunction>
    inline auto NuclearComposition::getWeightedSum(TFunction func) const
        -> decltype(func(std::declval<Code>())) {
      using ResultQuantity = decltype(func(std::declval<Code>()));

      auto const prod = [&](auto const compID, auto const fraction) {
        return func(compID) * fraction;
      };

      if constexpr (phys::units::is_quantity_v<ResultQuantity>) {
        return std::inner_product(
            components_.cbegin(), components_.cend(), numberFractions_.cbegin(),
            ResultQuantity::zero(), // .zero() is defined for quantity types only
            std::plus<ResultQuantity>(), prod);
      } else {
        return std::inner_product(
            components_.cbegin(), components_.cend(), numberFractions_.cbegin(),
            ResultQuantity(0), // in other cases we have to use a bare 0
            std::plus<ResultQuantity>(), prod);
      }
    }

    inline size_t NuclearComposition::getSize() const { return numberFractions_.size(); }

    inline std::vector<double> const& NuclearComposition::getFractions() const {
      return numberFractions_;
    }

    inline std::vector<Code> const& NuclearComposition::getComponents() const {
      return components_;
    }

    inline double const NuclearComposition::getAverageMassNumber() const {
      return avgMassNumber_;
    }

    template <class TRNG>
    inline Code NuclearComposition::sampleTarget(
        std::vector<CrossSectionType> const& sigma, TRNG&& randomStream) const {
      if (sigma.size() != numberFractions_.size()) {
        throw std::runtime_error("incompatible vector sigma as input");
      }

      auto zip_beg = boost::make_zip_iterator(
          boost::make_tuple(numberFractions_.cbegin(), sigma.cbegin()));
      auto zip_end = boost::make_zip_iterator(
          boost::make_tuple(numberFractions_.cend(), sigma.cend()));
      using zip_iter_type = decltype(zip_beg);

      auto const mult_func = [](zip_iter_type::value_type const& zipit) -> double {
        return zipit.get<0>() * zipit.get<1>().magnitude();
      };

      using transform_iter_type =
          boost::transform_iterator<decltype(mult_func), zip_iter_type, double, double>;

      auto trans_beg = transform_iter_type{zip_beg, mult_func};
      auto trans_end = transform_iter_type{zip_end, mult_func};

      std::discrete_distribution channelDist{trans_beg, trans_end};

      auto const iChannel = channelDist(randomStream);
      return components_[iChannel];
    }

    // Note: when this class ever modifies its internal data, the hash
    // must be updated, too!
    // the hash value is important to find tables, etc.
    inline size_t NuclearComposition::getHash() const { return hash_; }

    inline bool NuclearComposition::operator==(media::NuclearComposition const& v) const {
      return v.hash_ == hash_;
    }

    inline void NuclearComposition::updateHash() {
      std::vector<std::size_t> hashes;
      for (double ifrac : this->getFractions())
        hashes.push_back(std::hash<double>{}(ifrac));
      for (Code icode : this->getComponents())
        hashes.push_back(std::hash<int>{}(static_cast<int>(icode)));
      std::size_t h = std::hash<double>{}(this->getAverageMassNumber());
      for (std::size_t ih : hashes) h = h ^ (ih << 1);
      hash_ = h;
    }

  } // namespace media
} // namespace corsika
